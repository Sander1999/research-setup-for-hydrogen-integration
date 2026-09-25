#!/usr/bin/env python3
"""Export completed native cell CSV snapshots to standard VTU/PVD, no simulation."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path
from xml.etree import ElementTree as ET

PROJECT=Path(__file__).resolve().parents[1]
REQUIRED={"x_m","y_m","z_depth_m","dx_m","dy_m","dz_m","source_fraction"}


def read_cells(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        reader=csv.DictReader(handle);fields=reader.fieldnames or [];rows=list(reader)
    if REQUIRED-set(fields):raise ValueError(f"Missing geometry/source columns: {sorted(REQUIRED-set(fields))}")
    if not rows:raise ValueError("Empty native spatial snapshot")
    values={name:[float(row[name]) for row in rows] for name in fields}
    if not all(math.isfinite(value) for vals in values.values() for value in vals):raise ValueError("Native cell arrays must be finite numeric values")
    if any(v<=0 for k in ("dx_m","dy_m","dz_m") for v in values[k]):raise ValueError("Cell widths must be positive")
    if any(not 0<=v<=1 for v in values["source_fraction"]):raise ValueError("Source fractions outside [0,1]")
    return values


def _data(parent,name,values,kind="Float64",components=None):
    attrs={"type":kind,"format":"ascii"}
    if name is not None:attrs["Name"]=name
    if components is not None:attrs["NumberOfComponents"]=str(components)
    child=ET.SubElement(parent,"DataArray",attrs)
    child.text="\n"+" ".join(str(int(v)) if kind in ("Int64","UInt8") else format(v,".17g") for v in values)+"\n"


def write_vtu(path,values):
    n=len(values["x_m"]);points=[];lookup={};connectivity=[]
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for i in range(n):
        bounds=((values["x_m"][i]-values["dx_m"][i]/2,values["x_m"][i]+values["dx_m"][i]/2),
                (values["y_m"][i]-values["dy_m"][i]/2,values["y_m"][i]+values["dy_m"][i]/2),
                (-values["z_depth_m"][i]-values["dz_m"][i]/2,-values["z_depth_m"][i]+values["dz_m"][i]/2))
        for corner in corners:
            # Shared vertices avoid duplicate interior faces. Rounding only
            # geometry to 1e-10 m is far below the recorded 1e-8 m QA tolerance.
            point=tuple(round(bounds[d][corner[d]],10) for d in range(3))
            if point not in lookup:lookup[point]=len(points);points.append(point)
            connectivity.append(lookup[point])
    root=ET.Element("VTKFile",type="UnstructuredGrid",version="0.1",byte_order="LittleEndian")
    grid=ET.SubElement(root,"UnstructuredGrid");piece=ET.SubElement(grid,"Piece",NumberOfPoints=str(len(points)),NumberOfCells=str(n))
    pdata=ET.SubElement(piece,"Points");_data(pdata,None,[v for p in points for v in p],components=3)
    cells=ET.SubElement(piece,"Cells");_data(cells,"connectivity",connectivity,"Int64")
    _data(cells,"offsets",[8*(i+1) for i in range(n)],"Int64");_data(cells,"types",[12]*n,"UInt8")
    cell_data=ET.SubElement(piece,"CellData",Scalars="source_fraction")
    for key,vals in values.items():_data(cell_data,key,vals,"Int64" if key in {"cell","layer_id"} else "Float64")
    _data(cell_data,"elevation_m",[-v for v in values["z_depth_m"]])
    _data(cell_data,"formation_id",values.get("layer_id",[0]*n),"Int64")
    ET.indent(root);ET.ElementTree(root).write(path,encoding="utf-8",xml_declaration=True)
    return n,len(points)


def export_case(source_case,output,*,pilot=False):
    source_case=Path(source_case).resolve();output=Path(output).resolve()
    status=json.loads((source_case/"status.json").read_text())
    if status.get("status")!="complete":raise ValueError("Only completed native cases may be exported")
    summary=json.loads((source_case/"summary.json").read_text())
    parameters=json.loads((source_case/"parameters.json").read_text())
    final_day=float(summary["days"])
    snapshots={0.:source_case/"spatial_initial.csv"}
    for path in source_case.glob("spatial_day_*.csv"):
        match=re.fullmatch(r"spatial_day_([\d.eE+-]+)\.csv",path.name)
        if match:snapshots[float(match.group(1))]=path
    snapshots[final_day]=source_case/"spatial_final.csv"
    if min(snapshots)<0 or max(snapshots)>final_day+1e-8:raise ValueError("Snapshot times outside completed run")
    if output==source_case or output in source_case.parents:raise ValueError("Export must not replace native source data")
    vtkdir=output/"vtk";vtkdir.mkdir(parents=True,exist_ok=True)
    collection=ET.Element("VTKFile",type="Collection",version="0.1",byte_order="LittleEndian")
    coll=ET.SubElement(collection,"Collection");records=[];reference=None
    for index,(day,path) in enumerate(sorted(snapshots.items())):
        values=read_cells(path)
        geom={key:values[key] for key in REQUIRED if key!="source_fraction"}
        if reference is None:reference=geom
        elif geom!=reference:raise ValueError("Snapshot geometry changed; unsupported moving mesh")
        filename=f"snapshot_{index:03d}.vtu";ncells,npoints=write_vtu(vtkdir/filename,values)
        ET.SubElement(coll,"DataSet",timestep=format(day,".17g"),group="",part="0",file="vtk/"+filename)
        records.append({"time_days":day,"source_csv":path.name,"source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                        "vtu":"vtk/"+filename,"cells":ncells,"points":npoints,
                        "arrays":list(values)+["elevation_m","formation_id"]})
    ET.indent(collection);ET.ElementTree(collection).write(output/"source_model.pvd",encoding="utf-8",xml_declaration=True)
    result={"status":"exported","pilot_only":bool(pilot),"source_case_relative":os.path.relpath(source_case,output),
            "source_configuration_sha256":summary.get("configuration_sha256"),"frames":records,
            "pvd":"source_model.pvd","native_summary_days":final_day,
            "formation_labels":{str(i):layer.get("name",f"Formation {i}") for i,layer in enumerate(parameters.get("layers",[]))},
            "coordinate_convention":"Native depth positive down; VTK z=elevation=-depth. x/y and all cell arrays retained in native metres/order.",
            "array_policy":"Every numeric CSV column copied as cell data; elevation_m and formation_id are explicit aliases only.",
            "units":{"pressure_bar":"bar absolute","source_fraction":"cell volume fraction","damage":"empirical scalar damage proxy; not a resolved fracture",
                     "gas_saturation":"pore-fluid gas saturation","H2_kg_m3_pore_fluid":"total H2 kg per m3 pore fluid",
                     "permeability_used_md":"mD used in accepted native step","permeability_next_step_md":"mD after accepted damage; for next step",
                     "formation_id":"integer alias of native layer_id","conversion_fraction_accessible_rock":"fraction of accessible source converted",
                     "bulk_expansion_strain_diagnostic":"diagnostic strain; no simulated mechanical displacement"}}
    (output/"export_manifest.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",type=Path,default=PROJECT/"results/darts/natural_flow")
    parser.add_argument("--output",type=Path,default=PROJECT/"results/paraview")
    parser.add_argument("--pilot",action="store_true",help="Label an API-check export as pilot-only")
    args=parser.parse_args();result=export_case(args.case,args.output,pilot=args.pilot)
    print(json.dumps({"output":str(args.output.resolve()),"frames":len(result["frames"]),"cells":result["frames"][0]["cells"],"pilot_only":result["pilot_only"]}))
