#!/usr/bin/env python3
"""Independently read VTU/PVD XML and compare every numeric array with native CSV."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np


def verify(output):
    output=Path(output).resolve();manifest=json.loads((output/"export_manifest.json").read_text())
    case=(output/manifest["source_case_relative"]).resolve();records=[]
    pvd=ET.parse(output/manifest["pvd"]).findall("./Collection/DataSet")
    assert len(pvd)==len(manifest["frames"])
    max_center=max_width=0.
    for entry,frame in zip(pvd,manifest["frames"]):
        assert float(entry.attrib["timestep"])==frame["time_days"]
        assert entry.attrib["file"]==frame["vtu"]
        source=case/frame["source_csv"]
        assert hashlib.sha256(source.read_bytes()).hexdigest()==frame["source_sha256"]
        with source.open(newline="") as handle:
            reader=csv.DictReader(handle);names=reader.fieldnames;rows=list(reader)
        expected={name:np.array([float(row[name]) for row in rows]) for name in names}
        piece=ET.parse(output/frame["vtu"]).find("./UnstructuredGrid/Piece")
        arrays={item.attrib["Name"]:np.fromstring(item.text,sep=" ") for item in piece.findall("./CellData/DataArray")}
        for name,value in expected.items():np.testing.assert_array_equal(arrays[name],value,err_msg=name)
        np.testing.assert_array_equal(arrays["elevation_m"],-expected["z_depth_m"])
        np.testing.assert_array_equal(arrays["formation_id"],expected.get("layer_id",np.zeros(len(rows))))
        points=np.fromstring(piece.find("./Points/DataArray").text,sep=" ").reshape(-1,3)
        cell={item.attrib["Name"]:np.fromstring(item.text,sep=" ",dtype=int) for item in piece.findall("./Cells/DataArray")}
        np.testing.assert_array_equal(cell["types"],np.full(len(rows),12))
        np.testing.assert_array_equal(cell["offsets"],8*np.arange(1,len(rows)+1))
        indices=cell["connectivity"].reshape(-1,8)
        assert indices.min()>=0 and indices.max()<len(points)
        vertices=points[indices];centers=vertices.mean(axis=1);widths=np.ptp(vertices,axis=1)
        expected_centers=np.column_stack([expected["x_m"],expected["y_m"],-expected["z_depth_m"]])
        expected_widths=np.column_stack([expected["dx_m"],expected["dy_m"],expected["dz_m"]])
        np.testing.assert_allclose(centers,expected_centers,rtol=0,atol=1e-8)
        np.testing.assert_allclose(widths,expected_widths,rtol=0,atol=1e-8)
        # The three edges from node 0 must have positive cell orientation.
        signed=np.einsum('ij,ij->i',np.cross(vertices[:,1]-vertices[:,0],vertices[:,3]-vertices[:,0]),vertices[:,4]-vertices[:,0])
        assert np.all(signed>0)
        max_center=max(max_center,float(np.max(np.abs(centers-expected_centers))))
        max_width=max(max_width,float(np.max(np.abs(widths-expected_widths))))
        records.append({"time_days":frame["time_days"],"cells":len(rows),"arrays_checked":len(expected),"max_scalar_error":0.})
    report={"status":"passed","pilot_only":manifest["pilot_only"],"frames":records,
            "max_cell_center_error_m":max_center,"max_cell_width_error_m":max_width,
            "all_native_numeric_arrays_equal":True,"positive_hexahedron_orientation":True,
            "scope":"Export consistency only; native physical model is validated separately."}
    (output/"export_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1]/"results/paraview")
    args=parser.parse_args();print(json.dumps(verify(args.output),indent=2))
