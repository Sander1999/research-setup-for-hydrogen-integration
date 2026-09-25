#!/usr/bin/env pvpython
"""Build real ParaView 3-D PNGs and a reusable state from verified native VTU/PVD."""
from pathlib import Path
import argparse,json
from paraview.simple import *


def annotation(view,name,text,location="Upper Left Corner",size=15):
    source=Text(registrationName=name);source.Text=text
    display=Show(source,view);display.WindowLocation=location
    display.Color=[.08,.12,.18];display.FontSize=size
    return source


def color(display,view,field,title,bounds,*,log=False,categorical=None):
    ColorBy(display,("CELLS",field));lut=GetColorTransferFunction(field)
    lo,hi=map(float,bounds)
    if hi<=lo:hi=lo+max(abs(lo)*.01,1e-12)
    colors=[(.267,.005,.329),(.230,.322,.546),(.128,.567,.551),(.369,.789,.383),(.993,.906,.144)]
    if categorical:
        lut.InterpretValuesAsCategories=1
        lut.Annotations=[v for pair in sorted(categorical.items(),key=lambda kv:int(kv[0])) for v in pair]
        palette=[(.79,.70,.53),(.35,.59,.75),(.57,.66,.38),(.74,.47,.35),(.58,.52,.71),(.75,.75,.75),(.33,.71,.68)]
        lut.IndexedColors=[v for i in range(len(categorical)) for v in palette[i%len(palette)]]
    else:
        positions=[lo*(hi/lo)**(i/4) if log else lo+(hi-lo)*i/4 for i in range(5)]
        lut.RGBPoints=[v for p,rgb in zip(positions,colors) for v in (p,*rgb)]
        lut.ColorSpace="RGB";lut.RescaleTransferFunction(lo,hi);lut.UseLogScale=int(log)
    lut.AutomaticRescaleRangeMode="Never";display.SetScalarBarVisibility(view,True)
    bar=GetScalarBar(lut,view);bar.Title=title;bar.ComponentTitle=""
    bar.Orientation="Horizontal";bar.WindowLocation="Any Location";bar.Position=[.17,.04]
    bar.ScalarBarLength=.68;bar.TitleColor=[.08,.12,.18];bar.LabelColor=[.08,.12,.18]
    bar.TitleFontSize=12;bar.LabelFontSize=10
    return lut


def configure(view,bounds,last):
    view.UseColorPaletteForBackground=0;view.Background=[.975,.980,.988]
    view.OrientationAxesVisibility=1;view.OrientationAxesLabelColor=[.18,.23,.29]
    view.CameraParallelProjection=1;view.ViewTime=last
    x0,x1,y0,y1,z0,z1=bounds;cx=(x0+x1)/2;cy=(y0+y1)/2;cz=(z0+z1)/2
    sx,sy=x1-x0,y1-y0;extent=max(sx,sy,z1-z0)
    view.CameraFocalPoint=[cx,cy,cz];view.CameraPosition=[cx+1.3*extent,cy-1.6*extent,cz+1.2*extent]
    view.CameraViewUp=[0,0,1];view.CameraParallelScale=.72*extent
    view.AxesGrid.Visibility=1;view.AxesGrid.XTitle="";view.AxesGrid.YTitle="";view.AxesGrid.ZTitle=""
    if "LabelUniqueEdgesOnly" in view.AxesGrid.ListProperties():view.AxesGrid.LabelUniqueEdgesOnly=1
    # Interior coordinate ticks avoid two axis endpoints printing on top of
    # each other at the projected box corners. Titles use a separate note.
    for axis,lo,hi in (("X",x0,x1),("Y",y0,y1),("Z",z0,z1)):
        setattr(view.AxesGrid,axis+"AxisUseCustomLabels",1)
        setattr(view.AxesGrid,axis+"AxisLabels",[lo+(hi-lo)*i/4 for i in (1,2,3)])
    for axis in "XYZ":
        setattr(view.AxesGrid,axis+"TitleColor",[.18,.23,.29])
        setattr(view.AxesGrid,axis+"LabelColor",[.18,.23,.29])


def threshold(reader,name,field,low,high):
    source=Threshold(registrationName=name,Input=reader)
    source.Scalars=["CELLS",field];source.ThresholdMethod="Between"
    source.LowerThreshold=low;source.UpperThreshold=high
    return source


def verify_state(output):
    manifest=json.loads((output/"export_manifest.json").read_text())
    LoadState(str(output/"hydrogen_research.pvsm"))
    reader=next(source for (name,_),source in GetSources().items() if name.startswith("Native accepted"))
    frames=[];scene=GetAnimationScene()
    for frame in (manifest["frames"][0],manifest["frames"][-1]):
        day=frame["time_days"];scene.AnimationTime=day
        for view in GetViews():
            if hasattr(view,"ViewTime"):view.ViewTime=day
        reader.UpdatePipeline(day)
        assert reader.GetDataInformation().GetNumberOfCells()==frame["cells"]
        assert "source_fraction" in reader.CellData.keys()
        RenderAllViews()
        frames.append({"day":day,"cells":frame["cells"],"source_fraction_range":list(reader.CellData["source_fraction"].GetRange())})
    report={"status":"passed","fresh_process_state_reload":True,"first_and_final_frames":frames,
            "pilot_only":manifest["pilot_only"]}
    (output/"state_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report));return report


def build(output):
    output=output.resolve();manifest=json.loads((output/"export_manifest.json").read_text())
    checks=json.loads((output/"export_verification.json").read_text())
    if checks.get("status")!="passed":raise ValueError("Run verify_export.py before rendering")
    reader=PVDReader(registrationName="Native accepted DARTS snapshots | cell scalars",FileName=str(output/"source_model.pvd"))
    scene=GetAnimationScene();scene.UpdateAnimationUsingDataTimeSteps()
    ranges={};reader_checks=[]
    fields=("source_fraction","formation_id","H2_kg_m3_pore_fluid","gas_saturation","damage","permeability_used_md")
    for frame in manifest["frames"]:
        reader.UpdatePipeline(frame["time_days"])
        assert reader.GetDataInformation().GetNumberOfCells()==frame["cells"]
        for field in fields:
            if field not in reader.CellData.keys():raise ValueError(f"Native CSV lacks renderer field {field}")
            low,high=reader.CellData[field].GetRange()
            before=ranges.get(field,[low,high]);ranges[field]=[min(low,before[0]),max(high,before[1])]
        reader_checks.append({"day":frame["time_days"],"cells":frame["cells"]})
    last=manifest["frames"][-1]["time_days"];reader.UpdatePipeline(last);scene.AnimationTime=last
    bounds=list(reader.GetDataInformation().GetBounds());cx=(bounds[0]+bounds[1])/2;cy=(bounds[2]+bounds[3])/2;cz=(bounds[4]+bounds[5])/2
    source=threshold(reader,"Dunite hypothesis | actual fractional source cells","source_fraction",.001,1.)
    source.UpdatePipeline(last)
    outline=Outline(registrationName="Native domain outline",Input=reader)
    clip=Clip(registrationName="Interior half-domain | native cells",Input=reader)
    clip.ClipType="Plane";clip.ClipType.Origin=[cx,cy,cz];clip.ClipType.Normal=[0,1,0];clip.Invert=0
    clip.UpdatePipeline(last)
    layouts=[];views=[]
    for name in ("Geology and source","Hydrogen distribution","Damage and native permeability"):
        layout=CreateLayout(name=name);view=CreateView("RenderView");AssignViewToLayout(view,layout,0)
        layout.SetSize(1600,1000);configure(view,bounds,last);layouts.append(layout);views.append(view)
    right=CreateView("RenderView");layouts[2].SplitHorizontal(0,.5);AssignViewToLayout(right,layouts[2],2)
    layouts[2].SetSize(2000,1050);configure(right,bounds,last);views.append(right)
    geology,hydrogen,damage,perm=views
    for view in (damage,perm):view.CameraParallelScale*=1.12
    for view in views:
        box=Show(outline,view);box.DiffuseColor=[.34,.39,.46];box.LineWidth=1.1
        stamp=AnnotateTimeFilter(registrationName="Native accepted time",Input=reader)
        stamp.Format="Day {time:.0f}";display=Show(stamp,view);display.WindowLocation="Upper Right Corner"
        display.Color=[.08,.12,.18];display.FontSize=12
        if manifest["pilot_only"]:annotation(view,"PILOT marker","PILOT API CHECK — superseded geometry","Lower Left Corner",12)
        else:annotation(view,"Coordinates","x/y in metres; z = minus depth (m)","Lower Left Corner",10)

    host=Show(clip,geology);host.Representation="Surface";host.Opacity=.3
    formation_labels=manifest.get("formation_labels") or {"0":"Unlayered pilot host"}
    color(host,geology,"formation_id","Model formations",ranges["formation_id"],categorical=formation_labels)
    src=Show(source,geology);src.Representation="Surface With Edges";src.DiffuseColor=[.15,.37,.37]
    src.AmbientColor=[.15,.37,.37];src.EdgeColor=[.1,.18,.21];src.LineWidth=.4;ColorBy(src,None)
    annotation(geology,"Geology title","LAYERED HOST + HYPOTHETICAL DUNITE\nActual source-fraction cells; depths retained",size=16)

    conc_max=max(ranges["H2_kg_m3_pore_fluid"][1],1e-12)
    visible_low=max(1e-12,conc_max*.01)
    plume=threshold(reader,"H2-bearing cells | above 1 percent of series peak","H2_kg_m3_pore_fluid",visible_low,conc_max*(1+1e-10))
    plume.UpdatePipeline(last);plume_display=Show(plume,hydrogen);plume_display.Representation="Surface"
    color(plume_display,hydrogen,"H2_kg_m3_pore_fluid","Total H2 (kg / m3 pore fluid)",[0,conc_max])
    annotation(hydrogen,"Hydrogen title",f"NATIVE HYDROGEN DISTRIBUTION\nFree + dissolved H2; no surface arrival inferred\nDisplay > {visible_low:.3g} kg/m3 (1% of series peak)",size=15)

    d=Show(source,damage);d.Representation="Surface With Edges";d.EdgeColor=[.25,.25,.25];d.LineWidth=.3
    color(d,damage,"damage","Empirical damage proxy",[0,1])
    annotation(damage,"Damage title","SOURCE-ROCK DAMAGE\nCell proxy, not resolved fracture geometry",size=14)
    kd=Show(clip,perm);kd.Representation="Surface With Edges";kd.EdgeColor=[.25,.25,.25];kd.LineWidth=.3
    if ranges["permeability_used_md"][0]<=0:raise ValueError("Log permeability display requires positive native permeability")
    color(kd,perm,"permeability_used_md","Permeability used in accepted step (mD)",ranges["permeability_used_md"],log=True)
    annotation(perm,"Permeability title","NATIVE FLOW PERMEABILITY\nInterior half-domain; logarithmic colour scale",size=14)
    SetActiveSource(reader);SetActiveView(geology);RenderAllViews()
    paths={"state":str(output/"hydrogen_research.pvsm")}
    for layout,filename in zip(layouts,("geology_source.png","hydrogen_distribution.png","damage_permeability.png")):
        resolution=[2000,1050] if filename.startswith("damage") else [1600,1000]
        SaveScreenshot(str(output/filename),layout,ImageResolution=resolution);paths[filename]=str(output/filename)
    SetActiveView(geology);SaveState(paths["state"])
    report={"status":"rendered","runtime":"Standalone ParaView pvpython; no DARTS or project package imports",
            "pilot_only":manifest["pilot_only"],"reader_frames":reader_checks,"bounds_elevation_m":bounds,
            "global_colour_ranges":ranges,"hydrogen_display_threshold_kg_m3":visible_low,
            "vertical_exaggeration":1,"paths":paths}
    (output/"render_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report));return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1]/"results/paraview")
    parser.add_argument("--verify-state",action="store_true")
    args=parser.parse_args()
    verify_state(args.output) if args.verify_state else build(args.output)
