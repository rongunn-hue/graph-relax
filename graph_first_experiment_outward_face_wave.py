#!/usr/bin/env python3
"""Alternating, dynamically adjacent face-by-face vertical expansion."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_corridor_separation as corridor
import graph_first_experiment_outward_boundary_wave as wave
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

TOP_START="face:007";BOTTOM_START="face:008"

def faces_now(vertices,edges,coordinates):
    fs,d=inventory.enumerate_faces(vertices,edges,coordinates)
    inventory.add_classification_and_adjacency(fs,edges,d)
    return {f["face_id"]:f for f in fs}

def mean_y(face,co):return sum(co[v][1] for v in face["distinct_boundary_vertices"])/len(face["distinct_boundary_vertices"])

def next_face(fid,side,faces,co,ycenter):
    cy=mean_y(faces[fid],co);rows=[]
    for n in faces[fid]["neighboring_faces"]:
        f=faces[n["face_id"]]
        if not f["is_bounded"]:continue
        ny=mean_y(f,co)
        if (side=="TOP" and ny>cy+geometry.GEOMETRY_TOLERANCE) or (side=="BOTTOM" and ny<cy-geometry.GEOMETRY_TOLERANCE):
            rows.append((abs(ny-ycenter),f["face_id"],ny))
    rows.sort()
    return (rows[0][1] if rows else None),[{"face_id":f,"mean_y":y} for _,f,y in rows]

def face_gap(a,b,co):
    ea={tuple(x) for x in a["distinct_boundary_edges"]};eb={tuple(x) for x in b["distinct_boundary_edges"]};shared=ea&eb;samples=[]
    for xedge in sorted(ea-shared):
        for yedge in sorted(eb-shared):
            p,q=co[xedge[0]],co[xedge[1]];r,s=co[yedge[0]],co[yedge[1]]
            lo=max(min(p[0],q[0]),min(r[0],s[0]));hi=min(max(p[0],q[0]),max(r[0],s[0]))
            if hi-lo<=geometry.GEOMETRY_TOLERANCE:continue
            for x in sorted({lo,hi}):
                ya=wave.y_values(list(xedge),x,co);yb=wave.y_values(list(yedge),x,co)
                for u in ya:
                    for v in yb:
                        if abs(u-v)>geometry.GEOMETRY_TOLERANCE:samples.append((abs(u-v),x,xedge,yedge))
    if not samples:return None
    gap,x,e1,e2=min(samples)
    return {"full_gap":gap,"half_gap":gap/2,"limiting_x":x,"current_edge":list(e1),"next_edge":list(e2),"shared_edges":[list(x) for x in sorted(shared)]}

def render(vertices,edges,stages,out):
    pts=[p for s in stages for p in s["coordinates"].values()];xs=[p[0] for p in pts];ys=[p[1] for p in pts];m=80;mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys);pw=mxx-mnx+2*m;h=mxy-mny+2*m
    z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{mnx-m:.17g} {mny-m:.17g} {len(stages)*pw:.17g} {h:.17g}" width="{800*len(stages)}" height="800">',f'<rect x="{mnx-m:.17g}" y="{mny-m:.17g}" width="{len(stages)*pw:.17g}" height="{h:.17g}" fill="white"/>']
    for i,st in enumerate(stages):
        co={v:tuple(p) for v,p in st["coordinates"].items()};z.append(f'<g id="stage-{i}" transform="translate({i*pw:.17g} 0)">')
        for a,b in edges:z.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        moved=set(st.get("moved_vertices",[]));color="#00897b" if st.get("side")=="TOP" else "#8e24aa"
        for v in sorted(vertices):
            x,y=co[v];fill=color if v in moved else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            z.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in moved else 6}" fill="{fill}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-size="10">{vertices[v]["name"]}</text>')
        title="ORIGINAL" if i==0 else f'{st["side"]} {st["face_id"]}; moved={",".join(st["moved_vertices"])}; d={st["movement"]:.6g}'
        z.append(f'<text x="{mnx-m+10:.17g}" y="{mny-m+22:.17g}" font-size="15" font-weight="bold">{title}</text><text x="{mnx-m+10:.17g}" y="{mny-m+42:.17g}" font-size="12">crossings={st["metrics"]["proper_crossings"]}; violations={st["metrics"]["clearance_violations"]}</text></g>')
    z.append("</svg>");out.write_text("\n".join(z)+"\n")

def run(graph_path,base_path,outdir):
    base=json.loads(base_path.read_text());g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];co={v:tuple(p) for v,p in base["coordinates_before"].items()}
    ok,d,a,e,c,m,re=wave.valid(vertices,edges,co)
    if not ok or (len(vertices),len(edges),m["proper_crossings"])!=(34,44,0):raise AssertionError("bad baseline")
    ids=set(base["corridor"]["events_included"]);events=[x for x in base["before_crowding_events"] if x["violation_id"] in ids];yc=sum(x["crowding_point"][1] for x in events)/len(events)
    fs=faces_now(vertices,edges,co);states={f:"UNMOVED" for f,x in fs.items() if x["is_bounded"]};vstate={v:"UNMOVED" for v in vertices};top,bottom=TOP_START,BOTTOM_START;active={"TOP":True,"BOTTOM":True};trace=[];considered=[];stages=[{"stage":"ORIGINAL","coordinates":{v:list(p) for v,p in sorted(co.items())},"metrics":m}]
    for guard in range(2*len(states)+1):
        progressed=False
        for side in ("TOP","BOTTOM"):
            if not active[side]:continue
            if active["TOP"] and active["BOTTOM"] and top==bottom:
                trace.append({"side":"BOTH","face_id":top,"status":"STOPPED_WAVES_MET"});active={"TOP":False,"BOTTOM":False};break
            fid=top if side=="TOP" else bottom;opposite="MOVED_DOWN" if side=="TOP" else "MOVED_UP"
            if states[fid]==opposite:
                trace.append({"side":side,"face_id":fid,"status":"STOPPED_FACE_MOVED_OPPOSITE"});active[side]=False;continue
            fs=faces_now(vertices,edges,co);nxt,candidates=next_face(fid,side,fs,co,yc);adj=[x["face_id"] for x in fs[fid]["neighboring_faces"]];considered.append({"side":side,"face_id":fid,"next_face_id":nxt,"all_adjacent_faces":adj,"candidates":candidates})
            if nxt is None:
                trace.append({"side":side,"face_id":fid,"all_adjacent_faces":adj,"status":"STOPPED_NO_NEXT_BOUNDED_FACE"});active[side]=False;continue
            gap=face_gap(fs[fid],fs[nxt],co)
            if gap is None:
                trace.append({"side":side,"face_id":fid,"next_face_id":nxt,"status":"STOPPED_NO_VERTICAL_GAP"});active[side]=False;continue
            free=sorted(fs[fid]["free_vertices"]);movers=[v for v in free if vstate[v]=="UNMOVED"];sign=-1 if side=="TOP" else 1
            half_candidate=dict(co)
            for v in movers:half_candidate[v]=(co[v][0],co[v][1]+sign*gap["half_gap"])
            half_ok,half_diag,_,_,_,_,half_reconnect=wave.valid(vertices,edges,half_candidate)
            maximum,candidate,search=wave.maximum_valid_move(vertices,edges,co,movers,sign,gap["half_gap"]);actual=min(maximum,gap["half_gap"]);co=candidate;states[fid]="MOVED_UP" if side=="TOP" else "MOVED_DOWN"
            for v in movers:vstate[v]=states[fid]
            ok2,d2,a2,e2,c2,m2,re2=wave.valid(vertices,edges,co)
            if not ok2:raise AssertionError("invalid accepted state")
            failures=[]
            if not half_ok:
                if half_diag["proper_unrelated_edge_crossing_count"]:failures.append({"type":"PROPER_EDGE_CROSSING","details":half_diag["proper_unrelated_edge_crossings"]})
                if half_diag["coincident_vertex_pair_count"]:failures.append({"type":"COINCIDENT_VERTICES","details":half_diag["coincident_vertex_pairs"]})
                if half_diag["vertex_on_unrelated_edge_interior_count"]:failures.append({"type":"VERTEX_ON_UNRELATED_EDGE","details":half_diag["vertices_on_unrelated_edge_interiors"]})
                if not half_reconnect["valid"]:failures.append({"type":"CONNECTIVITY_OR_STALE_EDGE","details":half_reconnect["failures"]})
            if half_ok and abs(actual-gap["half_gap"])>geometry.GEOMETRY_TOLERANCE:raise AssertionError("hard-valid exact half move was reduced")
            if actual+geometry.GEOMETRY_TOLERANCE<gap["half_gap"] and not failures:raise AssertionError("reduction lacks hard-invalid target evidence")
            rec={"side":side,"face_id":fid,"next_face_id":nxt,"all_adjacent_faces":adj,"free_for_face_vertices":free,"moved_vertices":movers,**gap,"target_formula":"TARGET = full_gap / 2","target_uses_readability_clearance":False,"exact_half_hard_valid":half_ok,"exact_half_hard_failures":failures,"max_valid_move":maximum,"actual_move":actual,"reduction_reason":None if half_ok else failures,"search":search,"status":states[fid],"metrics_after":m2,"edge_reconnection":re2};trace.append(rec)
            stages.append({"stage":f"AFTER {side} {fid}","side":side,"face_id":fid,"moved_vertices":movers,"movement":actual,"coordinates":{v:list(p) for v,p in sorted(co.items())},"metrics":m2,"edge_reconnection":re2});progressed=True
            if side=="TOP":top=nxt
            else:bottom=nxt
        if not active["TOP"] and not active["BOTTOM"]:break
        if not progressed and (active["TOP"] or active["BOTTOM"]):raise AssertionError("active wave made no progress")
    else:raise AssertionError("did not terminate")
    top_moved=[x for x in trace if x.get("status")=="MOVED_UP"];bottom_moved=[x for x in trace if x.get("status")=="MOVED_DOWN"]
    c1=len(top_moved)>1 and "component:C1" in top_moved[1]["free_for_face_vertices"];c2=len(bottom_moved)>1 and "component:C2" in bottom_moved[1]["free_for_face_vertices"]
    report={"schema":"graph-relax.outward-face-wave.v1","source_coordinates":base["source_coordinates"],"Y_CENTER":yc,"vertex_count":34,"electrical_incidence_count":44,"initial_top_face":TOP_START,"initial_bottom_face":BOTTOM_START,"trace":trace,"faces_considered":considered,"face_states":states,"vertex_movement_states":vstate,"second_top_reaches_c1_face":c1,"second_bottom_reaches_c2_face":c2,"stages":stages,"final_coordinates":stages[-1]["coordinates"]}
    outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,stages,outdir/"iamp_outward_face_wave.svg");(outdir/"iamp_outward_face_wave_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");return report

def main():
    p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--base",type=Path,default=b/"iamp_corridor_separation_report.json");p.add_argument("--output",type=Path,default=b);q=p.parse_args();r=run(q.graph,q.base,q.output)
    print("initial",r["initial_top_face"],r["initial_bottom_face"])
    for x in r["trace"]:print(x["side"],x["face_id"],"->",x.get("next_face_id"),x["status"],"moved",x.get("moved_vertices"),"half",x.get("half_gap"),"max",x.get("max_valid_move"),"move",x.get("actual_move"),"crossings",x.get("metrics_after",{}).get("proper_crossings"),"violations",x.get("metrics_after",{}).get("clearance_violations"))
if __name__=="__main__":main()
