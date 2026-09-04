#!/usr/bin/env python3
"""Face:004 crowding-side to unused-interior directed-opening diagnostic."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_sequential_local_face_pulls as sequential
import placement_geometry as geometry
import planar_pull_ownership as ownership

FACE_ID="face:004"
REFERENCE_FRACTION=.10


def relation_key(item):return item["geometry_type"],item["item_a"],item["item_b"]


def item_touches_face(item,boundary_vertices,boundary_edges):
    if item.startswith("vertex:"):return item[len("vertex:"):] in boundary_vertices
    edge=tuple(sorted(item[len("edge:"):].split("|")))
    return edge in boundary_edges or bool(set(edge)&boundary_vertices)


def evaluate(vertices,edges,coordinates):
    diagnostics,analysis,events,metrics=sequential.geometry_metrics(vertices,edges,coordinates)
    segments=ownership.rebuild_edge_segments(coordinates,edges);ownership.assert_reconnected(segments,coordinates,edges)
    return analysis,events,{**metrics,"vertex_count":len(vertices),"electrical_incidence_count":len(edges),"stale_edge_endpoint_count":0,
        "usable_zero_crossing":metrics["proper_crossings"]==0 and metrics["coincident_vertices"]==0 and metrics["vertex_on_unrelated_edge_events"]==0}


def classify(before,after):
    if after>=geometry.READABILITY_CLEARANCE:return "DISAPPEARED"
    if after>before+geometry.GEOMETRY_TOLERANCE:return "IMPROVED"
    if after<before-geometry.GEOMETRY_TOLERANCE:return "WORSENED"
    return "UNCHANGED"


def event_report(original_events,baseline,final):
    ids={(e["violation_type"],e["item_a"],e["item_b"]):e["violation_id"] for e in original_events["events"]}
    before={relation_key(x):x for x in baseline["relations"]};after={relation_key(x):x for x in final["relations"]};baseline_bad=sorted(k for k,x in before.items() if x["distance"]<geometry.READABILITY_CLEARANCE)
    records=[]
    for index,k in enumerate(baseline_bad,1):
        records.append({"event_id":ids.get(k,f"PRE004:{index:02d}"),"item_a":before[k]["item_a"],"item_b":before[k]["item_b"],"baseline_clearance":before[k]["distance"],"resulting_clearance":after[k]["distance"],"classification":classify(before[k]["distance"],after[k]["distance"]),
            "disappeared":after[k]["distance"]>=geometry.READABILITY_CLEARANCE,"improved":after[k]["distance"]>before[k]["distance"]+geometry.GEOMETRY_TOLERANCE,"unchanged":abs(after[k]["distance"]-before[k]["distance"])<=geometry.GEOMETRY_TOLERANCE,"worsened":after[k]["distance"]<before[k]["distance"]-geometry.GEOMETRY_TOLERANCE})
    base=set(baseline_bad);new=[{"item_a":after[k]["item_a"],"item_b":after[k]["item_b"],"geometry_type":after[k]["geometry_type"],"clearance":after[k]["distance"],"deficit":after[k]["deficit"]} for k in sorted(k for k,x in after.items() if x["distance"]<geometry.READABILITY_CLEARANCE and k not in base)]
    return records,new


def panel(out,index,label,vertices,edges,coordinates,face_walk,pole,crowding_center,movable,movements,events,minx,miny,width,height,metrics,show_diagnostics):
    offset=index*width;out.append(f'<g transform="translate({offset-minx:.17g},{-miny:.17g})">')
    points=' '.join(f'{coordinates[v][0]:.17g},{coordinates[v][1]:.17g}' for v in face_walk);out.append(f'<polygon points="{points}" fill="#90caf9" fill-opacity=".22" stroke="#1565c0" stroke-width="3"/>')
    for a,b in edges:out.append(f'<line data-edge="{a}|{b}" x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
    if show_diagnostics:
        out.append(f'<line x1="{crowding_center[0]:.17g}" y1="{crowding_center[1]:.17g}" x2="{pole[0]:.17g}" y2="{pole[1]:.17g}" stroke="#6a1b9a" stroke-width="4" marker-end="url(#arrow)"/>')
        out.append(f'<circle cx="{crowding_center[0]:.17g}" cy="{crowding_center[1]:.17g}" r="9" fill="#ffca28" stroke="#000"/><text x="{crowding_center[0]+11:.17g}" y="{crowding_center[1]-10:.17g}" font-family="sans-serif" font-size="12">C004</text>')
        for movement in movements:
            v=movement["vertex"];old=movement["original_coordinate"];new=movement["new_coordinate"]
            out.append(f'<line x1="{old[0]:.17g}" y1="{old[1]:.17g}" x2="{new[0]:.17g}" y2="{new[1]:.17g}" stroke="#2e7d32" stroke-width="3" marker-end="url(#arrow-green)"/>')
    out.append(f'<circle cx="{pole[0]:.17g}" cy="{pole[1]:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/><text x="{pole[0]+11:.17g}" y="{pole[1]-10:.17g}" font-family="sans-serif" font-size="12">P004</text>')
    for v in sorted(vertices):
        x,y=coordinates[v];fill="#2e7d32" if v in movable else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315";out.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in movable else 6}" fill="{fill}" stroke="white"/><text x="{x+9:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]}</text>')
    for event in events:
        x,y=event["crowding_point"];out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000"/>')
    out.append(f'<text x="{minx+15:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="17" font-weight="bold">{label}</text><text x="{minx+15:.17g}" y="{miny+50:.17g}" font-family="sans-serif" font-size="13">crossings={metrics["proper_crossings"]}; violations={metrics["clearance_violations"]}</text></g>')


def render(vertices,edges,pre,inward,opened,walk,pole,crowding_center,movable,movements,pre_events,inward_events,opened_events,metrics,output):
    all_points=list(pre.values())+list(inward.values())+list(opened.values());xs=[p[0] for p in all_points];ys=[p[1] for p in all_points];margin=90;minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {3*width:.17g} {height:.17g}" width="2100" height="800">','<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#6a1b9a"/></marker><marker id="arrow-green" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#2e7d32"/></marker></defs>',f'<rect width="{3*width:.17g}" height="{height:.17g}" fill="white"/>']
    panel(out,0,"PRE-face:004",vertices,edges,pre,walk,pole,crowding_center,movable,[],pre_events,minx,miny,width,height,metrics["pre"],False)
    panel(out,1,"PREVIOUS RADIAL INWARD",vertices,edges,inward,walk,pole,crowding_center,movable,[],inward_events,minx,miny,width,height,metrics["inward"],False)
    panel(out,2,"CROWDING-DIRECTED OPENING",vertices,edges,opened,walk,pole,crowding_center,movable,movements,opened_events,minx,miny,width,height,metrics["opened"],True)
    out.append('</svg>');output.write_text('\n'.join(out)+'\n')


def run(graph_path,direct_path,sequential_path,inventory_path,prior_direction_path,original_event_path,output_dir):
    vertices,edges,_,_,_=crowding.load_inputs(graph_path,direct_path);seq=json.loads(sequential_path.read_text());pre_stage=next(x for x in seq["stages"] if x["stage"]=="AFTER face:003");pre={v:tuple(p) for v,p in pre_stage["coordinates"].items()}
    prior=json.loads(prior_direction_path.read_text());inward={v:tuple(p) for v,p in prior["case_a"]["coordinates"].items()};movable=prior["movable_vertices"]
    face=next(x for x in json.loads(inventory_path.read_text())["faces"] if x["face_id"]==FACE_ID);walk=face["ordered_boundary_walk"];boundary_vertices=set(face["distinct_boundary_vertices"]);boundary_edges={tuple(edge) for edge in face["distinct_boundary_edges"]}
    pre_analysis,pre_events,pre_metrics=evaluate(vertices,edges,pre)
    if len(vertices)!=34 or len(edges)!=44 or pre_metrics["proper_crossings"]!=0 or pre_metrics["clearance_violations"]!=6:raise AssertionError("pre-face:004 baseline mismatch")
    qualifying=[event for event in pre_events if item_touches_face(event["item_a"],boundary_vertices,boundary_edges) or item_touches_face(event["item_b"],boundary_vertices,boundary_edges)]
    if not qualifying:raise RuntimeError("face:004 has no measurable crowded-side events")
    center=(sum(e["crowding_point"][0] for e in qualifying)/len(qualifying),sum(e["crowding_point"][1] for e in qualifying)/len(qualifying))
    points=[pre[v] for v in walk];pole,pole_clearance,_=face_poles.safe_pole(points);raw=(pole[0]-center[0],pole[1]-center[1]);length=math.hypot(*raw);direction=(raw[0]/length,raw[1]/length)
    destinations={};movements=[]
    for v in movable:
        old=pre[v];radial=(REFERENCE_FRACTION*(pole[0]-old[0]),REFERENCE_FRACTION*(pole[1]-old[1]));magnitude=math.hypot(*radial);directed=(magnitude*direction[0],magnitude*direction[1]);new=(old[0]+directed[0],old[1]+directed[1]);destinations[v]=new
        if abs(math.hypot(*directed)-magnitude)>geometry.GEOMETRY_TOLERANCE:raise AssertionError("displacement magnitude changed")
        movements.append({"vertex":v,"original_coordinate":list(old),"old_10_percent_radial_displacement_vector":list(radial),"old_radial_displacement_magnitude":magnitude,"new_directed_opening_displacement_vector":list(directed),"new_directed_displacement_magnitude":math.hypot(*directed),"new_coordinate":list(new)})
    opened=ownership.plant_destinations(pre,destinations);segments=ownership.rebuild_edge_segments(opened,edges);ownership.assert_reconnected(segments,opened,edges)
    inward_analysis,inward_events,inward_metrics=evaluate(vertices,edges,inward);opened_analysis,opened_events,opened_metrics=evaluate(vertices,edges,opened)
    event_records,new_events=event_report(json.loads(original_event_path.read_text()),pre_analysis,opened_analysis)
    comparison={"PRE_004":pre_metrics,"RADIAL_INWARD":prior["case_a"]["metrics"],"RADIAL_OUTWARD":prior["case_b"]["metrics"],"DIRECTED_OPENING":opened_metrics}
    report={"schema":"graph-relax.face004-directed-opening.v1","source_geometry":f"{sequential_path}#AFTER face:003","face_id":FACE_ID,"face_boundary_walk":walk,"movable_vertices":movable,"qualifying_crowding_events":qualifying,"qualifying_event_count":len(qualifying),
        "crowding_side_representative_C004":list(center),"safe_interior_point_P004":list(pole),"P004_boundary_clearance":pole_clearance,"d_raw":list(raw),"d_raw_magnitude":length,"d_hat":list(direction),"movement_reference_fraction":REFERENCE_FRACTION,"movements":movements,
        "pre_metrics":pre_metrics,"directed_opening_metrics":opened_metrics,"baseline_event_results":event_records,"new_violations":new_events,"comparison":comparison,
        "final_coordinates":{v:list(opened[v]) for v in sorted(opened)},"complete_rebuilt_edge_segments":segments,"operation_valid":opened_metrics["usable_zero_crossing"]}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,pre,inward,opened,walk,pole,center,movable,movements,pre_events,inward_events,opened_events,{"pre":pre_metrics,"inward":inward_metrics,"opened":opened_metrics},output_dir/'iamp_face004_directed_opening.svg');(output_dir/'iamp_face004_directed_opening_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct-report',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--sequential-report',type=Path,default=b/'iamp_sequential_local_face_pulls_report.json');p.add_argument('--inventory',type=Path,default=b/'iamp_all_faces_report.json');p.add_argument('--prior-direction-report',type=Path,default=b/'iamp_face004_direction_test_report.json');p.add_argument('--original-events',type=Path,default=b/'iamp_crowding_points_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct_report,a.sequential_report,a.inventory,a.prior_direction_report,a.original_events,a.output);print(json.dumps({k:r[k] for k in ('crowding_side_representative_C004','safe_interior_point_P004','d_raw','d_hat','movements','comparison','baseline_event_results','new_violations','operation_valid')},indent=2))

if __name__=='__main__':main()
