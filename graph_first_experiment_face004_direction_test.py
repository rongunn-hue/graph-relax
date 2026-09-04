#!/usr/bin/env python3
"""Controlled inward-versus-exact-outward face:004 direction test."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_sequential_local_face_pulls as sequential
import placement_geometry as geometry
import planar_pull_ownership as ownership

FACE_ID="face:004"
S=.90


def relation_key(item):return item["geometry_type"],item["item_a"],item["item_b"]


def evaluate(vertices,edges,coordinates):
    diagnostics,analysis,events,metrics=sequential.geometry_metrics(vertices,edges,coordinates)
    segments=ownership.rebuild_edge_segments(coordinates,edges);ownership.assert_reconnected(segments,coordinates,edges)
    return analysis,events,{**metrics,"vertex_count":len(vertices),"electrical_incidence_count":len(edges),
        "stale_edge_endpoint_count":0,"usable_zero_crossing":metrics["proper_crossings"]==0 and metrics["coincident_vertices"]==0 and metrics["vertex_on_unrelated_edge_events"]==0}


def event_comparison(original_event_doc,pre_analysis,a_analysis,b_analysis):
    original_ids={(e["violation_type"],e["item_a"],e["item_b"]):e["violation_id"] for e in original_event_doc["events"]}
    pre={relation_key(x):x for x in pre_analysis["relations"]};am={relation_key(x):x for x in a_analysis["relations"]};bm={relation_key(x):x for x in b_analysis["relations"]}
    baseline_keys=sorted(k for k,x in pre.items() if x["distance"]<geometry.READABILITY_CLEARANCE)
    def classification(before,after):
        if after>=geometry.READABILITY_CLEARANCE:return "DISAPPEARED"
        if after>before+geometry.GEOMETRY_TOLERANCE:return "IMPROVED"
        if after<before-geometry.GEOMETRY_TOLERANCE:return "WORSENED"
        return "UNCHANGED"
    records=[]
    for index,k in enumerate(baseline_keys,1):
        records.append({"event_id":original_ids.get(k,f"PRE004:{index:02d}"),"item_a":pre[k]["item_a"],"item_b":pre[k]["item_b"],
            "clearance_before":pre[k]["distance"],"clearance_after_case_a":am[k]["distance"],"case_a_classification":classification(pre[k]["distance"],am[k]["distance"]),
            "clearance_after_case_b":bm[k]["distance"],"case_b_classification":classification(pre[k]["distance"],bm[k]["distance"])})
    base=set(baseline_keys)
    def new_events(mapping):
        return [{"item_a":mapping[k]["item_a"],"item_b":mapping[k]["item_b"],"geometry_type":mapping[k]["geometry_type"],"clearance":mapping[k]["distance"],"deficit":mapping[k]["deficit"]}
                for k in sorted(k for k,x in mapping.items() if x["distance"]<geometry.READABILITY_CLEARANCE and k not in base)]
    return records,new_events(am),new_events(bm)


def render(vertices,edges,pre,cases,pole,movable,output):
    xs=[p[0] for p in pre.values()];ys=[p[1] for p in pre.values()];margin=100;minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>']
    stages=[("PRE-face:004",pre,None)]+[(case["label"],{v:tuple(p) for v,p in case["coordinates"].items()},case) for case in cases]
    for index,(label,coordinates,case) in enumerate(stages):
        begin=index*3;timing=f'begin="{begin}s" dur="3s"' if index<len(stages)-1 else f'begin="{begin}s" fill="freeze"';out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" {timing}/>')
        for a,b in edges:out.append(f'<line data-edge="{a}|{b}" x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        if case:
            for v in movable:
                x0,y0=pre[v];x1,y1=coordinates[v];out.append(f'<line x1="{x0:.17g}" y1="{y0:.17g}" x2="{x1:.17g}" y2="{y1:.17g}" stroke="#2e7d32" stroke-width="3" marker-end="url(#arrow)"/>')
        out.append(f'<circle cx="{pole[0]:.17g}" cy="{pole[1]:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/>')
        for v in sorted(vertices):
            x,y=coordinates[v];fill="#2e7d32" if v in movable else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            out.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in movable else 6}" fill="{fill}" stroke="white"/><text x="{x+9:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]}</text>')
        events=[] if case is None else case["crowding_events"]
        for event in events:
            x,y=event["crowding_point"];out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000"/>')
        metrics=cases[0]["pre_metrics"] if case is None else case["metrics"]
        out.append(f'<text x="{minx+16:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="19" font-weight="bold">{label}</text><text x="{minx+16:.17g}" y="{miny+52:.17g}" font-family="sans-serif" font-size="14">crossings={metrics["proper_crossings"]}; violations={metrics["clearance_violations"]}</text></g>')
    out.insert(2,'<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#2e7d32"/></marker></defs>');out.append('</svg>');output.write_text('\n'.join(out)+'\n')


def run(graph_path,direct_path,sequential_path,crowding_path,output_dir):
    vertices,edges,_,_,_=crowding.load_inputs(graph_path,direct_path);seq=json.loads(sequential_path.read_text());pre_stage=next(x for x in seq["stages"] if x["stage"]=="AFTER face:003");a_expected=next(x for x in seq["stages"] if x["stage"]=="AFTER face:004")
    pre={v:tuple(p) for v,p in pre_stage["coordinates"].items()};pole=tuple(a_expected["pole_coordinate"]);movable=a_expected["new_vertices_moved"]
    inward={v:(pole[0]+S*(pre[v][0]-pole[0]),pole[1]+S*(pre[v][1]-pole[1])) for v in movable}
    outward={v:(pre[v][0]+(1-S)*(pre[v][0]-pole[0]),pre[v][1]+(1-S)*(pre[v][1]-pole[1])) for v in movable}
    a=ownership.plant_destinations(pre,inward);b=ownership.plant_destinations(pre,outward)
    pre_analysis,pre_events,pre_metrics=evaluate(vertices,edges,pre);a_analysis,a_events,a_metrics=evaluate(vertices,edges,a);b_analysis,b_events,b_metrics=evaluate(vertices,edges,b)
    max_a_delta=max(math.dist(a[v],tuple(a_expected["coordinates"][v])) for v in a)
    if max_a_delta>geometry.GEOMETRY_TOLERANCE:raise AssertionError("Case A did not reproduce sequential face:004")
    displacement_a=sum(math.dist(pre[v],a[v]) for v in movable);displacement_b=sum(math.dist(pre[v],b[v]) for v in movable)
    if abs(displacement_a-displacement_b)>geometry.GEOMETRY_TOLERANCE:raise AssertionError("direction reversal changed magnitude")
    event_records,new_a,new_b=event_comparison(json.loads(crowding_path.read_text()),pre_analysis,a_analysis,b_analysis)
    cases=[]
    for label,coordinates,events,metrics,new in (("CASE A — INWARD",a,a_events,a_metrics,new_a),("CASE B — OUTWARD",b,b_events,b_metrics,new_b)):
        changes=sequential.relation_changes(pre_analysis,a_analysis if coordinates is a else b_analysis)
        cases.append({"label":label,"coordinates":{v:list(coordinates[v]) for v in sorted(coordinates)},"metrics":metrics,"pre_metrics":pre_metrics,"crowding_events":events,
            "removed_violations":changes["removed"],"improved_violations":changes["improved"],"worsened_violations":changes["worsened"],"new_violations":new})
    report={"schema":"graph-relax.face004-direction-test.v1","source_geometry":f"{sequential_path}#AFTER face:003","face_id":FACE_ID,"pull_factor":S,"pole_coordinate":list(pole),"movable_vertices":movable,
        "ownership_before":pre_stage["ownership_table"],"pre_coordinates":{v:list(pre[v]) for v in sorted(pre)},"pre_metrics":pre_metrics,"case_a":cases[0],"case_b":cases[1],
        "event_comparison":event_records,"total_displacement_case_a":displacement_a,"total_displacement_case_b":displacement_b,"displacement_magnitudes_equal":abs(displacement_a-displacement_b)<=geometry.GEOMETRY_TOLERANCE,
        "case_a_reproduces_sequential_result":max_a_delta<=geometry.GEOMETRY_TOLERANCE,"maximum_case_a_coordinate_difference":max_a_delta}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,pre,cases,pole,movable,output_dir/'iamp_face004_direction_test.svg');(output_dir/'iamp_face004_direction_test_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    p=argparse.ArgumentParser();base=Path('output/graph_first');p.add_argument('--graph',type=Path,default=base/'iamp_graph.json');p.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json');p.add_argument('--sequential-report',type=Path,default=base/'iamp_sequential_local_face_pulls_report.json');p.add_argument('--crowding-report',type=Path,default=base/'iamp_crowding_points_report.json');p.add_argument('--output',type=Path,default=base);a=p.parse_args();r=run(a.graph,a.direct_report,a.sequential_report,a.crowding_report,a.output)
    for name in ('pre_metrics','case_a','case_b'):
        value=r[name];metrics=value if name=='pre_metrics' else value['metrics'];print(name,metrics)
    print('displacements',r['total_displacement_case_a'],r['total_displacement_case_b'])

if __name__=='__main__':main()
