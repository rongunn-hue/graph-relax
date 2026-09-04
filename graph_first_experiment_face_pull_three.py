#!/usr/bin/env python3
"""Independent and cumulative fixed 10% topology-free face-pole pulls."""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_face_pull_rmin1_free as free_pull
import placement_geometry as geometry

CONTRACTION_SCALE=.90
FACE_ORDER=("RMIN1","RA_SW","RMIN4")


def classify_face(walk,edges):
    boundary_edges={tuple(sorted((walk[i],walk[(i+1)%len(walk)]))) for i in range(len(walk))}
    incident=defaultdict(set)
    for edge in edges:
        edge=tuple(sorted(edge))
        for vertex in edge: incident[vertex].add(edge)
    records=[];free=[];anchors=[]
    for vertex in dict.fromkeys(walk):
        face_edges={edge for edge in boundary_edges if vertex in edge}
        is_free=len(incident[vertex])==2 and len(face_edges)==2 and incident[vertex]==face_edges
        classification="FREE" if is_free else "ANCHOR"
        (free if is_free else anchors).append(vertex)
        records.append({"vertex":vertex,"classification":classification,"degree":len(incident[vertex]),
                        "incident_edges":[list(edge) for edge in sorted(incident[vertex])],
                        "face_boundary_edges":[list(edge) for edge in sorted(face_edges)],
                        "additional_incident_edges":[list(edge) for edge in sorted(incident[vertex]-face_edges)]})
    return records,free,anchors


def apply_pull(coordinates,free,pole):
    result=dict(coordinates)
    for vertex in free:
        x,y=coordinates[vertex]
        result[vertex]=(pole[0]+CONTRACTION_SCALE*(x-pole[0]),pole[1]+CONTRACTION_SCALE*(y-pole[1]))
    return result


def unique_locations(events):
    points=[]
    for event in events:
        point=tuple(event["crowding_point"])
        if not any(math.dist(point,prior)<=geometry.GEOMETRY_TOLERANCE for prior in points):points.append(point)
    return len(points)


def validate_reconnection(coordinates,edges):
    rendered=[(edge,coordinates[edge[0]],coordinates[edge[1]]) for edge in edges]
    stale=sum(start!=coordinates[edge[0]] or end!=coordinates[edge[1]] for edge,start,end in rendered)
    if stale:raise AssertionError("stale edge endpoint")
    return stale


def measure(vertices,edges,coordinates,baseline):
    stale=validate_reconnection(coordinates,edges)
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    analysis=crowding.analyze(vertices,edges,coordinates)
    events=crowding_points.localize_events(analysis);clusters,_,_=crowding_points.cluster_events(events)
    changes,new_items=free_pull.original_violation_changes(baseline,analysis)
    return {"proper_crossings":diagnostics["proper_unrelated_edge_crossing_count"],"proper_crossing_list":diagnostics["proper_unrelated_edge_crossings"],
            "coincident_vertices":diagnostics["coincident_vertex_pair_count"],"coincident_vertex_pairs":diagnostics["coincident_vertex_pairs"],
            "vertex_on_unrelated_edge_events":diagnostics["vertex_on_unrelated_edge_interior_count"],"vertex_on_unrelated_edge_event_list":diagnostics["vertices_on_unrelated_edge_interiors"],
            "minimum_node_distance":diagnostics["minimum_node_distance"],"stale_edge_endpoint_count":stale,
            "clearance_violations":len(events),"unique_crowding_event_locations":unique_locations(events),
            "crowding_clusters":len(clusters),"worst_clearance":min((e["clearance"] for e in events),default=None),
            "largest_clearance_deficit":max((e["deficit"] for e in events),default=0.0),
            "original_violation_changes":changes,"violations_removed":sum(x["removed_as_violation"] for x in changes),
            "violations_improved":sum(x["classification"]=="IMPROVED" for x in changes),
            "violations_worsened":sum(x["classification"]=="DEGRADED" for x in changes),
            "new_violations":new_items,
            "abstract_graph_valid":diagnostics["proper_unrelated_edge_crossing_count"]==0 and diagnostics["coincident_vertex_pair_count"]==0 and diagnostics["vertex_on_unrelated_edge_interior_count"]==0}


def stage_record(label,face,coordinates,free,anchors,metrics,before=None,accepted=True):
    displacement=sum(math.dist(before[v],coordinates[v]) for v in free) if before else 0.0
    return {"stage":label,"face":face["visual_location"] if face else None,"face_id":face["face_id"] if face else None,
            "pole_coordinate":face["pole_coordinate"] if face else None,"free_vertices":free,"anchor_vertices":anchors,
            "accepted":accepted,"total_free_vertex_displacement":displacement,
            "coordinates":{v:list(coordinates[v]) for v in sorted(coordinates)},**metrics}


def render(vertices,edges,original,faces,stages,title):
    xs=[p[0] for p in original.values()];ys=[p[1] for p in original.values()];margin=100
    minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>']
    by_name={face["visual_location"]:face for face in faces}
    for index,stage in enumerate(stages):
        begin=index*2.5;timing=f'begin="{begin}s" dur="2.5s"' if index<len(stages)-1 else f'begin="{begin}s" fill="freeze"'
        out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" {timing}/>')
        co={v:tuple(p) for v,p in stage["coordinates"].items()};face=by_name.get(stage["face"]);boundary=list(dict.fromkeys(face["ordered_boundary_cycle"])) if face else []
        boundary_edges={frozenset((face["ordered_boundary_cycle"][i],face["ordered_boundary_cycle"][(i+1)%len(face["ordered_boundary_cycle"])])) for i in range(len(face["ordered_boundary_cycle"]))} if face else set()
        if face:
            out.append('<polygon points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in face["ordered_boundary_cycle"])+ '" fill="#90caf9" fill-opacity=".18"/>')
        for a,b in edges:
            color="#1565c0" if frozenset((a,b)) in boundary_edges else "#777";sw=4 if color=="#1565c0" else 1.5
            out.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="{color}" stroke-width="{sw}"/>')
        if face:
            px,py=face["pole_coordinate"]
            for v in stage["free_vertices"]:out.append(f'<line x1="{px:.17g}" y1="{py:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#2e7d32" stroke-dasharray="6 5"/>')
            out.append(f'<circle cx="{px:.17g}" cy="{py:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/>')
        for vertex in sorted(vertices):
            x,y=co[vertex];fill="#2e7d32" if vertex in stage["free_vertices"] else "#ef6c00" if vertex in stage["anchor_vertices"] else "#1565c0" if vertices[vertex]["type"]=="COMPONENT" else "#d84315";radius=9 if vertex in stage["free_vertices"] or vertex in stage["anchor_vertices"] else 6
            out.append(f'<circle data-vertex="{vertex}" cx="{x:.17g}" cy="{y:.17g}" r="{radius}" fill="{fill}" stroke="white"/><text x="{x+9:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[vertex]["name"]}</text>')
        status="VALID" if stage["abstract_graph_valid"] else "INVALID"
        out.append(f'<text x="{minx+16:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="19" font-weight="bold">{title}: {stage["stage"]} — {status}</text><text x="{minx+16:.17g}" y="{miny+52:.17g}" font-family="sans-serif" font-size="14">crossings={stage["proper_crossings"]}; violations={stage["clearance_violations"]}; green=FREE; orange=ANCHOR</text></g>')
    out.append('</svg>');return '\n'.join(out)+'\n'


def run(graph_path,direct_path,face_path,output_dir):
    vertices,edges,original,_,_=crowding.load_inputs(graph_path,direct_path);face_doc=json.loads(face_path.read_text());faces=face_doc["selected_faces"]
    by_name={face["visual_location"]:face for face in faces};baseline=crowding.analyze(vertices,edges,original);base_metrics=measure(vertices,edges,original,baseline)
    topology={};independent=[]
    independent_stages=[stage_record("ORIGINAL",None,original,[],[],base_metrics)]
    for name in FACE_ORDER:
        face=by_name[name];classification,free,anchors=classify_face(face["ordered_boundary_cycle"],edges);topology[name]={"face_id":face["face_id"],"ordered_boundary_walk":face["ordered_boundary_cycle"],"distinct_boundary_vertex_count":len(set(face["ordered_boundary_cycle"])),"area":face["area"],"classification":classification,"free_vertices":free,"anchor_vertices":anchors}
        candidate=apply_pull(original,free,tuple(face["pole_coordinate"]));metrics=measure(vertices,edges,candidate,baseline)
        record=stage_record(f"{name}-ONLY s=0.90",face,candidate,free,anchors,metrics,original);independent.append(record);independent_stages.append(record)
    current=dict(original);cumulative=[stage_record("ORIGINAL",None,current,[],[],base_metrics)];operations=[]
    for name in FACE_ORDER:
        face=by_name[name];classification,free,anchors=classify_face(face["ordered_boundary_cycle"],edges)
        if free!=topology[name]["free_vertices"] or anchors!=topology[name]["anchor_vertices"]:raise AssertionError("topological classification changed")
        before=dict(current);candidate=apply_pull(current,free,tuple(face["pole_coordinate"]));candidate_metrics=measure(vertices,edges,candidate,baseline)
        accepted=candidate_metrics["abstract_graph_valid"]
        if accepted:current=candidate
        accepted_metrics=measure(vertices,edges,current,baseline)
        record=stage_record(f"AFTER {name}",face,current,free,anchors,accepted_metrics,before,accepted);record["candidate_diagnostics"]=candidate_metrics;record["rejected_and_restored"]=not accepted
        operations.append(record);cumulative.append(record)
    final=cumulative[-1]
    report={"schema":"graph-relax.three-face-free-pulls.v1","source_coordinates":f"{direct_path}#final_coordinates","contraction_scale":CONTRACTION_SCALE,"readability_clearance":geometry.READABILITY_CLEARANCE,"vertex_count":len(vertices),"electrical_incidence_count":len(edges),"face_order":list(FACE_ORDER),"topological_classification":topology,"original":independent_stages[0],"independent_tests":independent,"independent_comparison":[{"face":r["face"],"boundary_vertex_count":topology[r["face"]]["distinct_boundary_vertex_count"],"free_vertex_count":len(r["free_vertices"]),"anchor_count":len(r["anchor_vertices"]),"face_area":topology[r["face"]]["area"],"clearance_violations_before":base_metrics["clearance_violations"],"clearance_violations_after":r["clearance_violations"],"violations_removed":r["violations_removed"],"new_violations":len(r["new_violations"]),"crossings_after":r["proper_crossings"],"total_free_vertex_displacement":r["total_free_vertex_displacement"]} for r in independent],"cumulative_stages":cumulative,"cumulative_operations":operations,"final":final}
    output_dir.mkdir(parents=True,exist_ok=True);(output_dir/'iamp_face_pull_three_independent.svg').write_text(render(vertices,edges,original,faces,independent_stages,'Independent face pulls'))
    (output_dir/'iamp_face_pull_three_cumulative.svg').write_text(render(vertices,edges,original,faces,cumulative,'Cumulative face pulls'))
    (output_dir/'iamp_face_pull_three_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    parser=argparse.ArgumentParser();base=Path('output/graph_first');parser.add_argument('--graph',type=Path,default=base/'iamp_graph.json');parser.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json');parser.add_argument('--face-report',type=Path,default=base/'iamp_face_poles_report.json');parser.add_argument('--output',type=Path,default=base);a=parser.parse_args();r=run(a.graph,a.direct_report,a.face_report,a.output)
    for name in FACE_ORDER:print(name,'FREE',r['topological_classification'][name]['free_vertices'],'ANCHOR',r['topological_classification'][name]['anchor_vertices'])
    for x in r['independent_comparison']:print(x)
    for x in r['cumulative_stages']:print(x['stage'],x['proper_crossings'],x['clearance_violations'],x.get('accepted'))

if __name__=='__main__':main()
