#!/usr/bin/env python3
"""Audit and, only if unambiguous, test one two-sided crowded-seam opening."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_sequential_local_face_pulls as sequential
import placement_geometry as geometry


def point_in_bounded_face(point,face,coordinates):
    walk=face["ordered_boundary_walk"];points=[coordinates[v] for v in walk];inside=False
    if min(geometry.point_segment_distance(point,points[i],points[(i+1)%len(points)]) for i in range(len(points)))<=geometry.GEOMETRY_TOLERANCE:return False
    for i,first in enumerate(points):
        second=points[(i+1)%len(points)]
        if (first[1]>point[1]) != (second[1]>point[1]):
            x=first[0]+(point[1]-first[1])*(second[0]-first[0])/(second[1]-first[1])
            if x>point[0]:inside=not inside
    return inside


def incident_faces(item,faces):
    if item.startswith("vertex:"):
        vertex=item[len("vertex:"):];return [f for f in faces if vertex in f["distinct_boundary_vertices"]]
    edge=tuple(sorted(item[len("edge:"):].split("|")));return [f for f in faces if edge in [tuple(e) for e in f["distinct_boundary_edges"]]]


def local_gap_faces(item,event_midpoint,faces,coordinates):
    incident=incident_faces(item,faces)
    containing=[f for f in incident if f["is_bounded"] and point_in_bounded_face(event_midpoint,f,coordinates)]
    return incident,containing


def render(vertices,edges,coordinates,faces,events,audits,output):
    xs=[p[0] for p in coordinates.values()];ys=[p[1] for p in coordinates.values()];margin=100;minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',f'<text x="{minx+16:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="19" font-weight="bold">NO ELIGIBLE TWO-SIDED SEAM — GRAPH UNCHANGED</text>']
    colors=["#bbdefb","#c8e6c9","#ffe0b2","#e1bee7","#b2dfdb","#f8bbd0","#dcedc8","#fff9c4","#d1c4e9","#b3e5fc","#ffccbc"]
    for index,face in enumerate(f for f in faces if f["is_bounded"] and f["simple_boundary"]):
        points=' '.join(f'{coordinates[v][0]:.17g},{coordinates[v][1]:.17g}' for v in face["ordered_boundary_walk"]);out.append(f'<polygon points="{points}" fill="{colors[index%len(colors)]}" fill-opacity=".16"/>')
    for a,b in edges:out.append(f'<line data-edge="{a}|{b}" x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
    for face in faces:
        if not face["is_bounded"] or not face["simple_boundary"]:continue
        point,_,_=face_poles.safe_pole([coordinates[v] for v in face["ordered_boundary_walk"]]);out.append(f'<text x="{point[0]:.17g}" y="{point[1]:.17g}" font-family="sans-serif" font-size="11" fill="#555">{face["face_id"]}</text>')
    for vertex in sorted(vertices):
        x,y=coordinates[vertex];fill="#1565c0" if vertices[vertex]["type"]=="COMPONENT" else "#d84315";out.append(f'<circle data-vertex="{vertex}" cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{fill}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[vertex]["name"]}</text>')
    for event,audit in zip(events,audits):
        p1,p2=event["pA"],event["pB"];x,y=event["crowding_point"];out.append(f'<line x1="{p1[0]:.17g}" y1="{p1[1]:.17g}" x2="{p2[0]:.17g}" y2="{p2[1]:.17g}" stroke="#d50000" stroke-width="2"/><circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="#d50000" stroke="white"/><text x="{x+8:.17g}" y="{y+14:.17g}" font-family="sans-serif" font-size="11" font-weight="bold" fill="#d50000">{event["violation_id"]}: {audit["eligibility"]}</text>')
    out.append('</svg>');output.write_text('\n'.join(out)+'\n')


def run(graph_path,direct_path,sequential_path,output_dir):
    vertices,edges,_,_,_=crowding.load_inputs(graph_path,direct_path);seq=json.loads(sequential_path.read_text());pre_stage=next(s for s in seq["stages"] if s["stage"]=="AFTER face:003");coordinates={v:tuple(p) for v,p in pre_stage["coordinates"].items()};original=dict(coordinates)
    diagnostics,analysis,events,metrics=sequential.geometry_metrics(vertices,edges,coordinates)
    if len(vertices)!=34 or len(edges)!=44 or metrics["proper_crossings"]!=0 or metrics["clearance_violations"]!=6:raise AssertionError("common pre-face:004 baseline mismatch")
    faces,dart_face=inventory.enumerate_faces(vertices,edges,coordinates);inventory.add_classification_and_adjacency(faces,edges,dart_face)
    owner={row["vertex"]:row["owner_face"] for row in pre_stage["ownership_table"] if row["owner_face"]}
    audits=[];eligible=[]
    for event in sorted(events,key=lambda e:(e["clearance"],e["item_a"],e["item_b"])):
        midpoint=tuple(event["crowding_point"]);incident_a,local_a=local_gap_faces(event["item_a"],midpoint,faces,coordinates);incident_b,local_b=local_gap_faces(event["item_b"],midpoint,faces,coordinates)
        reasons=[]
        if len(local_a)!=1:reasons.append(f"item A local face is {'unassigned' if not local_a else 'ambiguous'}")
        if len(local_b)!=1:reasons.append(f"item B local face is {'unassigned' if not local_b else 'ambiguous'}")
        face_a=local_a[0] if len(local_a)==1 else None;face_b=local_b[0] if len(local_b)==1 else None
        if face_a and face_b and face_a["face_id"]==face_b["face_id"]:reasons.append("both geometries open into the same local face")
        candidates_a=sorted(v for v in (face_a["free_vertices"] if face_a else []) if v not in owner);candidates_b=sorted(v for v in (face_b["free_vertices"] if face_b else []) if v not in owner)
        conflicts=sorted(set(candidates_a)&set(candidates_b));resolved_a=list(candidates_a);resolved_b=list(candidates_b)
        for vertex in conflicts:
            # Both candidate sides are at the same diagnostic layer; stable
            # face ID is the permanent ownership tie-break.
            winner=min(face_a["face_id"],face_b["face_id"])
            if winner==face_a["face_id"]:resolved_b.remove(vertex)
            else:resolved_a.remove(vertex)
        if face_a and face_b and face_a["face_id"]!=face_b["face_id"]:
            if not resolved_a:reasons.append("Face A has no unowned exclusive FREE vertex")
            if not resolved_b:reasons.append("Face B has no unowned exclusive FREE vertex")
        status="ELIGIBLE" if not reasons else "INELIGIBLE"
        audit={"event_id":event["violation_id"],"event_type":event["violation_type"],"clearance":event["clearance"],"item_a":event["item_a"],"item_b":event["item_b"],"pA":event["pA"],"pB":event["pB"],"event_midpoint":event["crowding_point"],
            "item_a_incident_faces":[f["face_id"] for f in incident_a],"item_b_incident_faces":[f["face_id"] for f in incident_b],"item_a_local_gap_faces":[f["face_id"] for f in local_a],"item_b_local_gap_faces":[f["face_id"] for f in local_b],
            "face_a":face_a["face_id"] if face_a else None,"face_b":face_b["face_id"] if face_b else None,"face_a_unowned_free_candidates":candidates_a,"face_b_unowned_free_candidates":candidates_b,"ownership_conflicts":conflicts,"resolved_side_a_movables":resolved_a,"resolved_side_b_movables":resolved_b,"eligibility":status,"ineligibility_reasons":reasons}
        audits.append(audit)
        if status=="ELIGIBLE":eligible.append((event,audit))
    selected=eligible[0] if eligible else None
    if selected is not None:raise AssertionError("this data unexpectedly produced an eligible event; movement implementation intentionally absent from this stopped diagnostic")
    if coordinates!=original:raise AssertionError("eligibility audit moved graph geometry")
    report={"schema":"graph-relax.single-seam-two-sided-opening.v1","source_geometry":f"{sequential_path}#AFTER face:003","status":"STOPPED — NO CURRENT EVENT SATISFIES TWO-SIDED ELIGIBILITY","primitive_tested":False,"selected_event":None,
        "vertex_count":len(vertices),"electrical_incidence_count":len(edges),"proper_crossings":metrics["proper_crossings"],"clearance_violations":metrics["clearance_violations"],"coordinates_unchanged":True,"persistent_owner_face_by_vertex":owner,"event_eligibility_audit":audits,
        "conclusion":"No baseline violation has two unambiguous, distinct local gap faces that both retain an exclusive unowned FREE vertex after persistent ownership."}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,coordinates,faces,sorted(events,key=lambda e:(e["clearance"],e["item_a"],e["item_b"])),audits,output_dir/'iamp_single_seam_two_sided_opening.svg');(output_dir/'iamp_single_seam_two_sided_opening_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct-report',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--sequential-report',type=Path,default=b/'iamp_sequential_local_face_pulls_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct_report,a.sequential_report,a.output);print(r['status']);
    for event in r['event_eligibility_audit']:print(event['event_id'],event['clearance'],event['face_a'],event['face_b'],event['ineligibility_reasons'])

if __name__=='__main__':main()
