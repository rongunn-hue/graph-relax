#!/usr/bin/env python3
"""Diagnostic contraction of only topologically free RMIN1-face vertices."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_face_pull_rmin1 as full_pull
import placement_geometry as geometry

SCALES=(1.0,.9,.8,.7,.6,.5)


def classify_boundary(boundary,edges):
    incident={node:[] for node in boundary}
    for edge in edges:
        for node in edge:
            if node in incident: incident[node].append(edge)
    records=[]; free=[]; anchors=[]
    for index,node in enumerate(boundary):
        expected={frozenset((node,boundary[index-1])),frozenset((node,boundary[(index+1)%len(boundary)]))}
        actual={frozenset(edge) for edge in incident[node]}
        external=sorted([list(edge) for edge in incident[node] if frozenset(edge) not in expected])
        kind="FREE" if len(actual)==2 and actual==expected else "ANCHOR"
        (free if kind=="FREE" else anchors).append(node)
        records.append({"vertex":node,"classification":kind,"degree":len(incident[node]),
                        "incident_edges":[list(edge) for edge in sorted(incident[node])],
                        "face_boundary_edges":[[node,boundary[index-1]],[node,boundary[(index+1)%len(boundary)]]],
                        "additional_incident_edges":external})
    return records,free,anchors


def snapshot(original,free,pole,scale):
    result=dict(original)
    for node in free:
        x,y=original[node]; result[node]=(pole[0]+scale*(x-pole[0]),pole[1]+scale*(y-pole[1]))
    return result


def unique_event_locations(events):
    points=[]
    for event in events:
        point=tuple(event["crowding_point"])
        if not any(math.dist(point,existing)<=geometry.GEOMETRY_TOLERANCE for existing in points): points.append(point)
    return len(points)


def original_violation_changes(baseline,current):
    old={full_pull.relation_key(item):item for item in baseline["relations"]}
    new={full_pull.relation_key(item):item for item in current["relations"]}
    old_bad=sorted(key for key,item in old.items() if item["distance"]<geometry.READABILITY_CLEARANCE)
    changes=[]
    for key in old_bad:
        delta=new[key]["distance"]-old[key]["distance"]
        classification=("IMPROVED" if delta>geometry.GEOMETRY_TOLERANCE else
                        "DEGRADED" if delta < -geometry.GEOMETRY_TOLERANCE else "UNCHANGED")
        changes.append({"item_a":old[key]["item_a"],"item_b":old[key]["item_b"],
                        "geometry_type":old[key]["geometry_type"],"original_clearance":old[key]["distance"],
                        "new_clearance":new[key]["distance"],"clearance_change":delta,"classification":classification,
                        "removed_as_violation":new[key]["distance"]>=geometry.READABILITY_CLEARANCE})
    current_bad={key for key,item in new.items() if item["distance"]<geometry.READABILITY_CLEARANCE}
    new_items=[]
    for key in sorted(current_bad-set(old_bad)):
        new_items.append({"item_a":new[key]["item_a"],"item_b":new[key]["item_b"],
                          "geometry_type":new[key]["geometry_type"],"clearance":new[key]["distance"],
                          "deficit":new[key]["deficit"]})
    return changes,new_items


def stage(vertices,edges,original,coordinates,boundary,free,anchors,pole,scale,baseline,initial_area):
    allowed=set(free)
    if any(coordinates[node]!=original[node] for node in set(vertices)-allowed): raise AssertionError("fixed vertex moved")
    # Complete reconnection is represented by rebuilding every segment here.
    segments=[{"edge":[a,b],"start":list(coordinates[a]),"end":list(coordinates[b])} for a,b in edges]
    stale=sum(segment["start"]!=list(coordinates[segment["edge"][0]]) or segment["end"]!=list(coordinates[segment["edge"][1]]) for segment in segments)
    if stale: raise AssertionError("stale edge endpoint")
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    analysis=crowding.analyze(vertices,edges,coordinates)
    events=crowding_points.localize_events(analysis); clusters,_,_=crowding_points.cluster_events(events)
    changes,new_items=original_violation_changes(baseline,analysis)
    points=[coordinates[node] for node in boundary]
    area=abs(face_poles.polygon_area(points)); perimeter=face_poles.polygon_perimeter(points)
    pole_clearance=min(geometry.point_segment_distance(pole,points[i],points[(i+1)%len(points)]) for i in range(len(points)))
    valid=(diagnostics["proper_unrelated_edge_crossing_count"]==0 and diagnostics["coincident_vertex_pair_count"]==0 and diagnostics["vertex_on_unrelated_edge_interior_count"]==0)
    movement=[]
    for node in free:
        before=original[node]; after=coordinates[node]
        movement.append({"vertex":node,"original_coordinate":list(before),"new_coordinate":list(after),
                         "displacement_distance":math.dist(before,after),"distance_to_pole_before":math.dist(before,pole),
                         "distance_to_pole_after":math.dist(after,pole)})
    return {"s":scale,"status":"VALID ZERO-CROSSING" if valid else "INVALID / CROSSING" if diagnostics["proper_unrelated_edge_crossing_count"] else "INVALID / GEOMETRIC DEGENERACY",
            "coordinates":{node:list(coordinates[node]) for node in sorted(coordinates)},"proper_crossings":diagnostics["proper_unrelated_edge_crossing_count"],
            "proper_crossing_list":diagnostics["proper_unrelated_edge_crossings"],"coincident_vertices":diagnostics["coincident_vertex_pair_count"],
            "coincident_vertex_pairs":diagnostics["coincident_vertex_pairs"],"vertex_on_unrelated_edge_events":diagnostics["vertex_on_unrelated_edge_interior_count"],
            "vertex_on_unrelated_edge_event_list":diagnostics["vertices_on_unrelated_edge_interiors"],"minimum_node_distance":diagnostics["minimum_node_distance"],
            "electrical_edge_count":len(segments),"stale_edge_endpoint_count":stale,
            "readability_clearance":geometry.READABILITY_CLEARANCE,"total_clearance_violations":len(events),
            "unique_crowding_event_locations":unique_event_locations(events),"spatial_crowding_cluster_count":len(clusters),
            "global_item_crowding_region_count":len(analysis["regions"]),"worst_clearance":min((e["clearance"] for e in events),default=None),
            "largest_clearance_deficit":max((e["deficit"] for e in events),default=0.0),
            "original_violation_clearance_changes":changes,"newly_created_violations":new_items,
            "original_violations_removed":sum(item["removed_as_violation"] for item in changes),
            "original_violations_improved":sum(item["classification"]=="IMPROVED" for item in changes),
            "original_violations_degraded":sum(item["classification"]=="DEGRADED" for item in changes),
            "face_area":area,"face_area_ratio":area/initial_area,"face_perimeter":perimeter,
            "pole_to_nearest_boundary_distance":pole_clearance,"free_vertex_movements":movement,
            "anchor_vertices":[{"vertex":node,"coordinate":list(coordinates[node]),"coordinate_unchanged":coordinates[node]==original[node]} for node in anchors]}


def render(vertices,edges,original,boundary,free,anchors,pole,stages):
    xs=[p[0] for p in original.values()];ys=[p[1] for p in original.values()];margin=100
    minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
         f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>']
    boundary_edges={frozenset((boundary[i],boundary[(i+1)%len(boundary)])) for i in range(len(boundary))}
    for index,item in enumerate(stages):
        begin=index*2.5; timing=(f'begin="{begin}s" dur="2.5s"' if index<len(stages)-1 else f'begin="{begin}s" fill="freeze"')
        out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" {timing}/>')
        coordinates={node:tuple(point) for node,point in item["coordinates"].items()};points=' '.join(f'{coordinates[n][0]:.17g},{coordinates[n][1]:.17g}' for n in boundary)
        out.append(f'<polygon points="{points}" fill="#90caf9" fill-opacity=".18"/>')
        for a,b in edges:
            color="#1565c0" if frozenset((a,b)) in boundary_edges else "#777";sw=4 if color=="#1565c0" else 1.5
            out.append(f'<line data-edge="{a}|{b}" x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="{color}" stroke-width="{sw}"/>')
        for node in free:
            x,y=coordinates[node];out.append(f'<line x1="{pole[0]:.17g}" y1="{pole[1]:.17g}" x2="{x:.17g}" y2="{y:.17g}" stroke="#2e7d32" stroke-dasharray="6 5"/>')
        for node in sorted(vertices):
            x,y=coordinates[node];fill="#2e7d32" if node in free else "#ef6c00" if node in anchors else "#1565c0" if vertices[node]["type"]=="COMPONENT" else "#d84315"
            radius=9 if node in free or node in anchors else 6
            out.append(f'<circle data-vertex="{node}" cx="{x:.17g}" cy="{y:.17g}" r="{radius}" fill="{fill}" stroke="white" stroke-width="1.5"/>')
            out.append(f'<text x="{x+10:.17g}" y="{y-9:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
        out.extend([f'<circle cx="{pole[0]:.17g}" cy="{pole[1]:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/>',
                    f'<text x="{pole[0]+12:.17g}" y="{pole[1]-10:.17g}" font-family="sans-serif" font-size="14" font-weight="bold">RMIN1 pole</text>',
                    f'<text x="{minx+18:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="20" font-weight="bold">s={item["s"]:.2f} — {item["status"]}</text>',
                    f'<text x="{minx+18:.17g}" y="{miny+52:.17g}" font-family="sans-serif" font-size="14">crossings={item["proper_crossings"]}; clearance violations={item["total_clearance_violations"]}; green=FREE; orange=ANCHOR</text>','</g>'])
    out.append('</svg>');return '\n'.join(out)+'\n'


def run(graph_path,direct_report_path,face_report_path,full_report_path,output_dir):
    vertices,edges,original,_,_=crowding.load_inputs(graph_path,direct_report_path)
    face_doc=json.loads(face_report_path.read_text());face=next(x for x in face_doc["selected_faces"] if x["visual_location"]=="RMIN1")
    boundary=list(dict.fromkeys(face["ordered_boundary_cycle"]));pole=tuple(face["pole_coordinate"])
    classification,free,anchors=classify_boundary(boundary,edges);baseline=crowding.analyze(vertices,edges,original);initial_area=abs(face_poles.polygon_area([original[n] for n in boundary]))
    stages=[]
    for scale in SCALES: stages.append(stage(vertices,edges,original,snapshot(original,free,pole,scale),boundary,free,anchors,pole,scale,baseline,initial_area))
    full=json.loads(full_report_path.read_text()) if full_report_path.exists() else None
    full_counts={item["s"]:item["proper_crossings"] for item in full["stages"]} if full else {}
    comparison=[{"s":item["s"],"full_boundary_crossings":full_counts.get(item["s"]),"restricted_free_crossings":item["proper_crossings"]} for item in stages]
    report={"schema":"graph-relax.face-pull-rmin1-free.v1","source_coordinates":f"{direct_report_path}#final_coordinates","source_face":f"{face_report_path}#{face['face_id']}",
            "vertex_count":len(vertices),"electrical_incidence_count":len(edges),"face_id":face["face_id"],"face_boundary_cycle":boundary,"pole_coordinate":list(pole),
            "boundary_classification":classification,"free_vertices":free,"anchor_vertices":anchors,"snapshots_direct_from_original":True,
            "complete_edges_rebuilt_each_snapshot":True,"full_boundary_comparison":comparison,"stages":stages}
    output_dir.mkdir(parents=True,exist_ok=True);(output_dir/'iamp_face_pull_rmin1_free.svg').write_text(render(vertices,edges,original,boundary,free,anchors,pole,stages))
    (output_dir/'iamp_face_pull_rmin1_free_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    parser=argparse.ArgumentParser();base=Path('output/graph_first');parser.add_argument('--graph',type=Path,default=base/'iamp_graph.json');parser.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json');parser.add_argument('--face-report',type=Path,default=base/'iamp_face_poles_report.json');parser.add_argument('--full-report',type=Path,default=base/'iamp_face_pull_rmin1_report.json');parser.add_argument('--output',type=Path,default=base);a=parser.parse_args()
    report=run(a.graph,a.direct_report,a.face_report,a.full_report,a.output);print('FREE',report['free_vertices']);print('ANCHOR',report['anchor_vertices'])
    for s in report['stages']:print(s['s'],s['status'],s['proper_crossings'],s['total_clearance_violations'],s['spatial_crowding_cluster_count'],s['face_area_ratio'])

if __name__=='__main__':main()
