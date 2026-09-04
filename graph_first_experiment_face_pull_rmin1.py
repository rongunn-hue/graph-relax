#!/usr/bin/env python3
"""Diagnostic snapshots of uniform contraction of the RMIN1 bounded face."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_face_poles as face_poles
import placement_geometry as geometry

SCALES=(1.0,.9,.8,.7,.6)


def relation_key(item):
    return item["geometry_type"],item["item_a"],item["item_b"]


def compact_relation(item):
    return {"item_a":item["item_a"],"item_b":item["item_b"],
            "geometry_type":item["geometry_type"],"distance":item["distance"],
            "deficit":item["deficit"]}


def compare_crowding(baseline,current):
    old={relation_key(item):item for item in baseline["relations"]}
    new={relation_key(item):item for item in current["relations"]}
    old_bad={key for key,item in old.items() if item["distance"]<geometry.READABILITY_CLEARANCE}
    new_bad={key for key,item in new.items() if item["distance"]<geometry.READABILITY_CLEARANCE}
    def change(key):
        return {"item_a":new[key]["item_a"],"item_b":new[key]["item_b"],
                "geometry_type":new[key]["geometry_type"],
                "clearance_before":old[key]["distance"],"clearance_after":new[key]["distance"],
                "clearance_change":new[key]["distance"]-old[key]["distance"]}
    return {
        "removed_violations":[change(key) for key in sorted(old_bad-new_bad)],
        "new_violations":[change(key) for key in sorted(new_bad-old_bad)],
        "existing_violations_improved":[change(key) for key in sorted(old_bad)
                                        if new[key]["distance"]>old[key]["distance"]+geometry.GEOMETRY_TOLERANCE],
        "existing_violations_worsened":[change(key) for key in sorted(old_bad)
                                        if new[key]["distance"]<old[key]["distance"]-geometry.GEOMETRY_TOLERANCE],
        "existing_violations_unchanged":[change(key) for key in sorted(old_bad)
                                         if abs(new[key]["distance"]-old[key]["distance"])<=geometry.GEOMETRY_TOLERANCE],
    }


def snapshot(original,boundary,pole,scale):
    result=dict(original)
    for node in boundary:
        x,y=original[node]
        result[node]=(pole[0]+scale*(x-pole[0]),pole[1]+scale*(y-pole[1]))
    return result


def face_area(boundary,coordinates):
    return abs(face_poles.polygon_area([coordinates[node] for node in boundary]))


def stage_metrics(vertices,edges,coordinates,boundary,original,fixed,baseline_crowding,initial_area,scale):
    for node in fixed:
        if coordinates[node]!=original[node]: raise AssertionError(f"nonboundary vertex moved: {node}")
    # Segments are always regenerated from this authoritative coordinate map.
    rendered=[(coordinates[first],coordinates[second]) for first,second in edges]
    if any(a!=coordinates[first] or b!=coordinates[second]
           for (first,second),(a,b) in zip(edges,rendered)):
        raise AssertionError("stale edge endpoint")
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    analysis=crowding.analyze(vertices,edges,coordinates)
    area=face_area(boundary,coordinates)
    violations=analysis["violations"]
    valid=(diagnostics["proper_unrelated_edge_crossing_count"]==0 and
           diagnostics["coincident_vertex_pair_count"]==0 and
           diagnostics["vertex_on_unrelated_edge_interior_count"]==0)
    return {"s":scale,"coordinates":{node:list(coordinates[node]) for node in sorted(coordinates)},
            "moved_vertices":sorted(boundary),"reconnected_edge_count":len(edges),
            "proper_crossings":diagnostics["proper_unrelated_edge_crossing_count"],
            "proper_crossing_list":diagnostics["proper_unrelated_edge_crossings"],
            "coincident_vertex_pairs":diagnostics["coincident_vertex_pair_count"],
            "coincident_vertex_pair_list":diagnostics["coincident_vertex_pairs"],
            "vertex_on_unrelated_edge_events":diagnostics["vertex_on_unrelated_edge_interior_count"],
            "vertex_on_unrelated_edge_event_list":diagnostics["vertices_on_unrelated_edge_interiors"],
            "minimum_node_distance":diagnostics["minimum_node_distance"],
            "status":"VALID ZERO-CROSSING" if valid else "INVALID / CROSSING" if diagnostics["proper_unrelated_edge_crossing_count"] else "INVALID / GEOMETRIC DEGENERACY",
            "clearance_violation_count":len(violations),"crowded_cluster_count":len(analysis["regions"]),
            "worst_clearance":min(item["distance"] for item in violations) if violations else None,
            "largest_clearance_deficit":max((item["deficit"] for item in violations),default=0.0),
            "face_area_before":initial_area,"face_area_after":area,"face_area_ratio":area/initial_area,
            "expected_area_ratio":scale*scale,"area_ratio_error":area/initial_area-scale*scale,
            "crowding_changes_relative_to_s_1":compare_crowding(baseline_crowding,analysis)}


def draw_graph(out,vertices,edges,coordinates,boundary,pole,stage,index):
    status_color="#2e7d32" if stage["status"]=="VALID ZERO-CROSSING" else "#c62828"
    boundary_edges={frozenset((boundary[i],boundary[(i+1)%len(boundary)])) for i in range(len(boundary))}
    points=' '.join(f'{coordinates[node][0]:.17g},{coordinates[node][1]:.17g}' for node in boundary)
    out.append(f'<polygon points="{points}" fill="#90caf9" fill-opacity="0.18" stroke="none"/>')
    for first,second in edges:
        color="#1565c0" if frozenset((first,second)) in boundary_edges else "#777"
        width="4" if color=="#1565c0" else "1.5"
        out.append(f'<line data-edge="{first}|{second}" x1="{coordinates[first][0]:.17g}" y1="{coordinates[first][1]:.17g}" x2="{coordinates[second][0]:.17g}" y2="{coordinates[second][1]:.17g}" stroke="{color}" stroke-width="{width}"/>')
    for node in sorted(vertices):
        x,y=coordinates[node]; color="#1565c0" if vertices[node]["type"]=="COMPONENT" else "#d84315"
        out.append(f'<circle data-vertex="{node}" cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/>')
        out.append(f'<text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
    out.extend([f'<circle cx="{pole[0]:.17g}" cy="{pole[1]:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/>',
                f'<text x="{pole[0]+13:.17g}" y="{pole[1]-10:.17g}" font-family="sans-serif" font-size="14" font-weight="bold" fill="#6a1b9a">RMIN1 face pole</text>',
                f'<text x="-350" y="-475" font-family="sans-serif" font-size="22" font-weight="bold">s={stage["s"]:.2f} — <tspan fill="{status_color}">{stage["status"]}</tspan></text>',
                f'<text x="-350" y="-445" font-family="sans-serif" font-size="14">crossings={stage["proper_crossings"]}; clearance violations={stage["clearance_violation_count"]}; minimum node distance={stage["minimum_node_distance"]:.6f}</text>'])


def render(vertices,edges,original,boundary,pole,stages):
    xs=[p[0] for p in original.values()];ys=[p[1] for p in original.values()]
    margin=100;minx=min(xs)-margin;miny=min(ys)-margin;width=max(xs)-min(xs)+2*margin;height=max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
         f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',
         '<text x="-350" y="-505" font-family="sans-serif" font-size="18">RMIN1 face contraction — independent direct snapshots</text>']
    for index,stage in enumerate(stages):
        begin=index*2.5
        if index<len(stages)-1:
            out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" begin="{begin}s" dur="2.5s"/>')
        else:
            out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" begin="{begin}s" fill="freeze"/>')
        coordinates={node:tuple(point) for node,point in stage["coordinates"].items()}
        draw_graph(out,vertices,edges,coordinates,boundary,pole,stage,index);out.append('</g>')
    out.append('</svg>');return '\n'.join(out)+'\n'


def run(graph_path,direct_report_path,face_report_path,output_dir):
    vertices,edges,original,_,_=crowding.load_inputs(graph_path,direct_report_path)
    face_report=json.loads(face_report_path.read_text())
    face=next(item for item in face_report["selected_faces"] if item["visual_location"]=="RMIN1")
    boundary=list(dict.fromkeys(face["ordered_boundary_cycle"]));pole=tuple(face["pole_coordinate"])
    fixed=set(vertices)-set(boundary);initial_area=face_area(boundary,original)
    baseline=crowding.analyze(vertices,edges,original)
    stages=[]
    for scale in SCALES:
        coordinates=snapshot(original,boundary,pole,scale)
        stages.append(stage_metrics(vertices,edges,coordinates,boundary,original,fixed,baseline,initial_area,scale))
    report={"schema":"graph-relax.face-pull-rmin1.v1","source_coordinates":f"{direct_report_path}#final_coordinates",
            "source_face":f"{face_report_path}#{face['face_id']}","vertex_count":len(vertices),
            "electrical_incidence_count":len(edges),"face_id":face["face_id"],"face_boundary_vertices":boundary,
            "fixed_vertices":sorted(fixed),"pole_coordinate":list(pole),"snapshots_direct_from_original":True,
            "all_edges_rebuilt_from_authoritative_endpoints":True,"stages":stages}
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/'iamp_face_pull_rmin1.svg').write_text(render(vertices,edges,original,boundary,pole,stages))
    (output_dir/'iamp_face_pull_rmin1_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser();base=Path('output/graph_first')
    parser.add_argument('--graph',type=Path,default=base/'iamp_graph.json');parser.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json')
    parser.add_argument('--face-report',type=Path,default=base/'iamp_face_poles_report.json');parser.add_argument('--output',type=Path,default=base);args=parser.parse_args()
    report=run(args.graph,args.direct_report,args.face_report,args.output)
    for stage in report['stages']:
        print(stage['s'],stage['status'],'crossings',stage['proper_crossings'],'min-node',stage['minimum_node_distance'],'violations',stage['clearance_violation_count'],'clusters',stage['crowded_cluster_count'],'area-ratio',stage['face_area_ratio'])

if __name__=='__main__':main()
