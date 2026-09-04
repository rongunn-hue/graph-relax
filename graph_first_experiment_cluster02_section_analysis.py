#!/usr/bin/env python3
"""Measurement-only separation-axis analysis of the eight-event cluster 02."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import placement_geometry as geometry


def edge_parameter(point, first, second):
    dx, dy = second[0]-first[0], second[1]-first[1]
    denominator = dx*dx+dy*dy
    if denominator == 0.0:
        return 0.0
    return max(0.0, min(1.0, ((point[0]-first[0])*dx+(point[1]-first[1])*dy)/denominator))


def choose_cluster(point_report):
    candidates = [cluster for cluster in point_report["clusters"]
                  if cluster["violation_count"] == 8]
    if len(candidates) != 1:
        raise AssertionError("expected exactly one eight-violation localized cluster")
    ids = set(candidates[0]["member_violation_ids"])
    return candidates[0], [event for event in point_report["events"]
                           if event["violation_id"] in ids]


def distribute(item, point, vector, weight, coordinates, accumulated, records):
    weighted = (vector[0]*weight, vector[1]*weight)
    if item.startswith("vertex:"):
        node = item[len("vertex:"):]
        accumulated[node][0] += weighted[0]; accumulated[node][1] += weighted[1]
        records.append({"vertex": node, "barycentric_weight": 1.0,
                        "weighted_vector": list(weighted)})
        return
    first, second = item[len("edge:"):].split("|")
    parameter = edge_parameter(point, coordinates[first], coordinates[second])
    for node, barycentric in ((first, 1.0-parameter), (second, parameter)):
        contribution = (weighted[0]*barycentric, weighted[1]*barycentric)
        accumulated[node][0] += contribution[0]; accumulated[node][1] += contribution[1]
        records.append({"vertex": node, "barycentric_weight": barycentric,
                        "segment_parameter": parameter,
                        "weighted_vector": list(contribution)})


def analyze(point_report, coordinates):
    cluster, source_events = choose_cluster(point_report)
    events = []
    doubled_x = doubled_y = total_weight = 0.0
    for event in source_events:
        first, second = event["pA"], event["pB"]
        dx, dy = first[0]-second[0], first[1]-second[1]
        distance = math.hypot(dx, dy)
        if distance <= geometry.GEOMETRY_TOLERANCE:
            raise AssertionError("violation has undefined separation axis")
        unit = (dx/distance, dy/distance)
        angle = math.atan2(unit[1], unit[0]) % math.pi
        weight = event["deficit"]
        doubled_x += weight*math.cos(2*angle)
        doubled_y += weight*math.sin(2*angle)
        total_weight += weight
        events.append({"violation_id": event["violation_id"],
                       "item_a": event["item_a"], "item_b": event["item_b"],
                       "pA": list(first), "pB": list(second),
                       "clearance": event["clearance"], "deficit": weight,
                       "nA": list(unit), "nB": [-unit[0], -unit[1]],
                       "separation_axis_radians": angle,
                       "separation_axis_degrees_modulo_180": math.degrees(angle),
                       "crowding_point": event["crowding_point"]})
    axis = (0.5*math.atan2(doubled_y, doubled_x)) % math.pi
    coherence = math.hypot(doubled_x, doubled_y)/total_weight
    normal = (math.cos(axis), math.sin(axis))
    centroid = cluster["crowding_point_centroid"]
    participating = set()
    for event in events:
        for item in (event["item_a"], event["item_b"]):
            if item.startswith("vertex:"):
                participating.add(item[len("vertex:"):])
            else:
                participating.update(item[len("edge:"):].split("|"))
    accumulated = {node: [0.0, 0.0] for node in sorted(participating)}
    contribution_records = []
    for event in events:
        distribute(event["item_a"], tuple(event["pA"]), tuple(event["nA"]),
                   event["deficit"], coordinates, accumulated, contribution_records)
        distribute(event["item_b"], tuple(event["pB"]), tuple(event["nB"]),
                   event["deficit"], coordinates, accumulated, contribution_records)
    vertices = []
    side_groups = {"SIDE_A": [], "SIDE_B": [], "ON_AXIS": []}
    for node in sorted(participating):
        point = coordinates[node]
        projection = (point[0]-centroid[0])*normal[0]+(point[1]-centroid[1])*normal[1]
        side = ("SIDE_A" if projection > geometry.GEOMETRY_TOLERANCE else
                "SIDE_B" if projection < -geometry.GEOMETRY_TOLERANCE else "ON_AXIS")
        vector = accumulated[node]
        magnitude = math.hypot(*vector)
        vector_angle = math.degrees(math.atan2(vector[1], vector[0])) % 360 if magnitude else None
        axial_component = vector[0]*normal[0]+vector[1]*normal[1]
        record = {"vertex": node, "side": side, "side_projection": projection,
                  "accumulated_vector": vector, "accumulated_vector_magnitude": magnitude,
                  "accumulated_vector_angle_degrees": vector_angle,
                  "accumulated_vector_axis_component": axial_component}
        vertices.append(record); side_groups[side].append(node)
    participating_edges = sorted({tuple(item[len("edge:"):].split("|"))
                                  for event in events for item in (event["item_a"], event["item_b"])
                                  if item.startswith("edge:")})
    edge_projections = []
    for first, second in participating_edges:
        edge_projections.append({"edge": [first, second],
                                 "endpoint_projections": {
                                     first: (coordinates[first][0]-centroid[0])*normal[0]+(coordinates[first][1]-centroid[1])*normal[1],
                                     second: (coordinates[second][0]-centroid[0])*normal[0]+(coordinates[second][1]-centroid[1])*normal[1]},
                                 "straddles_axis": (((coordinates[first][0]-centroid[0])*normal[0]+(coordinates[first][1]-centroid[1])*normal[1]) *
                                                    ((coordinates[second][0]-centroid[0])*normal[0]+(coordinates[second][1]-centroid[1])*normal[1]) < 0)})
    positive = sum(record["accumulated_vector_axis_component"] > geometry.GEOMETRY_TOLERANCE
                   for record in vertices)
    negative = sum(record["accumulated_vector_axis_component"] < -geometry.GEOMETRY_TOLERANCE
                   for record in vertices)
    zero = len(vertices)-positive-negative
    vector_direction_groups = {
        "POSITIVE_ALONG_AXIS": [record["vertex"] for record in vertices
                                if record["accumulated_vector_axis_component"] > geometry.GEOMETRY_TOLERANCE],
        "NEGATIVE_ALONG_AXIS": [record["vertex"] for record in vertices
                                if record["accumulated_vector_axis_component"] < -geometry.GEOMETRY_TOLERANCE],
        "NO_NET_TENDENCY": [record["vertex"] for record in vertices
                            if abs(record["accumulated_vector_axis_component"]) <= geometry.GEOMETRY_TOLERANCE],
    }
    # A zero-contribution participant is neutral evidence, not evidence against
    # the two nonzero, oppositely directed groups.  This classification uses no
    # fitted threshold: only the existing numerical geometry tolerance.
    structure = ("TWO_COHERENT_OPPOSING_SECTIONS" if positive and negative
                 else "ONE_COHERENT_MOVABLE_SECTION" if bool(positive) ^ bool(negative)
                 else "NO_SIMPLE_SECTION_STRUCTURE")
    return {"source_cluster": cluster["cluster_id"], "violation_count": len(events),
            "cluster_centroid": centroid, "events": events,
            "dominant_separation_axis_radians": axis,
            "dominant_separation_axis_degrees": math.degrees(axis),
            "axial_coherence": coherence, "dominant_axis_unit_vector": list(normal),
            "side_groups": side_groups, "vertices": vertices,
            "participating_edges": edge_projections,
            "barycentric_contributions": contribution_records,
            "accumulated_axis_sign_counts": {"positive": positive, "negative": negative, "zero": zero},
            "accumulated_vector_direction_groups": vector_direction_groups,
            "numerical_structure_classification": structure}


def render_svg(vertices, edges, coordinates, analysis):
    xs, ys = [p[0] for p in coordinates.values()], [p[1] for p in coordinates.values()]
    margin=100.; min_x,min_y=min(xs)-margin,min(ys)-margin
    width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.17g} {min_y:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
           f'<rect x="{min_x:.17g}" y="{min_y:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',
           f'<text x="{min_x+16:.17g}" y="{min_y+28:.17g}" font-family="sans-serif" font-size="18">Cluster 02 separation-axis analysis — diagnostic only</text>']
    for first,second in edges:
        parts.append(f'<line x1="{coordinates[first][0]:.17g}" y1="{coordinates[first][1]:.17g}" x2="{coordinates[second][0]:.17g}" y2="{coordinates[second][1]:.17g}" stroke="#aaa" stroke-width="1.5"/>')
    for node in sorted(vertices):
        x,y=coordinates[node]; color="#1565c0" if vertices[node]["type"]=="COMPONENT" else "#d84315"
        parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/>')
        parts.append(f'<text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
    for event in analysis["events"]:
        a,b,c=event["pA"],event["pB"],event["crowding_point"]
        parts.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#e91e63" stroke-width="2"/>')
        parts.append(f'<circle cx="{c[0]:.17g}" cy="{c[1]:.17g}" r="4" fill="#e91e63"/>')
        parts.append(f'<text x="{c[0]+6:.17g}" y="{c[1]-6:.17g}" font-family="sans-serif" font-size="9" fill="#ad1457">{event["violation_id"]}</text>')
    cx,cy=analysis["cluster_centroid"]; nx_,ny_=analysis["dominant_axis_unit_vector"]
    parts.append(f'<line x1="{cx-90*nx_:.17g}" y1="{cy-90*ny_:.17g}" x2="{cx+90*nx_:.17g}" y2="{cy+90*ny_:.17g}" stroke="#6a1b9a" stroke-width="3" stroke-dasharray="8 5"/>')
    parts.append(f'<circle cx="{cx:.17g}" cy="{cy:.17g}" r="7" fill="none" stroke="#6a1b9a" stroke-width="2"/>')
    max_magnitude=max(record["accumulated_vector_magnitude"] for record in analysis["vertices"])
    scale=45/max_magnitude if max_magnitude else 0
    by_vertex={record["vertex"]:record for record in analysis["vertices"]}
    for node,record in by_vertex.items():
        x,y=coordinates[node]; vx,vy=record["accumulated_vector"]
        color="#00897b" if record["side"]=="SIDE_A" else "#f57c00" if record["side"]=="SIDE_B" else "#555"
        parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="11" fill="none" stroke="{color}" stroke-width="2"/>')
        parts.append(f'<line x1="{x:.17g}" y1="{y:.17g}" x2="{x+scale*vx:.17g}" y2="{y+scale*vy:.17g}" stroke="{color}" stroke-width="3"/>')
    parts.append('</svg>'); return "\n".join(parts)+"\n"


def run(graph_path, direct_report_path, points_report_path, output_dir):
    vertices, edges, coordinates, _, _=crowding.load_inputs(graph_path,direct_report_path)
    original=dict(coordinates); point_report=json.loads(points_report_path.read_text())
    result=analyze(point_report,coordinates)
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    if coordinates!=original or diagnostics["proper_unrelated_edge_crossing_count"]!=0:
        raise AssertionError("analysis changed graph geometry")
    report={"schema":"graph-relax.cluster02-section-analysis.v1",
            "source_coordinates":f"{direct_report_path.as_posix()}#final_coordinates",
            "vertex_count":len(vertices),"electrical_incidence_count":len(edges),
            "proper_crossing_count":diagnostics["proper_unrelated_edge_crossing_count"],
            "coordinates_unchanged":coordinates==original,**result}
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/"iamp_cluster02_section_analysis.svg").write_text(render_svg(vertices,edges,coordinates,result),encoding="utf-8")
    (output_dir/"iamp_cluster02_section_analysis_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main():
    parser=argparse.ArgumentParser();base=Path("output/graph_first")
    parser.add_argument("--graph",type=Path,default=base/"iamp_graph.json")
    parser.add_argument("--direct-report",type=Path,default=base/"iamp_nearness_rotation_direct_report.json")
    parser.add_argument("--points-report",type=Path,default=base/"iamp_crowding_points_report.json")
    parser.add_argument("--output",type=Path,default=base);args=parser.parse_args()
    report=run(args.graph,args.direct_report,args.points_report,args.output)
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=="__main__":main()
