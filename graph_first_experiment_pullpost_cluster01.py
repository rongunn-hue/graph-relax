#!/usr/bin/env python3
"""One direct pull-post move for localized crowding cluster:01."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import placement_geometry as geometry


POST_DISTANCE = 10000.0
TARGET_CLUSTER = "cluster:01"


def violation_signature(event):
    return event["item_a"], event["item_b"], event["violation_type"]


def direct_pullpost_move(vertices, edges, coordinates, point_report):
    cluster = next(item for item in point_report["clusters"]
                   if item["cluster_id"] == TARGET_CLUSTER)
    if cluster["member_violation_ids"] != ["V01"]:
        raise AssertionError("cluster:01 is no longer the isolated V01 event")
    event = next(item for item in point_report["events"] if item["violation_id"] == "V01")
    if event["item_a"] != "vertex:component:TPB" or \
            event["item_b"] != "edge:component:R3|net:NODE_B":
        raise AssertionError("cluster:01 does not identify TPB versus R3–NODE_B")
    moving = "component:TPB"
    old = coordinates[moving]
    closest_vertex_point = tuple(event["pA"])
    closest_edge_point = tuple(event["pB"])
    if math.dist(old, closest_vertex_point) > geometry.GEOMETRY_TOLERANCE:
        raise AssertionError("saved closest vertex point is not TPB")
    dx, dy = old[0]-closest_edge_point[0], old[1]-closest_edge_point[1]
    current_clearance = math.hypot(dx, dy)
    if current_clearance <= geometry.GEOMETRY_TOLERANCE:
        raise ValueError("outward direction is undefined at zero clearance")
    direction = (dx/current_clearance, dy/current_clearance)
    requested_displacement = geometry.READABILITY_CLEARANCE-current_clearance
    if requested_displacement < 0:
        raise ValueError("target event already satisfies readability clearance")
    # Algebraically identical to V + (target-current)*d, but anchored at the
    # closest point so the requested point-to-segment distance is represented
    # without subtractive cancellation.
    numerical_target = geometry.READABILITY_CLEARANCE + geometry.GEOMETRY_TOLERANCE
    new = (closest_edge_point[0]+numerical_target*direction[0],
           closest_edge_point[1]+numerical_target*direction[1])
    displacement = math.dist(old, new)
    post = (old[0]+POST_DISTANCE*direction[0], old[1]+POST_DISTANCE*direction[1])
    edge = ("component:R3", "net:NODE_B")
    resulting_clearance = geometry.point_segment_distance(
        new, coordinates[edge[0]], coordinates[edge[1]])
    if resulting_clearance + geometry.GEOMETRY_TOLERANCE < geometry.READABILITY_CLEARANCE:
        raise AssertionError("direct pull-post calculation did not reach target clearance")
    planted = dict(coordinates)
    planted[moving] = new
    rebuilt = geometry.rebuild_straight_edge_segments(planted, edges)
    reconnection = geometry.validate_rebuilt_edge_segments(planted, edges, rebuilt)
    if not reconnection["valid"]:
        raise AssertionError("TPB incident edges were not reconnected")
    diagnostics = geometry.graph_geometry_diagnostics(planted, edges)
    valid = (diagnostics["proper_unrelated_edge_crossing_count"] == 0 and
             diagnostics["coincident_vertex_pair_count"] == 0 and
             diagnostics["vertex_on_unrelated_edge_interior_count"] == 0 and
             diagnostics["minimum_node_distance"] + geometry.GEOMETRY_TOLERANCE
             >= geometry.MIN_NODE_DISTANCE)
    if not valid:
        return None, {"status": "REJECTED_GRAPH_VALIDITY", "diagnostics": diagnostics}
    incident = [item for item in rebuilt if moving in item["edge"]]
    return planted, {"status": "ACCEPTED", "moving_vertex": moving,
                     "old_coordinate": list(old), "new_coordinate": list(new),
                     "closest_edge_point": list(closest_edge_point),
                     "movement_direction": list(direction),
                     "movement_angle_radians": math.atan2(direction[1], direction[0]),
                     "movement_angle_degrees": math.degrees(math.atan2(direction[1], direction[0])),
                     "pull_post_coordinate": list(post), "post_distance": POST_DISTANCE,
                     "displacement": displacement,
                     "exact_geometric_displacement_to_target": requested_displacement,
                     "clearance_before": current_clearance,
                     "clearance_after": resulting_clearance,
                     "numerical_legal_side_tolerance": geometry.GEOMETRY_TOLERANCE,
                     "reconnection_validation": reconnection,
                     "incident_edges_regenerated": incident,
                     "diagnostics": diagnostics}


def render_svg(vertices, edges, before, after, movement):
    xs = [point[0] for point in before.values()] + [point[0] for point in after.values()]
    ys = [point[1] for point in before.values()] + [point[1] for point in after.values()]
    margin = 100.0; min_x, min_y = min(xs)-margin, min(ys)-margin
    width, height = max(xs)-min(xs)+2*margin, max(ys)-min(ys)+2*margin
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.17g} {min_y:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
             '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#00897b"/></marker></defs>',
             f'<rect x="{min_x:.17g}" y="{min_y:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',
             f'<text x="{min_x+16:.17g}" y="{min_y+28:.17g}" font-family="sans-serif" font-size="18">Single pull-post de-crowding: cluster:01</text>']
    for first, second in edges:
        parts.append(f'<line x1="{before[first][0]:.17g}" y1="{before[first][1]:.17g}" '
                     f'x2="{before[second][0]:.17g}" y2="{before[second][1]:.17g}" '
                     'stroke="#bdbdbd" stroke-width="1" opacity="0.35"/>')
    for first, second in edges:
        parts.append(f'<line x1="{after[first][0]:.17g}" y1="{after[first][1]:.17g}" '
                     f'x2="{after[second][0]:.17g}" y2="{after[second][1]:.17g}" '
                     'stroke="#616161" stroke-width="1.5"/>')
    for node in sorted(vertices):
        x, y = after[node]
        color = "#1565c0" if vertices[node]["type"] == "COMPONENT" else "#d84315"
        parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/>')
        parts.append(f'<text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
    old, new = movement["old_coordinate"], movement["new_coordinate"]
    closest = movement["closest_edge_point"]
    direction = movement["movement_direction"]
    display_post = (new[0]+90*direction[0], new[1]+90*direction[1])
    parts.append(f'<circle cx="{old[0]:.17g}" cy="{old[1]:.17g}" r="8" fill="none" stroke="#e53935" stroke-width="2" stroke-dasharray="3 3"/>')
    parts.append(f'<line x1="{old[0]:.17g}" y1="{old[1]:.17g}" x2="{new[0]:.17g}" y2="{new[1]:.17g}" stroke="#00897b" stroke-width="3" marker-end="url(#arrow)"/>')
    parts.append(f'<line x1="{new[0]:.17g}" y1="{new[1]:.17g}" x2="{display_post[0]:.17g}" y2="{display_post[1]:.17g}" stroke="#00897b" stroke-width="1.5" stroke-dasharray="7 5" marker-end="url(#arrow)"/>')
    parts.append(f'<line x1="{old[0]:.17g}" y1="{old[1]:.17g}" x2="{closest[0]:.17g}" y2="{closest[1]:.17g}" stroke="#e53935" stroke-width="2"/>')
    final_closest = crowding.project_to_segment(tuple(new), after["component:R3"], after["net:NODE_B"])
    parts.append(f'<line x1="{new[0]:.17g}" y1="{new[1]:.17g}" x2="{final_closest[0]:.17g}" y2="{final_closest[1]:.17g}" stroke="#43a047" stroke-width="3"/>')
    parts.append(f'<circle cx="{new[0]:.17g}" cy="{new[1]:.17g}" r="10" fill="none" stroke="#43a047" stroke-width="3"/>')
    parts.append(f'<text x="{new[0]+12:.17g}" y="{new[1]+18:.17g}" font-family="sans-serif" font-size="12" fill="#2e7d32">final clearance 6.0; post at ({movement["pull_post_coordinate"][0]:.3f}, {movement["pull_post_coordinate"][1]:.3f})</text>')
    parts.append('</svg>')
    return "\n".join(parts)+"\n"


def run(graph_path, direct_report_path, crowding_points_report_path, output_dir):
    vertices, edges, coordinates, _, _ = crowding.load_inputs(graph_path, direct_report_path)
    point_report = json.loads(crowding_points_report_path.read_text(encoding="utf-8"))
    before_analysis = crowding.analyze(vertices, edges, coordinates)
    before_events = crowding_points.localize_events(before_analysis)
    before_clusters, _, _ = crowding_points.cluster_events(before_events)
    after, movement = direct_pullpost_move(vertices, edges, coordinates, point_report)
    if after is None:
        report = {"status": movement["status"], "movement": movement}
    else:
        after_analysis = crowding.analyze(vertices, edges, after)
        after_events = crowding_points.localize_events(after_analysis)
        after_clusters, _, _ = crowding_points.cluster_events(after_events)
        before_signatures = {violation_signature(event) for event in before_events}
        after_signatures = {violation_signature(event) for event in after_events}
        report = {"schema": "graph-relax.pullpost-cluster01.v1", "status": "ACCEPTED",
                  "source_coordinates": f"{direct_report_path.as_posix()}#final_coordinates",
                  "target_cluster": TARGET_CLUSTER,
                  "readability_clearance": geometry.READABILITY_CLEARANCE,
                  "post_distance": POST_DISTANCE,
                  "vertex_count": len(vertices), "electrical_incidence_count": len(edges),
                  "movement": movement,
                  "before": {"clearance_violation_count": len(before_events),
                             "crowded_cluster_count": len(before_clusters),
                             "proper_crossings": geometry.graph_geometry_diagnostics(coordinates, edges)["proper_unrelated_edge_crossing_count"],
                             "global_minimum_node_distance": geometry.node_spacing_diagnostics(coordinates)["minimum_node_distance"]},
                  "after": {"clearance_violation_count": len(after_events),
                            "crowded_cluster_count": len(after_clusters),
                            "proper_crossings": movement["diagnostics"]["proper_unrelated_edge_crossing_count"],
                            "global_minimum_node_distance": movement["diagnostics"]["minimum_node_distance"]},
                  "removed_violation_signatures": [list(item) for item in sorted(before_signatures-after_signatures)],
                  "new_violation_signatures": [list(item) for item in sorted(after_signatures-before_signatures)],
                  "after_crowding_events": after_events, "after_crowding_clusters": after_clusters,
                  "final_coordinates": {node: list(after[node]) for node in sorted(after)}}
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir/"iamp_pullpost_cluster01.svg").write_text(
            render_svg(vertices, edges, coordinates, after, movement), encoding="utf-8")
    (output_dir/"iamp_pullpost_cluster01_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(); base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base/"iamp_graph.json")
    parser.add_argument("--direct-report", type=Path, default=base/"iamp_nearness_rotation_direct_report.json")
    parser.add_argument("--crowding-points-report", type=Path, default=base/"iamp_crowding_points_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    print(json.dumps(run(args.graph, args.direct_report, args.crowding_points_report, args.output),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
