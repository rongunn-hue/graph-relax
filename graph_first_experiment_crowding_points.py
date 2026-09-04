#!/usr/bin/env python3
"""Localize and spatially group existing readability-clearance events."""

from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import networkx as nx

import graph_first_experiment_crowding_analysis as crowding
import placement_geometry as geometry


CROWDING_CLUSTER_DISTANCE = 60.0
COLORS = ["#d81b60", "#8e24aa", "#3949ab", "#00897b", "#7cb342",
          "#f4511e", "#6d4c41", "#039be5", "#c0ca33", "#5e35b1"]


def localize_events(analysis):
    ordered = sorted(analysis["violations"],
                     key=lambda item: (-item["deficit"], item["item_a"], item["item_b"]))
    events = []
    for index, violation in enumerate(ordered, 1):
        first, second = (tuple(violation["closest_points"][0]),
                         tuple(violation["closest_points"][1]))
        event = {"violation_id": f"V{index:02d}",
                 "item_a": violation["item_a"], "item_b": violation["item_b"],
                 "violation_type": violation["geometry_type"],
                 "clearance": violation["distance"], "deficit": violation["deficit"],
                 "pA": list(first), "pB": list(second),
                 "crowding_point": [(first[0]+second[0])/2, (first[1]+second[1])/2]}
        events.append(event)
    return events


def participating_items(events):
    vertices, edges = set(), set()
    for event in events:
        for item in (event["item_a"], event["item_b"]):
            if item.startswith("vertex:"):
                vertices.add(item[len("vertex:"):])
            else:
                edges.add(tuple(item[len("edge:"):].split("|")))
    return sorted(vertices), [list(edge) for edge in sorted(edges)]


def cluster_events(events):
    proximity = nx.Graph()
    proximity.add_nodes_from(event["violation_id"] for event in events)
    by_id = {event["violation_id"]: event for event in events}
    pair_distances = []
    for first, second in combinations(events, 2):
        distance = math.dist(first["crowding_point"], second["crowding_point"])
        pair = {"event_a": first["violation_id"], "event_b": second["violation_id"],
                "distance": distance,
                "within_cluster_distance": distance <= CROWDING_CLUSTER_DISTANCE}
        pair_distances.append(pair)
        if pair["within_cluster_distance"]:
            proximity.add_edge(pair["event_a"], pair["event_b"], distance=distance)
    components = sorted((sorted(component) for component in nx.connected_components(proximity)),
                        key=lambda members: tuple(members))
    clusters, event_cluster = [], {}
    for index, members in enumerate(components, 1):
        cluster_id = f"cluster:{index:02d}"
        member_events = [by_id[item] for item in members]
        points = [event["crowding_point"] for event in member_events]
        xs, ys = [point[0] for point in points], [point[1] for point in points]
        diameter = max((math.dist(a, b) for a, b in combinations(points, 2)), default=0.0)
        vertices, edges = participating_items(member_events)
        cluster = {"cluster_id": cluster_id, "member_violation_ids": members,
                   "violation_count": len(members),
                   "minimum_clearance": min(event["clearance"] for event in member_events),
                   "maximum_deficit": max(event["deficit"] for event in member_events),
                   "crowding_point_centroid": [sum(xs)/len(xs), sum(ys)/len(ys)],
                   "crowding_point_bounding_box": {"min_x": min(xs), "min_y": min(ys),
                                                    "max_x": max(xs), "max_y": max(ys)},
                   "crowding_point_diameter": diameter,
                   "participating_vertices": vertices, "participating_edges": edges}
        clusters.append(cluster)
        for member in members:
            event_cluster[member] = cluster_id
    clusters.sort(key=lambda item: (-item["maximum_deficit"], -item["violation_count"],
                                    item["cluster_id"]))
    for event in events:
        event["cluster_id"] = event_cluster[event["violation_id"]]
    pair_distances.sort(key=lambda item: (item["distance"], item["event_a"], item["event_b"]))
    nearest_neighbors = []
    for event in events:
        related = [item for item in pair_distances
                   if event["violation_id"] in {item["event_a"], item["event_b"]}]
        nearest = min(related, key=lambda item: (item["distance"], item["event_a"], item["event_b"]))
        nearest_neighbors.append({"event": event["violation_id"],
                                  "neighbor": (nearest["event_b"] if nearest["event_a"] == event["violation_id"]
                                               else nearest["event_a"]),
                                  "distance": nearest["distance"],
                                  "within_cluster_distance": nearest["within_cluster_distance"]})
    return clusters, pair_distances, nearest_neighbors


def render_svg(vertices, edges, coordinates, events, clusters):
    xs, ys = [point[0] for point in coordinates.values()], [point[1] for point in coordinates.values()]
    margin = 100.0; min_x, min_y = min(xs)-margin, min(ys)-margin
    width, height = max(xs)-min(xs)+2*margin, max(ys)-min(ys)+2*margin
    color_by_cluster = {cluster["cluster_id"]: COLORS[index % len(COLORS)]
                        for index, cluster in enumerate(clusters)}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.17g} {min_y:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
             f'<rect x="{min_x:.17g}" y="{min_y:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',
             f'<text x="{min_x+16:.17g}" y="{min_y+28:.17g}" font-family="sans-serif" font-size="18">IAMP localized crowding-event points; cluster distance 60.0</text>']
    for first, second in edges:
        parts.append(f'<line x1="{coordinates[first][0]:.17g}" y1="{coordinates[first][1]:.17g}" '
                     f'x2="{coordinates[second][0]:.17g}" y2="{coordinates[second][1]:.17g}" '
                     'stroke="#777" stroke-width="1.5"/>')
    for node in sorted(vertices):
        x, y = coordinates[node]
        color = "#1565c0" if vertices[node]["type"] == "COMPONENT" else "#d84315"
        parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/>')
        parts.append(f'<text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
    for event in events:
        color = color_by_cluster[event["cluster_id"]]
        first, second, point = event["pA"], event["pB"], event["crowding_point"]
        parts.append(f'<line x1="{first[0]:.17g}" y1="{first[1]:.17g}" '
                     f'x2="{second[0]:.17g}" y2="{second[1]:.17g}" '
                     f'stroke="{color}" stroke-width="1.2" opacity="0.65"/>')
        parts.append(f'<circle cx="{point[0]:.17g}" cy="{point[1]:.17g}" r="5" '
                     f'fill="{color}" stroke="white" stroke-width="1"/>')
        parts.append(f'<text x="{point[0]+7:.17g}" y="{point[1]-7:.17g}" '
                     f'font-family="sans-serif" font-size="10" fill="{color}">{event["violation_id"]}</text>')
    for cluster in clusters:
        x, y = cluster["crowding_point_centroid"]
        color = color_by_cluster[cluster["cluster_id"]]
        parts.append(f'<text x="{x+12:.17g}" y="{y+22:.17g}" font-family="sans-serif" '
                     f'font-size="14" font-weight="bold" fill="{color}">{cluster["cluster_id"]}</text>')
    parts.append('</svg>')
    return "\n".join(parts)+"\n"


def run(graph_path, direct_report_path, output_dir):
    vertices, edges, coordinates, _, _ = crowding.load_inputs(graph_path, direct_report_path)
    original = dict(coordinates)
    analysis = crowding.analyze(vertices, edges, coordinates)
    events = localize_events(analysis)
    clusters, pair_distances, nearest_neighbors = cluster_events(events)
    diagnostics = geometry.graph_geometry_diagnostics(coordinates, edges)
    if coordinates != original or diagnostics["proper_unrelated_edge_crossing_count"] != 0:
        raise AssertionError("crowding-point measurement changed graph geometry")
    report = {"schema": "graph-relax.localized-crowding-events.v1",
              "source_graph": graph_path.as_posix(),
              "source_coordinates": f"{direct_report_path.as_posix()}#final_coordinates",
              "readability_clearance": geometry.READABILITY_CLEARANCE,
              "crowding_cluster_distance": CROWDING_CLUSTER_DISTANCE,
              "vertex_count": len(vertices), "electrical_incidence_count": len(edges),
              "proper_crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
              "coordinates_unchanged": coordinates == original,
              "total_clearance_violations": len(events),
              "total_crowding_event_points": len(events),
              "spatial_cluster_count": len(clusters),
              "events": events, "clusters": clusters,
              "event_pair_distances": pair_distances,
              "proximity_graph_edges": [item for item in pair_distances
                                          if item["within_cluster_distance"]],
              "nearest_event_neighbors": nearest_neighbors}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"iamp_crowding_points.svg").write_text(
        render_svg(vertices, edges, coordinates, events, clusters), encoding="utf-8")
    (output_dir/"iamp_crowding_points_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser()
    base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base/"iamp_graph.json")
    parser.add_argument("--direct-report", type=Path,
                        default=base/"iamp_nearness_rotation_direct_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    report = run(args.graph, args.direct_report, args.output)
    print("violations", report["total_clearance_violations"], "points",
          report["total_crowding_event_points"], "clusters", report["spatial_cluster_count"])
    for cluster in report["clusters"]:
        print(json.dumps(cluster, sort_keys=True))
    print("nearest event-point neighbors")
    for item in report["nearest_event_neighbors"]:
        print(item["event"], item["neighbor"], item["distance"])


if __name__ == "__main__":
    main()
