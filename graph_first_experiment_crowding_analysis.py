#!/usr/bin/env python3
"""Measurement-only global readability-clearance analysis."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import networkx as nx

import placement_geometry as geometry


COLORS = ["#d81b60", "#8e24aa", "#3949ab", "#00897b", "#7cb342",
          "#f4511e", "#6d4c41", "#039be5", "#c0ca33", "#5e35b1"]


def load_inputs(graph_path, direct_report_path):
    graph_document = json.loads(graph_path.read_text(encoding="utf-8"))
    direct_report = json.loads(direct_report_path.read_text(encoding="utf-8"))
    vertices = {item["id"]: item for item in graph_document["vertices"]}
    edges = sorted((min(item["source"], item["target"]),
                    max(item["source"], item["target"]))
                   for item in graph_document["edges"])
    coordinates = {node: tuple(point)
                   for node, point in direct_report["final_coordinates"].items()}
    if set(vertices) != set(coordinates) or len(vertices) != 34 or len(edges) != 44:
        raise AssertionError("authoritative graph/final-coordinate mismatch")
    if direct_report["final_crossing_count"] != 0:
        raise AssertionError("input is not the final zero-crossing direct drawing")
    return vertices, edges, coordinates, graph_document, direct_report


def edge_id(edge):
    return f"edge:{edge[0]}|{edge[1]}"


def vertex_id(node):
    return f"vertex:{node}"


def project_to_segment(point, first, second):
    dx, dy = second[0]-first[0], second[1]-first[1]
    denominator = dx*dx+dy*dy
    if denominator == 0.0:
        return first
    fraction = ((point[0]-first[0])*dx+(point[1]-first[1])*dy)/denominator
    fraction = max(0.0, min(1.0, fraction))
    return first[0]+fraction*dx, first[1]+fraction*dy


def closest_segment_points(a, b, c, d):
    candidates = [(a, project_to_segment(a, c, d)),
                  (b, project_to_segment(b, c, d)),
                  (project_to_segment(c, a, b), c),
                  (project_to_segment(d, a, b), d)]
    return min(candidates, key=lambda pair: math.dist(*pair))


def all_relations(vertices, edges, coordinates):
    relations = []
    for first, second in combinations(sorted(vertices), 2):
        relations.append({"item_a": vertex_id(first), "item_b": vertex_id(second),
                          "geometry_type": "vertex-vertex",
                          "distance": math.dist(coordinates[first], coordinates[second]),
                          "closest_points": [coordinates[first], coordinates[second]]})
    for node in sorted(vertices):
        for edge in edges:
            if node in edge:
                continue
            projection = project_to_segment(coordinates[node],
                                            coordinates[edge[0]], coordinates[edge[1]])
            relations.append({"item_a": vertex_id(node), "item_b": edge_id(edge),
                              "geometry_type": "vertex-edge",
                              "distance": math.dist(coordinates[node], projection),
                              "closest_points": [coordinates[node], projection]})
    for first_index, first in enumerate(edges):
        for second in edges[first_index+1:]:
            if set(first) & set(second):
                continue
            closest = closest_segment_points(coordinates[first[0]], coordinates[first[1]],
                                             coordinates[second[0]], coordinates[second[1]])
            relations.append({"item_a": edge_id(first), "item_b": edge_id(second),
                              "geometry_type": "edge-edge",
                              "distance": geometry.segment_segment_distance(
                                  coordinates[first[0]], coordinates[first[1]],
                                  coordinates[second[0]], coordinates[second[1]]),
                              "closest_points": closest})
    for relation in relations:
        relation["deficit"] = max(0.0, geometry.READABILITY_CLEARANCE-relation["distance"])
    return relations


def item_geometry(item, coordinates):
    if item.startswith("vertex:"):
        point = coordinates[item[len("vertex:"):]]
        return [point], point
    first, second = item[len("edge:"):].split("|")
    points = [coordinates[first], coordinates[second]]
    return points, ((points[0][0]+points[1][0])/2, (points[0][1]+points[1][1])/2)


def analyze(vertices, edges, coordinates):
    relations = all_relations(vertices, edges, coordinates)
    violations = [item for item in relations if item["distance"] < geometry.READABILITY_CLEARANCE]
    crowding = nx.Graph()
    for relation in violations:
        crowding.add_edge(relation["item_a"], relation["item_b"])
    components = sorted((sorted(component) for component in nx.connected_components(crowding)),
                        key=lambda members: tuple(members))
    item_to_region, regions = {}, []
    for index, members in enumerate(components, 1):
        region_id = f"region:{index:03d}"
        member_set = set(members)
        region_violations = [item for item in violations
                             if item["item_a"] in member_set and item["item_b"] in member_set]
        geometry_points, representatives = [], []
        for item in members:
            points, representative = item_geometry(item, coordinates)
            geometry_points.extend(points); representatives.append(representative)
            item_to_region[item] = region_id
        xs, ys = [point[0] for point in geometry_points], [point[1] for point in geometry_points]
        centroid = (sum(point[0] for point in representatives)/len(representatives),
                    sum(point[1] for point in representatives)/len(representatives))
        diameter = max((math.dist(a, b) for a, b in combinations(geometry_points, 2)), default=0.0)
        vertex_items = sorted(item[len("vertex:"):] for item in members if item.startswith("vertex:"))
        edge_items = sorted([item[len("edge:"):].split("|") for item in members
                             if item.startswith("edge:")])
        regions.append({"region_id": region_id, "participating_vertices": vertex_items,
                        "participating_edges": edge_items,
                        "violation_count": len(region_violations),
                        "minimum_clearance": min(item["distance"] for item in region_violations),
                        "maximum_deficit": max(item["deficit"] for item in region_violations),
                        "centroid": list(centroid),
                        "bounding_box": {"min_x": min(xs), "min_y": min(ys),
                                         "max_x": max(xs), "max_y": max(ys)},
                        "width": max(xs)-min(xs), "height": max(ys)-min(ys),
                        "bounding_box_diagonal": math.hypot(max(xs)-min(xs), max(ys)-min(ys)),
                        "graph_drawing_diameter": diameter,
                        "participating_vertex_count": len(vertex_items),
                        "participating_edge_count": len(edge_items)})
    for relation in violations:
        relation["region_id"] = item_to_region[relation["item_a"]]
    regions.sort(key=lambda item: (-item["maximum_deficit"], -item["violation_count"],
                                  item["region_id"]))
    vertex_clearance, edge_clearance = {}, {}
    for node in sorted(vertices):
        item = vertex_id(node)
        nearest = min((relation for relation in relations
                       if item in {relation["item_a"], relation["item_b"]}),
                      key=lambda relation: (relation["distance"], relation["item_a"], relation["item_b"]))
        vertex_clearance[node] = {"clearance": nearest["distance"],
                                  "deficit": nearest["deficit"],
                                  "nearest_unrelated_item": (nearest["item_b"] if nearest["item_a"] == item
                                                             else nearest["item_a"])}
    for edge in edges:
        item = edge_id(edge)
        nearest = min((relation for relation in relations
                       if item in {relation["item_a"], relation["item_b"]}),
                      key=lambda relation: (relation["distance"], relation["item_a"], relation["item_b"]))
        edge_clearance["|".join(edge)] = {"edge": list(edge), "clearance": nearest["distance"],
                                          "deficit": nearest["deficit"],
                                          "nearest_unrelated_item": (nearest["item_b"] if nearest["item_a"] == item
                                                                     else nearest["item_a"])}
    type_counts = {kind: sum(item["geometry_type"] == kind for item in violations)
                   for kind in ("vertex-vertex", "vertex-edge", "edge-edge")}
    return {"relations": relations, "violations": violations, "regions": regions,
            "vertex_clearance": vertex_clearance, "edge_clearance": edge_clearance,
            "violation_type_counts": type_counts}


def svg(vertices, edges, coordinates, analysis):
    xs, ys = [p[0] for p in coordinates.values()], [p[1] for p in coordinates.values()]
    margin = 100; min_x, min_y = min(xs)-margin, min(ys)-margin
    width, height = max(xs)-min(xs)+2*margin, max(ys)-min(ys)+2*margin
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x} {min_y} {width} {height}" width="1400" height="1000">',
             f'<rect x="{min_x}" y="{min_y}" width="{width}" height="{height}" fill="white"/>',
             f'<text x="{min_x+16}" y="{min_y+28}" font-family="sans-serif" font-size="18">IAMP global 6.0-unit readability-clearance analysis</text>']
    for first, second in edges:
        parts.append(f'<line x1="{coordinates[first][0]:.17g}" y1="{coordinates[first][1]:.17g}" '
                     f'x2="{coordinates[second][0]:.17g}" y2="{coordinates[second][1]:.17g}" '
                     'stroke="#777" stroke-width="1.5"/>')
    region_by_id = {item["region_id"]: item for item in analysis["regions"]}
    color_by_id = {item["region_id"]: COLORS[index % len(COLORS)]
                   for index, item in enumerate(analysis["regions"])}
    for index, region in enumerate(analysis["regions"]):
        color = color_by_id[region["region_id"]]
        for edge in region["participating_edges"]:
            first, second = edge
            parts.append(f'<line x1="{coordinates[first][0]:.17g}" y1="{coordinates[first][1]:.17g}" '
                         f'x2="{coordinates[second][0]:.17g}" y2="{coordinates[second][1]:.17g}" '
                         f'stroke="{color}" stroke-width="4" opacity="0.45"/>')
        for node in region["participating_vertices"]:
            x, y = coordinates[node]
            parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="11" fill="none" '
                         f'stroke="{color}" stroke-width="3" opacity="0.8"/>')
        x, y = region["centroid"]
        parts.append(f'<text x="{x:.17g}" y="{y:.17g}" font-family="sans-serif" font-size="15" '
                     f'font-weight="bold" fill="{color}">{region["region_id"]}</text>')
    for node in sorted(vertices):
        x, y = coordinates[node]
        color = "#1565c0" if vertices[node]["type"] == "COMPONENT" else "#d84315"
        parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/>')
        parts.append(f'<text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[node]["name"]}</text>')
    worst = min(analysis["violations"], key=lambda item: item["distance"])
    first, second = worst["closest_points"]
    wx, wy = (first[0]+second[0])/2, (first[1]+second[1])/2
    parts.append(f'<circle cx="{wx:.17g}" cy="{wy:.17g}" r="15" fill="none" stroke="#000" stroke-width="4"/>')
    parts.append(f'<text x="{wx+18:.17g}" y="{wy-18:.17g}" font-family="sans-serif" font-size="13" font-weight="bold">WORST {worst["distance"]:.6f}</text>')
    parts.append('</svg>')
    return "\n".join(parts)+"\n"


def run(graph_path, direct_report_path, output_dir):
    vertices, edges, coordinates, graph_document, direct_report = load_inputs(
        graph_path, direct_report_path)
    original_coordinates = dict(coordinates)
    analysis = analyze(vertices, edges, coordinates)
    diagnostics = geometry.graph_geometry_diagnostics(coordinates, edges)
    if coordinates != original_coordinates or diagnostics["proper_unrelated_edge_crossing_count"] != 0:
        raise AssertionError("measurement altered or invalidated input geometry")
    violations = sorted(analysis["violations"],
                        key=lambda item: (item["distance"], item["item_a"], item["item_b"]))
    report = {"schema": "graph-relax.global-readability-clearance.v1",
              "source_graph": graph_path.as_posix(),
              "source_coordinates": f"{direct_report_path.as_posix()}#final_coordinates",
              "readability_clearance": geometry.READABILITY_CLEARANCE,
              "vertex_count": len(vertices), "electrical_incidence_count": len(edges),
              "coordinates_unchanged": coordinates == original_coordinates,
              "proper_crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
              "vertex_clearance": analysis["vertex_clearance"],
              "edge_clearance": analysis["edge_clearance"],
              "crowded_vertex_count": sum(item["deficit"] > 0
                                           for item in analysis["vertex_clearance"].values()),
              "crowded_edge_count": sum(item["deficit"] > 0
                                         for item in analysis["edge_clearance"].values()),
              "worst_vertex_clearance": min(item["clearance"] for item in analysis["vertex_clearance"].values()),
              "worst_edge_clearance": min(item["clearance"] for item in analysis["edge_clearance"].values()),
              "largest_clearance_deficit": max(item["deficit"] for item in violations),
              "violation_type_counts": analysis["violation_type_counts"],
              "crowded_region_count": len(analysis["regions"]),
              "crowded_regions": analysis["regions"],
              "worst_20_violations": violations[:20],
              "all_clearance_violations": violations}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"iamp_crowding_analysis.svg").write_text(
        svg(vertices, edges, coordinates, analysis), encoding="utf-8")
    (output_dir/"iamp_crowding_analysis_report.json").write_text(
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
    print(f"READABILITY_CLEARANCE {report['readability_clearance']}")
    print(json.dumps(report["violation_type_counts"], sort_keys=True))
    print(f"crowded regions {report['crowded_region_count']}")
    for region in report["crowded_regions"]:
        print(region["region_id"], region["minimum_clearance"], region["maximum_deficit"],
              region["violation_count"], region["participating_vertices"],
              region["participating_edges"], region["centroid"], region["bounding_box"])
    print("worst 20")
    for item in report["worst_20_violations"]:
        print(item["item_a"], item["item_b"], item["geometry_type"],
              item["distance"], item["deficit"], item["region_id"])


if __name__ == "__main__":
    main()
