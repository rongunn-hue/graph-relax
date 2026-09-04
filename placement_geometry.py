"""Shared geometric invariants for all graph placement methods.

Coordinates are abstract drawing units. Uniform spacing enforcement changes
only global scale about an explicit anchor; it cannot alter angles, distance
ratios, aspect ratio, cyclic order, or crossing topology.
"""

from __future__ import annotations

import math
from itertools import combinations
from typing import Iterable, Mapping


MIN_NODE_DISTANCE = 1.0
GEOMETRY_TOLERANCE = 1e-9
EDGE_STROKE_WIDTH = 1.5
MIN_VISIBLE_EDGE_GAP = 4.5
REQUIRED_CENTERLINE_CLEARANCE = EDGE_STROKE_WIDTH + MIN_VISIBLE_EDGE_GAP
READABILITY_CLEARANCE = 6.0


def rebuild_straight_edge_segments(
        coordinates: Mapping[str, tuple[float, float]],
        edges: Iterable[tuple[str, str]]) -> list[dict]:
    """Materialize every edge from the authoritative endpoint coordinates."""
    points = _coordinates(coordinates)
    normalized_edges = sorted((min(first, second), max(first, second))
                              for first, second in edges)
    return [{"edge": [first, second],
             "segment_start": points[first], "segment_end": points[second]}
            for first, second in normalized_edges]


def validate_rebuilt_edge_segments(
        coordinates: Mapping[str, tuple[float, float]],
        edges: Iterable[tuple[str, str]],
        segments: Iterable[dict],
        tolerance: float = GEOMETRY_TOLERANCE) -> dict:
    """Assert that materialized segments exactly represent graph connectivity."""
    points = _coordinates(coordinates)
    expected = sorted((min(first, second), max(first, second))
                      for first, second in edges)
    materialized = list(segments)
    actual = [tuple(item["edge"]) for item in materialized]
    failures = []
    if actual != expected:
        failures.append("EDGE_ENDPOINT_IDENTITIES_CHANGED")
    if len(materialized) != len(expected):
        failures.append("EDGE_COUNT_CHANGED")
    for item in materialized:
        first, second = item["edge"]
        if first not in points or second not in points:
            failures.append(f"UNKNOWN_ENDPOINT:{first}:{second}")
            continue
        if math.dist(tuple(item["segment_start"]), points[first]) > tolerance:
            failures.append(f"STALE_SEGMENT_START:{first}:{second}")
        if math.dist(tuple(item["segment_end"]), points[second]) > tolerance:
            failures.append(f"STALE_SEGMENT_END:{first}:{second}")
    return {
        "valid": not failures,
        "graph_edge_count": len(expected),
        "materialized_edge_count": len(materialized),
        "electrical_incidence_count": len(expected),
        "endpoint_identities_unchanged": actual == expected,
        "failures": failures,
    }


def _coordinates(coordinates: Mapping[str, tuple[float, float]]):
    return {node: (float(point[0]), float(point[1]))
            for node, point in sorted(coordinates.items())}


def node_spacing_diagnostics(coordinates: Mapping[str, tuple[float, float]],
                             tolerance: float = GEOMETRY_TOLERANCE) -> dict:
    points = _coordinates(coordinates)
    pairs = [(first, second, math.dist(points[first], points[second]))
             for first, second in combinations(points, 2)]
    minimum = min((distance for _, _, distance in pairs), default=None)
    coincident = [[first, second] for first, second, distance in pairs
                  if distance <= tolerance]
    return {
        "minimum_node_distance": minimum,
        "coincident_vertex_pair_count": len(coincident),
        "coincident_vertex_pairs": coincident,
    }


def enforce_minimum_node_spacing(
        coordinates: Mapping[str, tuple[float, float]],
        minimum_distance: float = MIN_NODE_DISTANCE,
        anchor: tuple[float, float] = (0.0, 0.0),
        tolerance: float = GEOMETRY_TOLERANCE) -> tuple[dict[str, tuple[float, float]], dict]:
    """Uniformly scale points about *anchor* to enforce minimum separation.

    Coincident vertices make scale-only enforcement mathematically impossible.
    In that case coordinates are returned unchanged and the report explicitly
    marks the placement as failed; no division or local displacement occurs.
    """
    if minimum_distance <= 0:
        raise ValueError("minimum_distance must be positive")
    points = _coordinates(coordinates)
    raw = node_spacing_diagnostics(points, tolerance)
    base_report = {
        "requested_minimum_node_distance": float(minimum_distance),
        "raw_minimum_node_distance": raw["minimum_node_distance"],
        "anchor": [float(anchor[0]), float(anchor[1])],
        "coincident_vertex_pair_count": raw["coincident_vertex_pair_count"],
        "coincident_vertex_pairs": raw["coincident_vertex_pairs"],
    }
    if len(points) < 2:
        report = {**base_report, "status": "SATISFIED", "applied_uniform_scale_factor": 1.0,
                  "final_minimum_node_distance": raw["minimum_node_distance"]}
        return points, report
    if raw["coincident_vertex_pair_count"]:
        report = {**base_report, "status": "FAILED_COINCIDENT_VERTICES",
                  "applied_uniform_scale_factor": None,
                  "final_minimum_node_distance": raw["minimum_node_distance"]}
        return points, report
    current_minimum = raw["minimum_node_distance"]
    if current_minimum is None or current_minimum + tolerance >= minimum_distance:
        scale = 1.0
    else:
        scale = minimum_distance / current_minimum
    ax, ay = float(anchor[0]), float(anchor[1])
    scaled = {node: (ax + scale * (x - ax), ay + scale * (y - ay))
              for node, (x, y) in points.items()}
    final = node_spacing_diagnostics(scaled, tolerance)
    report = {
        **base_report,
        "status": "SATISFIED",
        "applied_uniform_scale_factor": scale,
        "final_minimum_node_distance": final["minimum_node_distance"],
    }
    return scaled, report


def resolve_coincident_vertices(
        coordinates: Mapping[str, tuple[float, float]],
        minimum_distance: float = MIN_NODE_DISTANCE,
        tolerance: float = GEOMETRY_TOLERANCE) -> tuple[dict[str, tuple[float, float]], dict]:
    """Resolve coincident clusters without graph- or domain-specific knowledge.

    Each connected coincidence cluster is replaced by a lexically ordered
    regular polygon centered on the original shared point. Its adjacent chord
    is exactly *minimum_distance*. This is deterministic symmetry breaking;
    the subsequent uniform-spacing stage remains responsible for separation
    from vertices outside the cluster.
    """
    if minimum_distance <= 0:
        raise ValueError("minimum_distance must be positive")
    points = _coordinates(coordinates)
    coincidence_graph: dict[str, set[str]] = {node: set() for node in points}
    for first, second in combinations(points, 2):
        if math.dist(points[first], points[second]) <= tolerance:
            coincidence_graph[first].add(second)
            coincidence_graph[second].add(first)
    visited = set()
    clusters = []
    for node in points:
        if node in visited or not coincidence_graph[node]:
            continue
        stack, members = [node], []
        visited.add(node)
        while stack:
            current = stack.pop()
            members.append(current)
            for neighbor in sorted(coincidence_graph[current], reverse=True):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        clusters.append(sorted(members))
    clusters.sort()
    resolved = dict(points)
    records = []
    for members in clusters:
        center_x = sum(points[node][0] for node in members) / len(members)
        center_y = sum(points[node][1] for node in members) / len(members)
        radius = minimum_distance / (2.0 * math.sin(math.pi / len(members)))
        for index, node in enumerate(members):
            angle = 2.0 * math.pi * index / len(members)
            resolved[node] = (center_x + radius * math.cos(angle),
                              center_y + radius * math.sin(angle))
        records.append({"vertices": members, "original_point": [center_x, center_y],
                        "polygon_radius": radius})
    remaining = node_spacing_diagnostics(resolved, tolerance)
    return resolved, {
        "status": "RESOLVED" if clusters else "NOT_NEEDED",
        "method": "lexically ordered regular polygon at each coincident cluster",
        "cluster_count": len(clusters),
        "clusters": records,
        "remaining_coincident_vertex_pair_count": remaining["coincident_vertex_pair_count"],
        "remaining_coincident_vertex_pairs": remaining["coincident_vertex_pairs"],
    }


def cross(first, second, third) -> float:
    return ((second[0] - first[0]) * (third[1] - first[1]) -
            (second[1] - first[1]) * (third[0] - first[0]))


def proper_segment_crossing(a, b, c, d,
                            tolerance: float = GEOMETRY_TOLERANCE) -> bool:
    values = (cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b))
    epsilon = tolerance * max(1.0, math.dist(a, b), math.dist(c, d))
    return (((values[0] > epsilon and values[1] < -epsilon) or
             (values[0] < -epsilon and values[1] > epsilon)) and
            ((values[2] > epsilon and values[3] < -epsilon) or
             (values[2] < -epsilon and values[3] > epsilon)))


def point_on_segment_interior(point, first, second,
                              tolerance: float = GEOMETRY_TOLERANCE) -> bool:
    length = math.dist(first, second)
    if length <= tolerance or abs(cross(first, second, point)) > tolerance * max(1.0, length):
        return False
    fraction = ((point[0] - first[0]) * (second[0] - first[0]) +
                (point[1] - first[1]) * (second[1] - first[1])) / (length * length)
    return tolerance < fraction < 1.0 - tolerance


def point_segment_distance(point, first, second) -> float:
    dx, dy = second[0] - first[0], second[1] - first[1]
    denominator = dx*dx + dy*dy
    if denominator == 0.0:
        return math.dist(point, first)
    fraction = ((point[0]-first[0])*dx + (point[1]-first[1])*dy) / denominator
    fraction = max(0.0, min(1.0, fraction))
    projection = (first[0] + fraction*dx, first[1] + fraction*dy)
    return math.dist(point, projection)


def segment_segment_distance(a, b, c, d,
                             tolerance: float = GEOMETRY_TOLERANCE) -> float:
    if proper_segment_crossing(a, b, c, d, tolerance):
        return 0.0
    # Inclusive endpoint/contact intersections also have zero separation.
    if (point_segment_distance(a, c, d) <= tolerance or
            point_segment_distance(b, c, d) <= tolerance or
            point_segment_distance(c, a, b) <= tolerance or
            point_segment_distance(d, a, b) <= tolerance):
        return 0.0
    return min(point_segment_distance(a, c, d), point_segment_distance(b, c, d),
               point_segment_distance(c, a, b), point_segment_distance(d, a, b))


def graph_geometry_diagnostics(
        coordinates: Mapping[str, tuple[float, float]],
        edges: Iterable[tuple[str, str]],
        tolerance: float = GEOMETRY_TOLERANCE) -> dict:
    points = _coordinates(coordinates)
    normalized_edges = sorted((min(first, second), max(first, second))
                              for first, second in edges)
    spacing = node_spacing_diagnostics(points, tolerance)
    crossings = []
    for index, (first, second) in enumerate(normalized_edges):
        for third, fourth in normalized_edges[index + 1:]:
            if {first, second} & {third, fourth}:
                continue
            if proper_segment_crossing(points[first], points[second],
                                       points[third], points[fourth], tolerance):
                crossings.append({"edge_a": [first, second], "edge_b": [third, fourth]})
    on_edges = []
    for node in points:
        for first, second in normalized_edges:
            if node not in {first, second} and point_on_segment_interior(
                    points[node], points[first], points[second], tolerance):
                on_edges.append({"vertex": node, "edge": [first, second]})
    return {
        **spacing,
        "proper_unrelated_edge_crossing_count": len(crossings),
        "proper_unrelated_edge_crossings": crossings,
        "vertex_on_unrelated_edge_interior_count": len(on_edges),
        "vertices_on_unrelated_edge_interiors": on_edges,
    }
