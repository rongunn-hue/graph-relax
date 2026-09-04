#!/usr/bin/env python3
"""Analytic event-angle endpoint rotation for the spaced nearness drawing."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import placement_geometry as geometry
import graph_first_experiment_endpoint_rotation as sweep


TAU = 2.0 * math.pi
EVENT_DEDUP_TOLERANCE = 1e-11
NUMERICAL_SAFETY_FACTOR = 16.0


def normalize(angle):
    return angle % TAU


def signed_delta(angle, current):
    delta = (angle - current + math.pi) % TAU - math.pi
    if math.isclose(delta, -math.pi, abs_tol=EVENT_DEDUP_TOLERANCE):
        return math.pi
    return delta


def deduplicate_angles(angles):
    ordered = sorted(normalize(value) for value in angles)
    result = []
    for value in ordered:
        if not result or abs(value - result[-1]) > EVENT_DEDUP_TOLERANCE:
            result.append(value)
    if len(result) > 1 and abs((result[0] + TAU) - result[-1]) <= EVENT_DEDUP_TOLERANCE:
        result.pop()
    return result


def circle_circle_event_angles(pivot, radius, center, other_radius):
    dx, dy = center[0] - pivot[0], center[1] - pivot[1]
    distance = math.hypot(dx, dy)
    tolerance = geometry.GEOMETRY_TOLERANCE
    if distance <= tolerance:
        if math.isclose(radius, other_radius, abs_tol=tolerance):
            return [], "COINCIDENT_CIRCLES"
        return [], "CONCENTRIC_NO_INTERSECTION"
    if distance > radius + other_radius + tolerance or \
            distance < abs(radius - other_radius) - tolerance:
        return [], "NO_INTERSECTION"
    along = (radius * radius - other_radius * other_radius + distance * distance) / (2 * distance)
    height2 = radius * radius - along * along
    if height2 < -tolerance:
        return [], "NO_INTERSECTION"
    base_x = pivot[0] + along * dx / distance
    base_y = pivot[1] + along * dy / distance
    if abs(height2) <= tolerance:
        return [math.atan2(base_y - pivot[1], base_x - pivot[0])], "TANGENT"
    height = math.sqrt(max(0.0, height2))
    offset_x, offset_y = -dy * height / distance, dx * height / distance
    return [math.atan2(base_y + offset_y - pivot[1], base_x + offset_x - pivot[0]),
            math.atan2(base_y - offset_y - pivot[1], base_x - offset_x - pivot[0])], \
           "TWO_INTERSECTIONS"


def circle_segment_event_angles(pivot, radius, first, second):
    dx, dy = second[0] - first[0], second[1] - first[1]
    fx, fy = first[0] - pivot[0], first[1] - pivot[1]
    a = dx * dx + dy * dy
    if a <= geometry.GEOMETRY_TOLERANCE ** 2:
        return []
    b = 2.0 * (fx * dx + fy * dy)
    c = fx * fx + fy * fy - radius * radius
    discriminant = b * b - 4.0 * a * c
    if discriminant < -geometry.GEOMETRY_TOLERANCE:
        return []
    roots = [-b / (2 * a)] if abs(discriminant) <= geometry.GEOMETRY_TOLERANCE else [
        (-b - math.sqrt(max(0.0, discriminant))) / (2 * a),
        (-b + math.sqrt(max(0.0, discriminant))) / (2 * a)]
    angles = []
    for fraction in roots:
        if -geometry.GEOMETRY_TOLERANCE <= fraction <= 1 + geometry.GEOMETRY_TOLERANCE:
            x, y = first[0] + fraction * dx, first[1] + fraction * dy
            angles.append(math.atan2(y - pivot[1], x - pivot[0]))
    return angles


def event_candidates(graph, positions, moving, pivot, other_edge):
    pivot_point = positions[pivot]
    radius = math.dist(positions[moving], pivot_point)
    current = math.atan2(positions[moving][1] - pivot_point[1],
                         positions[moving][0] - pivot_point[0])
    angles = [current]
    event_counts = {"crossing_ray": 0, "spacing": 0, "vertex_edge": 0}
    degeneracies = []
    for node in other_edge:
        angles.append(math.atan2(positions[node][1] - pivot_point[1],
                                 positions[node][0] - pivot_point[0]))
        event_counts["crossing_ray"] += 1
    for node in sorted(graph):
        if node == moving:
            continue
        found, status = circle_circle_event_angles(
            pivot_point, radius, positions[node], geometry.MIN_NODE_DISTANCE)
        angles.extend(found)
        event_counts["spacing"] += len(found)
        if status == "COINCIDENT_CIRCLES":
            degeneracies.append({"type": status, "vertex": node})
    for first, second in sorted((min(a, b), max(a, b)) for a, b in graph.edges):
        if moving in {first, second}:
            continue
        found = circle_segment_event_angles(
            pivot_point, radius, positions[first], positions[second])
        angles.extend(found)
        event_counts["vertex_edge"] += len(found)
    events = deduplicate_angles(angles)
    drawing_scale = max(math.dist(pivot_point, point) for point in positions.values())
    epsilon = NUMERICAL_SAFETY_FACTOR * geometry.GEOMETRY_TOLERANCE * \
        max(1.0, drawing_scale) / max(radius, geometry.GEOMETRY_TOLERANCE)
    candidates = []
    for event in events:
        for candidate in (event - epsilon, event, event + epsilon):
            delta = signed_delta(candidate, current)
            if abs(delta) > EVENT_DEDUP_TOLERANCE:
                candidates.append(delta)
    candidates = sorted(set(round(value, 15) for value in candidates),
                        key=lambda value: (abs(value), value))
    return candidates, event_counts, degeneracies, current, radius, epsilon


def analytic_best_move(graph, positions, diagnostics, stats):
    before = diagnostics["proper_unrelated_edge_crossing_count"]
    options = []
    for crossing in diagnostics["proper_unrelated_edge_crossings"]:
        target = sweep.crossing_key(crossing)
        edge_a, edge_b = tuple(crossing["edge_a"]), tuple(crossing["edge_b"])
        motions = [(edge_a[0], edge_a[1], edge_b), (edge_a[1], edge_a[0], edge_b),
                   (edge_b[0], edge_b[1], edge_a), (edge_b[1], edge_b[0], edge_a)]
        for moving, pivot, other_edge in motions:
            candidates, counts, degeneracies, old_angle, radius, epsilon = event_candidates(
                graph, positions, moving, pivot, other_edge)
            stats["generated_event_angle_candidates"] += len(candidates)
            stats["degenerate_events"].extend(
                {"moving": moving, "pivot": pivot, **item} for item in degeneracies)
            for delta_radians in candidates:
                stats["candidate_evaluations"] += 1
                signed_degrees = math.degrees(delta_radians)
                candidate, candidate_diagnostics, useful = sweep.evaluate(
                    graph, positions, moving, pivot, signed_degrees, target)
                count = candidate_diagnostics["proper_unrelated_edge_crossing_count"]
                if not useful or count >= before:
                    continue
                displacement = math.dist(positions[moving], candidate[moving])
                options.append({
                    "target": target, "moving": moving, "pivot": pivot,
                    "other_edge": other_edge, "old_angle": old_angle,
                    "new_angle": normalize(old_angle + delta_radians), "radius": radius,
                    "signed_radians": delta_radians, "signed_degrees": signed_degrees,
                    "displacement": displacement, "positions": candidate,
                    "diagnostics": candidate_diagnostics, "event_counts": counts,
                    "angular_epsilon_radians": epsilon,
                    "selection_key": (count, -(before-count), abs(delta_radians), displacement,
                                      moving, pivot, delta_radians),
                })
    return min(options, key=lambda item: item["selection_key"]) if options else None


def analytic_search(graph, initial_positions):
    positions = dict(initial_positions)
    diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
    stats = {"candidate_evaluations": 0, "generated_event_angle_candidates": 0,
             "degenerate_events": []}
    moves, states = [], [dict(positions)]
    history = [diagnostics["proper_unrelated_edge_crossings"]]
    while diagnostics["proper_unrelated_edge_crossing_count"]:
        option = analytic_best_move(graph, positions, diagnostics, stats)
        if option is None:
            break
        complete = geometry.graph_geometry_diagnostics(option["positions"], graph.edges)
        if complete != option["diagnostics"]:
            raise AssertionError("analytic candidate full validation mismatch")
        before_keys = {sweep.crossing_key(item) for item in
                       diagnostics["proper_unrelated_edge_crossings"]}
        after_keys = {sweep.crossing_key(item) for item in
                      complete["proper_unrelated_edge_crossings"]}
        new_length = math.dist(option["positions"][option["moving"]],
                               option["positions"][option["pivot"]])
        if not math.isclose(new_length, option["radius"], rel_tol=1e-12, abs_tol=1e-9):
            raise AssertionError("analytic rotation changed pivot-edge radius")
        moves.append({
            "iteration": len(moves)+1,
            "targeted_crossing": {"edge_a": list(option["target"][0]),
                                  "edge_b": list(option["target"][1])},
            "moving_vertex": option["moving"], "pivot_vertex": option["pivot"],
            "preserved_radius": option["radius"],
            "old_angle_radians": option["old_angle"], "new_angle_radians": option["new_angle"],
            "signed_rotation_radians": option["signed_radians"],
            "signed_rotation_degrees": option["signed_degrees"],
            "absolute_rotation_degrees": abs(option["signed_degrees"]),
            "moving_vertex_displacement": option["displacement"],
            "crossing_count_before": diagnostics["proper_unrelated_edge_crossing_count"],
            "crossing_count_after": complete["proper_unrelated_edge_crossing_count"],
            "crossings_removed": [{"edge_a": list(key[0]), "edge_b": list(key[1])}
                                  for key in sorted(before_keys-after_keys)],
            "crossings_newly_created": [{"edge_a": list(key[0]), "edge_b": list(key[1])}
                                        for key in sorted(after_keys-before_keys)],
            "minimum_node_distance_after": complete["minimum_node_distance"],
            "event_counts": option["event_counts"],
            "angular_epsilon_radians": option["angular_epsilon_radians"],
            "start_position": list(positions[option["moving"]]),
            "end_position": list(option["positions"][option["moving"]]),
            "pivot_position": list(positions[option["pivot"]]),
            "preserved_edge_length_before": option["radius"],
        })
        positions, diagnostics = option["positions"], complete
        states.append(dict(positions))
        history.append(diagnostics["proper_unrelated_edge_crossings"])
    return positions, diagnostics, moves, states, history, stats


def verify_sweep_stages(graph, start, sweep_report):
    positions = dict(start)
    stages = []
    for index in range(len(sweep_report["moves"])+1):
        diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
        expected = sweep_report["crossing_history"][index]
        if diagnostics["proper_unrelated_edge_crossings"] != expected:
            raise AssertionError(f"stored sweep crossing history mismatch at stage {index}")
        stages.append({"stage": "start" if index == 0 else f"after_move_{index}",
                       "crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
                       "crossings": diagnostics["proper_unrelated_edge_crossings"]})
        if index < len(sweep_report["moves"]):
            move = sweep_report["moves"][index]
            positions[move["moving_vertex"]] = tuple(move["end_position"])
    return stages


def run(graph_path, coordinate_path, mds_report_path, sweep_report_path, output_dir):
    graph, start, initial, spacing = sweep.load_start(graph_path, coordinate_path, mds_report_path)
    stored_sweep = json.loads(sweep_report_path.read_text(encoding="utf-8"))
    sweep_stages = verify_sweep_stages(graph, start, stored_sweep)

    sweep_stats = {"candidate_evaluations": 0}
    then = time.perf_counter()
    _, sweep_final, sweep_moves, _, _ = sweep.run_search(graph, start, sweep_stats)
    sweep_runtime = time.perf_counter() - then
    then = time.perf_counter()
    final, diagnostics, moves, states, history, analytic_stats = analytic_search(graph, start)
    analytic_runtime = time.perf_counter() - then
    report = {
        "schema": "graph-relax.analytic-event-endpoint-rotation.v1",
        "sweep_stage_verification": sweep_stages,
        "analytic_initial_crossing_count": initial["proper_unrelated_edge_crossing_count"],
        "analytic_initial_crossings": initial["proper_unrelated_edge_crossings"],
        "analytic_final_crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
        "analytic_final_crossings": diagnostics["proper_unrelated_edge_crossings"],
        "analytic_moves": moves,
        "analytic_crossing_history": history,
        "analytic_stopping_reason": ("ZERO_CROSSINGS" if not diagnostics["proper_unrelated_edge_crossing_count"]
                                     else "NO_ADMISSIBLE_IMPROVING_ANALYTIC_EVENT_MOVE"),
        "analytic_final_geometry": diagnostics,
        "comparison": {
            "sweep": {"initial_crossings": initial["proper_unrelated_edge_crossing_count"],
                      "final_crossings": sweep_final["proper_unrelated_edge_crossing_count"],
                      "accepted_moves": len(sweep_moves),
                      "candidate_evaluations": sweep_stats["candidate_evaluations"],
                      "runtime_seconds": sweep_runtime,
                      "total_absolute_angular_movement_degrees": sum(abs(item["signed_rotation_degrees"])
                                                                     for item in sweep_moves),
                      "largest_absolute_angular_movement_degrees": max(abs(item["signed_rotation_degrees"])
                                                                       for item in sweep_moves),
                      "total_euclidean_node_displacement": sum(item["moving_vertex_displacement"]
                                                               for item in sweep_moves)},
            "analytic": {"initial_crossings": initial["proper_unrelated_edge_crossing_count"],
                         "final_crossings": diagnostics["proper_unrelated_edge_crossing_count"],
                         "accepted_moves": len(moves),
                         "candidate_evaluations": analytic_stats["candidate_evaluations"],
                         "generated_event_angle_candidates": analytic_stats["generated_event_angle_candidates"],
                         "runtime_seconds": analytic_runtime,
                         "total_absolute_angular_movement_degrees": sum(abs(item["signed_rotation_degrees"])
                                                                        for item in moves),
                         "largest_absolute_angular_movement_degrees": max(abs(item["signed_rotation_degrees"])
                                                                          for item in moves),
                         "total_euclidean_node_displacement": sum(item["moving_vertex_displacement"]
                                                                  for item in moves)},
        },
        "analytic_event_diagnostics": analytic_stats,
        "final_coordinates": {node: list(final[node]) for node in sorted(final)},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    viewbox = sweep.viewbox_for(states)
    (output_dir / "iamp_nearness_rotation_analytic_start.svg").write_text(
        sweep.static_svg(graph, start, "Analytic endpoint rotation start", viewbox), encoding="utf-8")
    (output_dir / "iamp_nearness_rotation_analytic_end.svg").write_text(
        sweep.static_svg(graph, final, "Analytic endpoint rotation end", viewbox), encoding="utf-8")
    (output_dir / "iamp_nearness_rotation_analytic_motion.svg").write_text(
        sweep.motion_svg(graph, states, moves, viewbox), encoding="utf-8")
    sweep.write_json(output_dir / "iamp_nearness_rotation_analytic_report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base/"iamp_graph.json")
    parser.add_argument("--coordinates", type=Path, default=base/"iamp_nearness_mds.json")
    parser.add_argument("--mds-report", type=Path, default=base/"iamp_nearness_mds_report.json")
    parser.add_argument("--sweep-report", type=Path, default=base/"iamp_nearness_rotation_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    print(json.dumps(run(args.graph, args.coordinates, args.mds_report,
                         args.sweep_report, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
