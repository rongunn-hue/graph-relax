#!/usr/bin/env python3
"""Experiment 4B: degree-seeded 3-D relaxation and rigid projection search."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
from pathlib import Path

import graph_relax as gr


DEPTH = min(gr.WIDTH, gr.HEIGHT)
VIEW_DIRECTIONS = 128
ROLLS_PER_DIRECTION = 8
DIRECTIONS_3D = tuple(
    direction for direction in itertools.product((-1, 0, 1), repeat=3)
    if direction != (0, 0, 0)
)


def degree_seeded_positions_3d(circuit):
    xy = gr.degree_seeded_positions(circuit)
    positions = {ref: (point[0], point[1], 0.0) for ref, point in xy.items()}
    xs, ys = [p[0] for p in positions.values()], [p[1] for p in positions.values()]
    zmax = 0.5 * max(max(xs) - min(xs), max(ys) - min(ys))
    return positions, zmax


def terminal_positions_3d(circuit, positions):
    result = {}
    for ref in circuit.refs:
        pins = circuit.pins[ref]
        for index, pin in enumerate(pins):
            dx, dy = gr.local_terminal(index, len(pins))
            x, y, z = positions[ref]
            result[f"{ref}.{pin}"] = (x + dx, y + dy, z)
    return result


def net_segments_3d(circuit, terminals):
    segments = []
    for net in circuit.nets:
        points = [terminals[t] for t in net.terminals]
        if len(points) <= 1:
            continue
        centroid = tuple(sum(p[axis] for p in points) / len(points) for axis in range(3))
        for terminal, point in zip(net.terminals, points):
            segments.append((net.name, terminal, point, centroid))
    return segments


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def norm(a):
    return math.sqrt(dot(a, a))


def normalize(a):
    length = norm(a)
    return tuple(x / length for x in a)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def segment_distance_3d(p1, q1, p2, q2):
    """Exact closest distance between two finite 3-D segments."""
    u, v, w = subtract(q1, p1), subtract(q2, p2), subtract(p1, p2)
    a, b, c, d, e = dot(u, u), dot(u, v), dot(v, v), dot(u, w), dot(v, w)
    denominator = a * c - b * b
    small = 1e-12
    s_num, s_den = denominator, denominator
    t_num, t_den = denominator, denominator
    if denominator < small:
        s_num, s_den = 0.0, 1.0
        t_num, t_den = e, c
    else:
        s_num = b * e - c * d
        t_num = a * e - b * d
        if s_num < 0.0:
            s_num, t_num, t_den = 0.0, e, c
        elif s_num > s_den:
            s_num, t_num, t_den = s_den, e + b, c
    if t_num < 0.0:
        t_num = 0.0
        if -d < 0.0:
            s_num, s_den = 0.0, 1.0
        elif -d > a:
            s_num, s_den = 1.0, 1.0
        else:
            s_num, s_den = -d, a
    elif t_num > t_den:
        t_num = t_den
        if -d + b < 0.0:
            s_num, s_den = 0.0, 1.0
        elif -d + b > a:
            s_num, s_den = 1.0, 1.0
        else:
            s_num, s_den = -d + b, a
    sc = 0.0 if abs(s_num) < small else s_num / s_den
    tc = 0.0 if abs(t_num) < small else t_num / t_den
    delta = add(w, subtract(scale(u, sc), scale(v, tc)))
    return norm(delta)


def box_bounds_3d(point):
    x, y, z = point
    return (x - gr.WIDTH / 2, y - gr.HEIGHT / 2, z - DEPTH / 2,
            x + gr.WIDTH / 2, y + gr.HEIGHT / 2, z + DEPTH / 2)


def bounds_gap_3d(a, b):
    gaps = (max(b[i] - a[i + 3], a[i] - b[i + 3], 0.0) for i in range(3))
    return math.sqrt(sum(gap * gap for gap in gaps))


def overlap_volume_3d(a, b):
    return math.prod(max(0.0, min(a[i + 3], b[i + 3]) - max(a[i], b[i]))
                     for i in range(3))


def point_box_distance_3d(point, bounds):
    gaps = (max(bounds[i] - point[i], point[i] - bounds[i + 3], 0.0)
            for i in range(3))
    return math.sqrt(sum(gap * gap for gap in gaps))


def segment_box_distance_3d(a, b, bounds):
    """Exact piecewise-quadratic minimum distance from segment to AABB."""
    direction = subtract(b, a)
    breaks = {0.0, 1.0}
    for axis in range(3):
        if abs(direction[axis]) >= 1e-12:
            for face in (bounds[axis], bounds[axis + 3]):
                value = (face - a[axis]) / direction[axis]
                if 0.0 < value < 1.0:
                    breaks.add(value)
    ordered = sorted(breaks)
    best_squared = math.inf
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2.0
        quadratic = linear = constant = 0.0
        for axis in range(3):
            coordinate = a[axis] + direction[axis] * midpoint
            if coordinate < bounds[axis]:
                c, m = bounds[axis] - a[axis], -direction[axis]
            elif coordinate > bounds[axis + 3]:
                c, m = a[axis] - bounds[axis + 3], direction[axis]
            else:
                c, m = 0.0, 0.0
            quadratic += m * m
            linear += 2.0 * c * m
            constant += c * c
        candidates = [left, right]
        if quadratic > 0.0:
            candidates.append(max(left, min(right, -linear / (2.0 * quadratic))))
        for parameter in candidates:
            squared = quadratic * parameter * parameter + linear * parameter + constant
            best_squared = min(best_squared, squared)
    return math.sqrt(max(0.0, best_squared))


def segment_bounds_3d(a, b):
    return tuple(min(a[i], b[i]) for i in range(3)) + \
           tuple(max(a[i], b[i]) for i in range(3))


def evaluate_3d(circuit, positions):
    terminals = terminal_positions_3d(circuit, positions)
    segments = net_segments_3d(circuit, terminals)
    connection_length = sum(norm(subtract(a, b)) for _, _, a, b in segments)
    length_raw = sum(dot(subtract(a, b), subtract(a, b)) for _, _, a, b in segments)
    boxes = {ref: box_bounds_3d(positions[ref]) for ref in circuit.refs}
    overlaps = clearance_violations = 0
    overlap_raw = clearance_raw = 0.0
    for index, ref_a in enumerate(circuit.refs):
        for ref_b in circuit.refs[index + 1:]:
            volume = overlap_volume_3d(boxes[ref_a], boxes[ref_b])
            if volume > 0.0:
                overlaps += 1
                overlap_raw += 1.0 + (volume / (gr.WIDTH * gr.HEIGHT * DEPTH)) ** 2
            else:
                gap = bounds_gap_3d(boxes[ref_a], boxes[ref_b])
                if gap < gr.BODY_CLEARANCE:
                    clearance_violations += 1
                    clearance_raw += ((gr.BODY_CLEARANCE - gap) / gr.BODY_CLEARANCE) ** 2
    collisions = 0
    crowding_raw = 0.0
    segment_boxes = [segment_bounds_3d(a, b) for _, _, a, b in segments]
    for index, (net_a, terminal, a, b) in enumerate(segments):
        owner = terminal.rsplit(".", 1)[0]
        for ref in circuit.refs:
            if ref == owner or bounds_gap_3d(segment_boxes[index], boxes[ref]) >= gr.NET_BODY_CLEARANCE:
                continue
            distance = segment_box_distance_3d(a, b, boxes[ref])
            if distance < gr.NET_BODY_CLEARANCE:
                crowding_raw += ((gr.NET_BODY_CLEARANCE - distance) / gr.NET_BODY_CLEARANCE) ** 2
        for offset, (net_b, _, c, d) in enumerate(segments[index + 1:], index + 1):
            if net_a == net_b or bounds_gap_3d(segment_boxes[index], segment_boxes[offset]) >= gr.NET_NET_CLEARANCE:
                continue
            distance = segment_distance_3d(a, b, c, d)
            if distance < gr.NET_NET_CLEARANCE:
                collisions += 1
                crowding_raw += ((gr.NET_NET_CLEARANCE - distance) / gr.NET_NET_CLEARANCE) ** 2
    raw = {"length": length_raw, "overlap": overlap_raw,
           "clearance": clearance_raw, "crossing": float(collisions),
           "crowding": crowding_raw}
    energy = {name: raw[name] * gr.WEIGHTS[name] for name in gr.WEIGHTS}
    return {"connection_length_3d": connection_length,
            "component_overlaps_3d": overlaps,
            "clearance_violations_3d": clearance_violations,
            "segment_collision_proximity_events_3d": collisions,
            "energy": energy, "objective": sum(energy.values())}


def anneal_3d(circuit, start, zmax, *, iterations=gr.ANNEAL_ITERATIONS):
    rng = random.Random(gr.ANNEAL_SEED)
    state = dict(start)
    energy = evaluate_3d(circuit, state)["objective"]
    best, best_energy = dict(state), energy
    downhill = uphill = rejected = 0
    normalized_directions = tuple(normalize(direction) for direction in DIRECTIONS_3D)
    for iteration in range(iterations):
        fraction = iteration / (iterations - 1)
        temperature = gr.INITIAL_TEMPERATURE * (gr.FINAL_TEMPERATURE / gr.INITIAL_TEMPERATURE) ** fraction
        distance = gr.INITIAL_MOVE_DISTANCE * (gr.FINAL_MOVE_DISTANCE / gr.INITIAL_MOVE_DISTANCE) ** fraction
        ref = circuit.refs[rng.randrange(len(circuit.refs))]
        direction = normalized_directions[rng.randrange(len(normalized_directions))]
        magnitude = rng.uniform(0.0, distance)
        old = state[ref]
        proposed = add(old, scale(direction, magnitude))
        state[ref] = (proposed[0], proposed[1], max(-zmax, min(zmax, proposed[2])))
        candidate = evaluate_3d(circuit, state)["objective"]
        delta = candidate - energy
        if delta <= 0.0:
            accepted, downhill = True, downhill + 1
        elif rng.random() < math.exp(-delta / temperature):
            accepted, uphill = True, uphill + 1
        else:
            accepted, rejected = False, rejected + 1
        if accepted:
            energy = candidate
            if energy < best_energy:
                best, best_energy = dict(state), energy
        else:
            state[ref] = old
    return best, {"accepted_downhill_moves": downhill, "accepted_uphill_moves": uphill,
                  "rejected_moves": rejected, "best_energy": best_energy}


def quench_3d(circuit, start, zmax):
    state = dict(start)
    current = evaluate_3d(circuit, state)["objective"]
    directions = tuple(normalize(direction) for direction in DIRECTIONS_3D)
    step, sweeps = gr.INITIAL_STEP, 0
    while sweeps < gr.MAX_SWEEPS:
        sweeps += 1
        accepted = 0
        for ref in circuit.refs:
            point = state[ref]
            best_point, best_energy = point, current
            for direction in directions:
                proposed = add(point, scale(direction, step))
                proposed = (proposed[0], proposed[1], max(-zmax, min(zmax, proposed[2])))
                state[ref] = proposed
                energy = evaluate_3d(circuit, state)["objective"]
                if energy < best_energy - 1e-9:
                    best_point, best_energy = proposed, energy
            state[ref] = best_point
            if best_energy < current - 1e-9:
                current, accepted = best_energy, accepted + 1
        if accepted == 0:
            if step <= gr.MIN_STEP:
                return state, sweeps, "no improving move at minimum step"
            step = max(gr.MIN_STEP, step / 2.0)
    return state, sweeps, "maximum sweeps reached"


def orientations(direction_count=VIEW_DIRECTIONS, rolls=ROLLS_PER_DIRECTION):
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    result = []
    for direction_index in range(direction_count):
        z = 1.0 - 2.0 * (direction_index + 0.5) / direction_count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        phi = golden_angle * direction_index
        view = (radius * math.cos(phi), radius * math.sin(phi), z)
        reference = (0.0, 0.0, 1.0) if abs(z) < 0.9 else (0.0, 1.0, 0.0)
        base_u = normalize(cross(reference, view))
        base_v = cross(view, base_u)
        for roll_index in range(rolls):
            roll = 2.0 * math.pi * roll_index / rolls
            u = add(scale(base_u, math.cos(roll)), scale(base_v, math.sin(roll)))
            v = add(scale(base_u, -math.sin(roll)), scale(base_v, math.cos(roll)))
            result.append({"direction_index": direction_index, "roll_index": roll_index,
                           "view": view, "roll_radians": roll, "u": u, "v": v})
    return result


def project_centers(positions, orientation):
    return {ref: (dot(point, orientation["u"]), dot(point, orientation["v"]))
            for ref, point in positions.items()}


def projection_search(circuit, positions, orientation_set=None):
    candidates = []
    best_positions = best_metrics = best_orientation = None
    best_index = None
    for index, orientation in enumerate(orientation_set or orientations()):
        projected = project_centers(positions, orientation)
        metrics = gr.evaluate(circuit, projected)
        candidates.append({"index": index, "direction_index": orientation["direction_index"],
                           "roll_index": orientation["roll_index"],
                           "objective": gr.rounded(metrics["objective"]),
                           "crossings": metrics["net_crossings"]})
        if best_metrics is None or metrics["objective"] < best_metrics["objective"] - 1e-9:
            best_positions, best_metrics = projected, metrics
            best_orientation, best_index = orientation, index
    return best_positions, best_metrics, best_orientation, best_index, candidates


def coordinate_document_3d(circuit, positions, zmax):
    terminals = terminal_positions_3d(circuit, positions)
    return {"circuit": circuit.name, "zmax": gr.rounded(zmax), "components": {
        ref: {"x": gr.rounded(positions[ref][0]), "y": gr.rounded(positions[ref][1]),
              "z": gr.rounded(positions[ref][2]), "width": gr.WIDTH,
              "height": gr.HEIGHT, "depth": DEPTH,
              "terminals": {pin: {"x": gr.rounded(terminals[f"{ref}.{pin}"][0]),
                                    "y": gr.rounded(terminals[f"{ref}.{pin}"][1]),
                                    "z": gr.rounded(terminals[f"{ref}.{pin}"][2])}
                            for pin in circuit.pins[ref]}}
        for ref in circuit.refs}}


def rounded_3d_metrics(metrics):
    result = {key: (gr.rounded(value) if isinstance(value, float) else value)
              for key, value in metrics.items() if key != "energy"}
    result["energy"] = {key: gr.rounded(value) for key, value in metrics["energy"].items()}
    return result


def extents(positions):
    return {axis: {"min": gr.rounded(min(p[i] for p in positions.values())),
                   "max": gr.rounded(max(p[i] for p in positions.values()))}
            for i, axis in enumerate(("x", "y", "z"))}


def canonical_projection(positions, axes):
    return {ref: (point[axes[0]], point[axes[1]]) for ref, point in positions.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    circuit = gr.parse_circuit(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start, zmax = degree_seeded_positions_3d(circuit)
    best, anneal_stats = anneal_3d(circuit, start, zmax)
    final, quench_sweeps, reason = quench_3d(circuit, best, zmax)
    final_3d_metrics = evaluate_3d(circuit, final)
    coordinates_3d_path = args.output_dir / "degree_seeded_3d_coordinates.json"
    gr.write_json(coordinates_3d_path, coordinate_document_3d(circuit, final, zmax))
    metrics_3d = {
        "zmax": gr.rounded(zmax), "component_depth": DEPTH,
        "collision_predicate": "minimum distance between unrelated 3-D net segments < 8",
        "movement_directions": [list(direction) for direction in DIRECTIONS_3D],
        "coordinate_extents": extents(final), "metrics": rounded_3d_metrics(final_3d_metrics),
        "annealing": {"seed": gr.ANNEAL_SEED, "initial_temperature": gr.INITIAL_TEMPERATURE,
                      "final_temperature": gr.FINAL_TEMPERATURE, "cooling_schedule": "geometric",
                      "initial_move_distance": gr.INITIAL_MOVE_DISTANCE,
                      "final_move_distance": gr.FINAL_MOVE_DISTANCE,
                      "move_distance_schedule": "geometric", "iterations": gr.ANNEAL_ITERATIONS,
                      "accepted_downhill_moves": anneal_stats["accepted_downhill_moves"],
                      "accepted_uphill_moves": anneal_stats["accepted_uphill_moves"],
                      "rejected_moves": anneal_stats["rejected_moves"],
                      "best_energy_encountered": gr.rounded(anneal_stats["best_energy"])},
        "quench": {"final_energy": gr.rounded(final_3d_metrics["objective"]),
                   "sweeps": quench_sweeps, "stopping_reason": reason},
        "coordinate_sha256": hashlib.sha256(coordinates_3d_path.read_bytes()).hexdigest(),
    }
    gr.write_json(args.output_dir / "degree_seeded_3d_metrics.json", metrics_3d)
    best_2d, best_metrics, orientation, best_index, candidates = projection_search(circuit, final)
    best_coordinates_path = args.output_dir / "best_3d_projection_coordinates.json"
    gr.write_json(best_coordinates_path, gr.coordinate_document(circuit, best_2d))
    gr.write_svg(args.output_dir / "best_3d_projection.svg", circuit, best_2d)
    baseline = json.loads((args.output_dir / "degree_seeded_metrics.json").read_text(encoding="utf-8"))
    baseline_metrics = baseline["states"]["experiment_4a_degree_seeded_annealed_and_quenched"]
    changes = {
        "connection_length_percent": gr.rounded(100.0 * (
            best_metrics["total_connection_length"] - baseline_metrics["total_connection_length"])
            / baseline_metrics["total_connection_length"]),
        "overlaps_percent": None if baseline_metrics["component_overlaps"] == 0 else gr.rounded(
            100.0 * (best_metrics["component_overlaps"] - baseline_metrics["component_overlaps"])
            / baseline_metrics["component_overlaps"]),
        "clearance_violations_percent": gr.rounded(100.0 * (
            best_metrics["clearance_violations"] - baseline_metrics["clearance_violations"])
            / baseline_metrics["clearance_violations"]),
        "crossings_percent": gr.rounded(100.0 * (
            best_metrics["net_crossings"] - baseline_metrics["net_crossings"])
            / baseline_metrics["net_crossings"]),
        "crowding_percent": gr.rounded(100.0 * (
            best_metrics["energy"]["crowding"] - baseline_metrics["energy"]["crowding"])
            / baseline_metrics["energy"]["crowding"]),
        "objective_percent": gr.rounded(100.0 * (
            best_metrics["objective"] - baseline_metrics["objective"])
            / baseline_metrics["objective"]),
    }
    projection_metrics = {
        "orientations_tested": len(candidates), "best_projection_index": best_index,
        "best_orientation": {"direction_index": orientation["direction_index"],
                             "roll_index": orientation["roll_index"],
                             "view_direction": [gr.rounded(x) for x in orientation["view"]],
                             "roll_radians": gr.rounded(orientation["roll_radians"])},
        "metrics": gr.rounded_metrics(best_metrics), "experiment_4a": baseline_metrics,
        "percentage_changes_from_experiment_4a": changes,
        "projections_beating_experiment_4a_objective": sum(
            candidate["objective"] < baseline_metrics["objective"] for candidate in candidates),
        "projections_beating_experiment_4a_crossings": sum(
            candidate["crossings"] < baseline_metrics["net_crossings"] for candidate in candidates),
        "coordinate_sha256": hashlib.sha256(best_coordinates_path.read_bytes()).hexdigest(),
    }
    gr.write_json(args.output_dir / "best_3d_projection_metrics.json", projection_metrics)
    gr.write_json(args.output_dir / "3d_projection_search.json", {
        "direction_count": VIEW_DIRECTIONS, "rolls_per_direction": ROLLS_PER_DIRECTION,
        "total_projections": len(candidates),
        "sampling": "Fibonacci sphere: z=1-2(i+0.5)/N, phi=i*pi*(3-sqrt(5)); rolls=2*pi*r/R",
        "enumeration": "direction index outermost, roll index innermost",
        "candidates": candidates})
    for filename, axes in (("3d_view_xy.svg", (0, 1)), ("3d_view_xz.svg", (0, 2)),
                           ("3d_view_yz.svg", (1, 2))):
        gr.write_svg(args.output_dir / filename, circuit, canonical_projection(final, axes))
    print(json.dumps({"three_d": metrics_3d, "projection": projection_metrics},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
