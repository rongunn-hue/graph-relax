#!/usr/bin/env python3
"""Experiment 4C: closeness-organized 2-D initialization."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
from pathlib import Path

import graph_relax as gr


INITIALIZATION_SWEEPS = 100
OVERLAP_RESOLUTION_PASSES = 100
MOVEMENT_TOLERANCE = 1e-6
COMPONENT_DIAGONAL = math.hypot(gr.WIDTH, gr.HEIGHT)
RADIAL_BAND_HALF_WIDTH = (gr.DEGREE_RING_SPACING - COMPONENT_DIAGONAL) / 4.0


def component_adjacency(circuit):
    adjacency = {ref: set() for ref in circuit.refs}
    for net in circuit.nets:
        refs = sorted({terminal.rsplit(".", 1)[0] for terminal in net.terminals})
        if len(refs) <= 1:
            continue
        for index, first in enumerate(refs):
            for second in refs[index + 1:]:
                adjacency[first].add(second)
                adjacency[second].add(first)
    return {ref: tuple(sorted(neighbors)) for ref, neighbors in adjacency.items()}


def shortest_paths(adjacency, source):
    distances = {source: 0}
    queue = collections.deque([source])
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current]:
            if neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    return distances


def graph_analysis(circuit):
    adjacency = component_adjacency(circuit)
    all_distances = {ref: shortest_paths(adjacency, ref) for ref in circuit.refs}
    count = len(circuit.refs)
    table = {}
    for ref in circuit.refs:
        distances = all_distances[ref]
        reachable = len(distances)
        distance_sum = sum(distances.values())
        if reachable <= 1 or distance_sum == 0:
            closeness = 0.0
        else:
            closeness = ((reachable - 1) / distance_sum) * ((reachable - 1) / (count - 1))
        table[ref] = {"degree": len(adjacency[ref]), "distance_sum": distance_sum,
                      "closeness": closeness, "reachable_components": reachable}
    root = min(circuit.refs, key=lambda ref: (-table[ref]["closeness"],
                                               table[ref]["distance_sum"],
                                               -table[ref]["degree"], ref))
    layers = all_distances[root]
    return adjacency, all_distances, table, root, layers, len(layers) == count


def layer_radii(layers):
    populations = collections.Counter(layers.values())
    radii = {0: 0.0}
    previous = 0.0
    for layer in range(1, max(populations) + 1):
        count = populations[layer]
        chord_radius = (gr.DEGREE_RING_SPACING / (2.0 * math.sin(math.pi / count))
                        if count > 1 else 0.0)
        radius = max(previous + gr.DEGREE_RING_SPACING, chord_radius)
        radii[layer] = radius
        previous = radius
    return radii


def clamp_to_band(point, current, nominal_radius):
    low = nominal_radius - RADIAL_BAND_HALF_WIDTH
    high = nominal_radius + RADIAL_BAND_HALF_WIDTH
    radius = math.hypot(*point)
    if radius < 1e-12:
        current_radius = math.hypot(*current)
        direction = ((current[0] / current_radius, current[1] / current_radius)
                     if current_radius > 0 else (1.0, 0.0))
    else:
        direction = (point[0] / radius, point[1] / radius)
    constrained_radius = max(low, min(high, radius))
    return direction[0] * constrained_radius, direction[1] * constrained_radius


def signed_angle_difference(first, second):
    return (second - first + math.pi) % (2.0 * math.pi) - math.pi


def resolve_layer_spacing(positions, layers):
    maximum_move = 0.0
    for _ in range(OVERLAP_RESOLUTION_PASSES):
        changed = False
        for layer in range(1, max(layers.values()) + 1):
            refs = sorted(ref for ref, value in layers.items() if value == layer)
            for index, first in enumerate(refs):
                for second in refs[index + 1:]:
                    p1, p2 = positions[first], positions[second]
                    r1, r2 = math.hypot(*p1), math.hypot(*p2)
                    cosine = max(-1.0, min(1.0,
                        (r1 * r1 + r2 * r2 - gr.DEGREE_RING_SPACING ** 2) / (2.0 * r1 * r2)))
                    required = math.acos(cosine)
                    angle1, angle2 = math.atan2(p1[1], p1[0]), math.atan2(p2[1], p2[0])
                    difference = signed_angle_difference(angle1, angle2)
                    if abs(difference) + 1e-12 >= required:
                        continue
                    sign = 1.0 if difference >= 0.0 else -1.0
                    correction = (required - abs(difference)) / 2.0
                    new_angle1 = angle1 - sign * correction
                    new_angle2 = angle2 + sign * correction
                    new1 = (r1 * math.cos(new_angle1), r1 * math.sin(new_angle1))
                    new2 = (r2 * math.cos(new_angle2), r2 * math.sin(new_angle2))
                    maximum_move = max(maximum_move, math.dist(p1, new1), math.dist(p2, new2))
                    positions[first], positions[second] = new1, new2
                    changed = True
        if not changed:
            break
    return maximum_move


def closeness_initial_positions(circuit):
    adjacency, distances, table, root, layers, connected = graph_analysis(circuit)
    if not connected:
        raise ValueError("IAMP component graph is disconnected; ROOT cannot layer all components")
    radii = layer_radii(layers)
    positions = {root: (0.0, 0.0)}
    for layer in range(1, max(layers.values()) + 1):
        refs = sorted(ref for ref, value in layers.items() if value == layer)
        for index, ref in enumerate(refs):
            angle = -math.pi / 2.0 + 2.0 * math.pi * index / len(refs)
            positions[ref] = (radii[layer] * math.cos(angle), radii[layer] * math.sin(angle))
    sweeps = 0
    for sweep in range(1, INITIALIZATION_SWEEPS + 1):
        maximum_move = 0.0
        for ref in sorted(component for component in circuit.refs if component != root):
            desired = (sum(positions[n][0] for n in adjacency[ref]) / len(adjacency[ref]),
                       sum(positions[n][1] for n in adjacency[ref]) / len(adjacency[ref]))
            constrained = clamp_to_band(desired, positions[ref], radii[layers[ref]])
            maximum_move = max(maximum_move, math.dist(positions[ref], constrained))
            positions[ref] = constrained
        maximum_move = max(maximum_move, resolve_layer_spacing(positions, layers))
        positions[root] = (0.0, 0.0)
        sweeps = sweep
        if maximum_move < MOVEMENT_TOLERANCE:
            break
    return positions, {"adjacency": adjacency, "all_distances": distances,
                       "centrality": table, "root": root, "layers": layers,
                       "connected": connected, "layer_radii": radii,
                       "initialization_sweeps": sweeps}


def release_initialization(positions):
    """Return optimizer input containing coordinates and no graph constraints."""
    return dict(positions)


def percentage_change(final, baseline):
    return None if baseline == 0 else gr.rounded(100.0 * (final - baseline) / baseline)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    circuit = gr.parse_circuit(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    initial, analysis = closeness_initial_positions(circuit)
    initial_metrics = gr.evaluate(circuit, initial)
    gr.write_svg(args.output_dir / "closeness_initial.svg", circuit, initial)
    graph_document = {
        "connected": analysis["connected"], "adjacency": analysis["adjacency"],
        "degree": {ref: analysis["centrality"][ref]["degree"] for ref in circuit.refs},
        "all_pairs_shortest_path_distances": analysis["all_distances"],
        "distance_sums": {ref: analysis["centrality"][ref]["distance_sum"] for ref in circuit.refs},
        "closeness": {ref: analysis["centrality"][ref]["closeness"] for ref in circuit.refs},
        "root": analysis["root"], "layers": analysis["layers"],
        "layer_radii": analysis["layer_radii"],
        "initialization_sweeps": analysis["initialization_sweeps"],
        "initialization_sweep_limit": INITIALIZATION_SWEEPS,
        "initialization_coordinates": {ref: {"x": gr.rounded(point[0]),
                                               "y": gr.rounded(point[1])}
                                       for ref, point in initial.items()},
    }
    gr.write_json(args.output_dir / "closeness_graph.json", graph_document)
    optimizer_start = release_initialization(initial)
    best, stats = gr.anneal(circuit, optimizer_start)
    final, quench_sweeps, reason = gr.relax(circuit, best)
    final_metrics = gr.evaluate(circuit, final)
    gr.write_svg(args.output_dir / "closeness_annealed.svg", circuit, final)
    coordinates_path = args.output_dir / "closeness_coordinates.json"
    gr.write_json(coordinates_path, gr.coordinate_document(circuit, final))
    baseline_one = json.loads((args.output_dir / "metrics.json").read_text(encoding="utf-8"))
    baseline_two = json.loads((args.output_dir / "annealed_metrics.json").read_text(encoding="utf-8"))
    baseline_three = json.loads((args.output_dir / "junction_annealed_metrics.json").read_text(encoding="utf-8"))
    baseline_four = json.loads((args.output_dir / "degree_seeded_metrics.json").read_text(encoding="utf-8"))
    four_a = baseline_four["states"]["experiment_4a_degree_seeded_annealed_and_quenched"]
    changes = {
        "total_connection_length": {"absolute": gr.rounded(final_metrics["total_connection_length"] - four_a["total_connection_length"]),
                                    "percent": percentage_change(final_metrics["total_connection_length"], four_a["total_connection_length"])},
        "component_overlaps": {"absolute": final_metrics["component_overlaps"] - four_a["component_overlaps"],
                               "percent": percentage_change(final_metrics["component_overlaps"], four_a["component_overlaps"])},
        "clearance_violations": {"absolute": final_metrics["clearance_violations"] - four_a["clearance_violations"],
                                 "percent": percentage_change(final_metrics["clearance_violations"], four_a["clearance_violations"])},
        "net_crossings": {"absolute": final_metrics["net_crossings"] - four_a["net_crossings"],
                          "percent": percentage_change(final_metrics["net_crossings"], four_a["net_crossings"])},
        "crowding": {"absolute": gr.rounded(final_metrics["energy"]["crowding"] - four_a["energy"]["crowding"]),
                     "percent": percentage_change(final_metrics["energy"]["crowding"], four_a["energy"]["crowding"])},
        "objective": {"absolute": gr.rounded(final_metrics["objective"] - four_a["objective"]),
                      "percent": percentage_change(final_metrics["objective"], four_a["objective"])},
    }
    centrality_table = [{"component": ref, "degree": analysis["centrality"][ref]["degree"],
                         "distance_sum": analysis["centrality"][ref]["distance_sum"],
                         "closeness": gr.rounded(analysis["centrality"][ref]["closeness"]),
                         "layer": analysis["layers"][ref]}
                        for ref in sorted(circuit.refs)]
    report = {
        "root": analysis["root"], "root_metrics": next(
            row for row in centrality_table if row["component"] == analysis["root"]),
        "maximum_layer": max(analysis["layers"].values()),
        "centrality_table": centrality_table,
        "closeness_initial": gr.rounded_metrics(initial_metrics),
        "states": {
            "original_initial": baseline_one["initial"],
            "experiment_1_final": baseline_one["final"],
            "experiment_2_annealed_and_quenched": baseline_two["states"]["experiment_2_annealed_and_quenched"],
            "experiment_3_junction_annealed_and_quenched": baseline_three["states"]["experiment_3_junction_annealed_and_quenched"],
            "experiment_4a_degree_seeded_annealed_and_quenched": four_a,
            "experiment_4c_closeness_annealed_and_quenched": gr.rounded_metrics(final_metrics),
        },
        "changes_from_experiment_4a": changes,
        "annealing": {"seed": gr.ANNEAL_SEED, "initial_temperature": gr.INITIAL_TEMPERATURE,
                      "final_temperature": gr.FINAL_TEMPERATURE, "cooling_schedule": "geometric",
                      "initial_move_distance": gr.INITIAL_MOVE_DISTANCE,
                      "final_move_distance": gr.FINAL_MOVE_DISTANCE,
                      "move_distance_schedule": "geometric", "iterations": gr.ANNEAL_ITERATIONS,
                      "accepted_downhill_moves": stats["accepted_downhill_moves"],
                      "accepted_uphill_moves": stats["accepted_uphill_moves"],
                      "rejected_moves": stats["rejected_moves"],
                      "best_energy_encountered": gr.rounded(stats["best_energy"])},
        "quench": {"final_energy": gr.rounded(final_metrics["objective"]),
                   "sweeps": quench_sweeps, "stopping_reason": reason},
        "coordinate_sha256": hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
    }
    gr.write_json(args.output_dir / "closeness_metrics.json", report)
    print(json.dumps(report, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
