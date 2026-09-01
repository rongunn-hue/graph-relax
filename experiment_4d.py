#!/usr/bin/env python3
"""Experiment 4D: neighborhood-affinity 2-D initialization."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import experiment_4c as e4c
import graph_relax as gr


INITIALIZATION_SWEEPS = 100
MOVEMENT_TOLERANCE = 1e-6


def affinity_values(adjacency, root):
    nonroot = sorted(ref for ref in adjacency if ref != root)
    values = {ref: {} for ref in nonroot}
    for index, first in enumerate(nonroot):
        first_neighbors = set(adjacency[first]) - {root}
        for second in nonroot[index + 1:]:
            second_neighbors = set(adjacency[second]) - {root}
            direct = int(second in adjacency[first])
            shared = len(first_neighbors & second_neighbors)
            affinity = direct + shared
            if affinity:
                values[first][second] = affinity
                values[second][first] = affinity
    return values


def provisional_positions(circuit, root):
    refs = sorted(ref for ref in circuit.refs if ref != root)
    radius = gr.DEGREE_RING_SPACING / (2.0 * math.sin(math.pi / len(refs)))
    positions = {root: (0.0, 0.0)}
    for index, ref in enumerate(refs):
        angle = -math.pi / 2.0 + 2.0 * math.pi * index / len(refs)
        positions[ref] = (radius * math.cos(angle), radius * math.sin(angle))
    return positions, radius


def clamp_radial(point, current, low, high):
    radius = math.hypot(*point)
    if radius < 1e-12:
        old_radius = math.hypot(*current)
        unit = ((current[0] / old_radius, current[1] / old_radius)
                if old_radius else (1.0, 0.0))
    else:
        unit = (point[0] / radius, point[1] / radius)
    radius = max(low, min(high, radius))
    return unit[0] * radius, unit[1] * radius


def affinity_proposals(current, affinities, root, low, high):
    """Compute all proposals from the same current-state snapshot."""
    proposals = {root: (0.0, 0.0)}
    for ref in sorted(affinities):
        related = affinities[ref]
        if related:
            total = sum(related.values())
            desired = (sum(weight * current[other][0] for other, weight in related.items()) / total,
                       sum(weight * current[other][1] for other, weight in related.items()) / total)
            proposals[ref] = clamp_radial(desired, current[ref], low, high)
        else:
            proposals[ref] = current[ref]
    return proposals


def guarantee_nonoverlap_ring(positions, root, radius):
    """Project a conflicted cyclic order onto its nearest evenly spaced ring."""
    refs = sorted((ref for ref in positions if ref != root),
                  key=lambda ref: (math.atan2(positions[ref][1], positions[ref][0]), ref))
    count = len(refs)
    angles = [math.atan2(positions[ref][1], positions[ref][0]) for ref in refs]
    gaps = [angles[i + 1] - angles[i] for i in range(count - 1)] + \
           [angles[0] + 2.0 * math.pi - angles[-1]]
    cut = (max(range(count), key=lambda index: (gaps[index], -index)) + 1) % count
    refs = refs[cut:] + refs[:cut]
    unwrapped = []
    previous = None
    for ref in refs:
        angle = math.atan2(positions[ref][1], positions[ref][0])
        if previous is not None:
            while angle <= previous:
                angle += 2.0 * math.pi
        unwrapped.append(angle)
        previous = angle
    separation = 2.0 * math.pi / count
    offset = sum(angle - index * separation for index, angle in enumerate(unwrapped)) / count
    maximum_move = 0.0
    for index, ref in enumerate(refs):
        angle = offset + index * separation
        point = (radius * math.cos(angle), radius * math.sin(angle))
        maximum_move = max(maximum_move, math.dist(positions[ref], point))
        positions[ref] = point
    return maximum_move


def affinity_initial_positions(circuit):
    adjacency, distances, centrality, root, layers, connected = e4c.graph_analysis(circuit)
    affinities = affinity_values(adjacency, root)
    provisional, radius = provisional_positions(circuit, root)
    band_low, band_high = radius - gr.DEGREE_RING_SPACING / 2.0, radius + gr.DEGREE_RING_SPACING / 2.0
    positions = dict(provisional)
    artificial_layers = {ref: (0 if ref == root else 1) for ref in circuit.refs}
    sweeps = 0
    stopping_reason = "initialization sweep limit reached"
    for sweep in range(1, INITIALIZATION_SWEEPS + 1):
        proposals = affinity_proposals(positions, affinities, root, band_low, band_high)
        maximum_move = max(math.dist(positions[ref], proposals[ref]) for ref in circuit.refs)
        positions = proposals
        maximum_move = max(maximum_move, e4c.resolve_layer_spacing(positions, artificial_layers))
        refs = sorted(ref for ref in circuit.refs if ref != root)
        if any(math.dist(positions[first], positions[second]) < gr.DEGREE_RING_SPACING - 1e-9
               for index, first in enumerate(refs) for second in refs[index + 1:]):
            maximum_move = max(maximum_move,
                               guarantee_nonoverlap_ring(positions, root, radius))
        positions[root] = (0.0, 0.0)
        sweeps = sweep
        if maximum_move < MOVEMENT_TOLERANCE:
            stopping_reason = "no meaningful initialization movement"
            break
    return positions, provisional, {
        "adjacency": adjacency, "all_distances": distances, "centrality": centrality,
        "root": root, "layers": layers, "connected": connected,
        "affinities": affinities, "provisional_radius": radius,
        "radial_band": {"minimum": band_low, "maximum": band_high},
        "initialization_sweeps": sweeps, "stopping_reason": stopping_reason,
    }


def release_initialization(positions):
    return dict(positions)


def change(final, baseline):
    absolute = final - baseline
    return {"absolute": gr.rounded(absolute) if isinstance(absolute, float) else absolute,
            "percent": None if baseline == 0 else gr.rounded(100.0 * absolute / baseline)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    circuit = gr.parse_circuit(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    initial, provisional, analysis = affinity_initial_positions(circuit)
    initial_metrics = gr.evaluate(circuit, initial)
    gr.write_svg(args.output_dir / "affinity_initial.svg", circuit, initial)
    nonzero_affinities = []
    for first in sorted(analysis["affinities"]):
        for second, affinity in sorted(analysis["affinities"][first].items()):
            if first < second:
                first_neighbors = set(analysis["adjacency"][first]) - {analysis["root"]}
                second_neighbors = set(analysis["adjacency"][second]) - {analysis["root"]}
                nonzero_affinities.append({
                    "a": first, "b": second,
                    "direct": int(second in analysis["adjacency"][first]),
                    "shared": len(first_neighbors & second_neighbors), "affinity": affinity})
    graph_document = {
        "connected": analysis["connected"], "adjacency": analysis["adjacency"],
        "centrality": analysis["centrality"], "root": analysis["root"],
        "nonzero_affinities": nonzero_affinities,
        "provisional_radius": analysis["provisional_radius"],
        "radial_band": analysis["radial_band"],
        "initialization_sweep_limit": INITIALIZATION_SWEEPS,
        "initialization_sweeps": analysis["initialization_sweeps"],
        "initialization_stopping_reason": analysis["stopping_reason"],
        "provisional_coordinates": {ref: {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}
                                    for ref, point in provisional.items()},
        "affinity_initialization_coordinates": {
            ref: {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}
            for ref, point in initial.items()},
    }
    gr.write_json(args.output_dir / "affinity_graph.json", graph_document)
    best, stats = gr.anneal(circuit, release_initialization(initial))
    final, quench_sweeps, quench_reason = gr.relax(circuit, best)
    final_metrics = gr.evaluate(circuit, final)
    gr.write_svg(args.output_dir / "affinity_annealed.svg", circuit, final)
    coordinates_path = args.output_dir / "affinity_coordinates.json"
    gr.write_json(coordinates_path, gr.coordinate_document(circuit, final))
    one = json.loads((args.output_dir / "metrics.json").read_text())
    two = json.loads((args.output_dir / "annealed_metrics.json").read_text())
    three = json.loads((args.output_dir / "junction_annealed_metrics.json").read_text())
    four_a_report = json.loads((args.output_dir / "degree_seeded_metrics.json").read_text())
    four_c_report = json.loads((args.output_dir / "closeness_metrics.json").read_text())
    four_a = four_a_report["states"]["experiment_4a_degree_seeded_annealed_and_quenched"]
    four_c = four_c_report["states"]["experiment_4c_closeness_annealed_and_quenched"]
    changes = {
        "total_connection_length": change(final_metrics["total_connection_length"], four_c["total_connection_length"]),
        "component_overlaps": change(final_metrics["component_overlaps"], four_c["component_overlaps"]),
        "clearance_violations": change(final_metrics["clearance_violations"], four_c["clearance_violations"]),
        "net_crossings": change(final_metrics["net_crossings"], four_c["net_crossings"]),
        "crowding": change(final_metrics["energy"]["crowding"], four_c["energy"]["crowding"]),
        "objective": change(final_metrics["objective"], four_c["objective"]),
    }
    report = {
        "root": analysis["root"],
        "initialization": {"sweep_limit": INITIALIZATION_SWEEPS,
                           "sweeps": analysis["initialization_sweeps"],
                           "stopping_reason": analysis["stopping_reason"],
                           "provisional_radius": gr.rounded(analysis["provisional_radius"]),
                           "radial_band": {key: gr.rounded(value) for key, value in analysis["radial_band"].items()}},
        "affinity_initial": gr.rounded_metrics(initial_metrics),
        "states": {
            "original_initial": one["initial"], "experiment_1_final": one["final"],
            "experiment_2_annealed_and_quenched": two["states"]["experiment_2_annealed_and_quenched"],
            "experiment_3_junction_annealed_and_quenched": three["states"]["experiment_3_junction_annealed_and_quenched"],
            "experiment_4a_degree_seeded_annealed_and_quenched": four_a,
            "experiment_4c_closeness_annealed_and_quenched": four_c,
            "experiment_4d_affinity_annealed_and_quenched": gr.rounded_metrics(final_metrics)},
        "changes_from_experiment_4c": changes,
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
                   "sweeps": quench_sweeps, "stopping_reason": quench_reason},
        "coordinate_sha256": hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
    }
    gr.write_json(args.output_dir / "affinity_metrics.json", report)
    print(json.dumps(report, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
