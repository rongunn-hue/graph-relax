#!/usr/bin/env python3
"""Experiment 4S: uniform single-pole contraction of an abstract planar graph."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


POLE = (4868.294118, -570.470588)
H_SCALE = 1823.809091
HEIGHT_RATIOS = (0.0, 0.10, 0.25, 0.50, 1.0, 2.0, 4.0, 9.0, 19.0, 99.0)
SVG_HEIGHT_RATIOS = (0.0, 0.50, 1.0, 2.0, 4.0, 9.0, 19.0, 99.0)
IDENTITY_TOLERANCE = 2e-8


def contraction_factor(height_ratio):
    return 1.0 / (1.0 + height_ratio)


def contracted_positions(initial, height_ratio):
    if height_ratio == 0.0:
        return {vertex: initial[vertex] for vertex in sorted(initial)}
    scale = contraction_factor(height_ratio)
    return {
        vertex: (POLE[0] + scale * (initial[vertex][0] - POLE[0]),
                 POLE[1] + scale * (initial[vertex][1] - POLE[1]))
        for vertex in sorted(initial)
    }


def pairwise_distances(positions):
    vertices = sorted(positions)
    return {(first, second): math.dist(positions[first], positions[second])
            for index, first in enumerate(vertices) for second in vertices[index + 1:]}


def extents(positions):
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    return max(xs) - min(xs), max(ys) - min(ys)


def evaluate(graph, embedding, initial, edges, baseline_pairs, baseline_metrics, height_ratio):
    scale = contraction_factor(height_ratio)
    positions = contracted_positions(initial, height_ratio)
    validation = abstract_graph.validate_abstract(graph, embedding, positions, initial, edges)
    width, height = extents(positions)
    area = width * height
    edge_length = sum(math.dist(positions[a], positions[b]) for a, b in edges)
    pairs = pairwise_distances(positions)
    pair_ratios = [pairs[pair] / baseline_pairs[pair] for pair in sorted(pairs)]
    maximum_pairwise = max(pairs.values())
    minimum_nonzero = min(distance for distance in pairs.values() if distance > 0)
    ratios = {
        "width": width / baseline_metrics["width"],
        "height": height / baseline_metrics["height"],
        "area": area / baseline_metrics["area"],
        "total_edge_length": edge_length / baseline_metrics["total_edge_length"],
    }
    errors = {
        "width_ratio_error": abs(ratios["width"] - scale),
        "height_ratio_error": abs(ratios["height"] - scale),
        "area_ratio_error": abs(ratios["area"] - scale * scale),
        "edge_length_ratio_error": abs(ratios["total_edge_length"] - scale),
        "maximum_pairwise_ratio_error": max(abs(ratio - scale) for ratio in pair_ratios),
    }
    identities_hold = max(errors.values()) <= IDENTITY_TOLERANCE
    if not identities_hold:
        raise RuntimeError(f"uniform-similarity identity failed at h/H={height_ratio}: {errors}")
    record = {
        "h_over_H": gr.rounded(height_ratio), "h": gr.rounded(height_ratio * H_SCALE),
        "contraction_factor": gr.rounded(scale),
        "width": gr.rounded(width), "height": gr.rounded(height), "area": gr.rounded(area),
        "total_edge_length": gr.rounded(edge_length),
        "maximum_pairwise_vertex_distance": gr.rounded(maximum_pairwise),
        "minimum_nonzero_pairwise_vertex_distance": gr.rounded(minimum_nonzero),
        "width_ratio": gr.rounded(ratios["width"]), "height_ratio": gr.rounded(ratios["height"]),
        "area_ratio": gr.rounded(ratios["area"]),
        "total_edge_length_ratio": gr.rounded(ratios["total_edge_length"]),
        "maximum_pairwise_ratio_error": gr.rounded(errors["maximum_pairwise_ratio_error"]),
        "similarity_identities_verified": identities_hold,
        "unrelated_edge_crossings": validation["unrelated_edge_crossings"],
        "coincident_vertex_pairs": validation["coincident_vertex_pairs"],
        "vertices_on_unrelated_edges": validation["vertices_on_unrelated_edges"],
        "cyclic_order_violations": validation["cyclic_order_violations"],
        "cyclic_order_orientation": validation["cyclic_order_orientation"],
        "abstract_graph_valid": validation["abstract_graph_valid"],
    }
    return record, positions, validation


def shared_viewport(initial):
    points = list(initial.values()) + [POLE]
    left, right = min(p[0] for p in points) - 80, max(p[0] for p in points) + 80
    top, bottom = min(p[1] for p in points) - 80, max(p[1] for p in points) + 80
    return left, top, right - left, bottom - top


def write_state_svg(path, positions, edges, viewport, height_ratio):
    left, top, width, height = viewport
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>Single-pole contraction h/H={height_ratio:g}</title>',
             '<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.extend(['</g>', f'<circle cx="{POLE[0]:.6f}" cy="{POLE[1]:.6f}" r="7" fill="#7b1fa2"/>',
                  f'<text x="{POLE[0]+10:.6f}" y="{POLE[1]+18:.6f}" font-family="monospace" font-size="15" fill="#7b1fa2">P</text>', '</svg>'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_overlay_svg(path, states, edges, viewport):
    left, top, width, height = viewport
    colors = ("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a", "#ff7c43", "#ffa600")
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Single-pole contraction overlay</title>']
    for (ratio, positions), color in zip(states, colors):
        lines.append(f'<g stroke="{color}" stroke-width="1" fill="{color}" opacity="0.46"><title>h/H={ratio:g}</title>')
        for first, second in edges:
            lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
        for vertex in sorted(positions):
            lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="2"/>')
        lines.append('</g>')
    lines.extend([f'<circle cx="{POLE[0]:.6f}" cy="{POLE[1]:.6f}" r="7" fill="#111"/>',
                  f'<text x="{POLE[0]+10:.6f}" y="{POLE[1]+18:.6f}" font-family="monospace" font-size="15">P</text>', '</svg>'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path, records):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def run(input_path, output_dir):
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    baseline_pairs = pairwise_distances(initial)
    baseline_width, baseline_height = extents(initial)
    baseline_metrics = {
        "width": baseline_width, "height": baseline_height,
        "area": baseline_width * baseline_height,
        "total_edge_length": sum(math.dist(initial[a], initial[b]) for a, b in edges),
    }
    sweep, payloads = [], {}
    for height_ratio in HEIGHT_RATIOS:
        payload = evaluate(graph, embedding, initial, edges, baseline_pairs, baseline_metrics, height_ratio)
        sweep.append(payload[0])
        payloads[height_ratio] = payload
    if not all(record["abstract_graph_valid"] and record["similarity_identities_verified"] for record in sweep):
        raise RuntimeError("a finite contraction state failed abstract validity or similarity identities")

    viewport = shared_viewport(initial)
    svg_paths = []
    overlay_states = []
    for height_ratio in SVG_HEIGHT_RATIOS:
        positions = payloads[height_ratio][1]
        name = f'planar_single_pole_{height_ratio:.2f}'.replace('.', '_') + '.svg'
        path = output_dir / name
        write_state_svg(path, positions, edges, viewport, height_ratio)
        svg_paths.append(path)
        overlay_states.append((height_ratio, positions))
    overlay_path = output_dir / "planar_single_pole_overlay.svg"
    write_overlay_svg(overlay_path, overlay_states, edges, viewport)

    coordinates_path = output_dir / "planar_single_pole_coordinates.json"
    gr.write_json(coordinates_path, {
        "model": "pure abstract graph uniform similarity contraction",
        "pole": abstract_graph.point_doc(POLE), "H": H_SCALE,
        "edges": [list(edge) for edge in edges],
        "states": [{"h_over_H": record["h_over_H"], "contraction_factor": record["contraction_factor"],
                    "vertices": {vertex: abstract_graph.point_doc(payloads[ratio][1][vertex]) for vertex in sorted(initial)}}
                   for ratio, record in zip(HEIGHT_RATIOS, sweep)],
    })
    metrics_path = output_dir / "planar_single_pole_metrics.json"
    gr.write_json(metrics_path, {"model": "abstract graph only", "physical_geometry_participated": False,
                                 "pole": abstract_graph.point_doc(POLE), "H": H_SCALE, "sweep": sweep})
    validation_path = output_dir / "planar_single_pole_validation.json"
    gr.write_json(validation_path, {"physical_validation_participated": False, "vertex_count": len(initial),
                                    "edge_count": len(edges), "all_states_valid": True,
                                    "all_similarity_identities_verified": True,
                                    "states": {str(ratio): payloads[ratio][2] for ratio in HEIGHT_RATIOS}})
    csv_path = output_dir / "planar_single_pole_sweep.csv"
    write_csv(csv_path, sweep)
    paths = [coordinates_path, metrics_path, validation_path, csv_path, overlay_path] + svg_paths
    return {"runtime_seconds": gr.rounded(time.perf_counter() - started), "sweep": sweep,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
