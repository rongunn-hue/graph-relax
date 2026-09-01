#!/usr/bin/env python3
"""Experiment 4U: exact equal-height two-pole contraction of an abstract graph."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from itertools import combinations
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


SWEEP_STEPS = 64
BISECTION_TOLERANCE = 1e-6
MAX_BISECTION_ITERATIONS = 32
HEIGHT_BACKOFF = 1e-6
GEOMETRIC_TOLERANCE = 1e-10
VISUAL_FRACTIONS = (0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
                    0.70, 0.80, 0.90, 0.95, 1.00)


def pair_model(initial, pole_a, pole_b):
    a, b = initial[pole_a], initial[pole_b]
    distance = math.dist(a, b)
    ux, uy = (b[0] - a[0]) / distance, (b[1] - a[1]) / distance
    nx, ny = -uy, ux
    coordinates = {}
    for vertex in sorted(initial):
        dx, dy = initial[vertex][0] - a[0], initial[vertex][1] - a[1]
        coordinates[vertex] = (dx * ux + dy * uy, dx * nx + dy * ny)
    moving = [vertex for vertex in sorted(initial) if vertex not in (pole_a, pole_b)]
    theoretical_bound, owner = min((abs(coordinates[vertex][1]), vertex) for vertex in moving)
    return {"pole_a": pole_a, "pole_b": pole_b, "a": a, "b": b, "distance": distance,
            "u": (ux, uy), "n": (nx, ny), "coordinates": coordinates,
            "theoretical_bound": theoretical_bound, "bound_owner": owner}


def positions_at_height(initial, model, height):
    if height < 0:
        return None
    if height == 0:
        return {vertex: initial[vertex] for vertex in sorted(initial)}
    a = model["a"]
    ux, uy = model["u"]
    nx, ny = model["n"]
    result = {}
    for vertex in sorted(initial):
        if vertex in (model["pole_a"], model["pole_b"]):
            result[vertex] = initial[vertex]
            continue
        parallel, perpendicular = model["coordinates"][vertex]
        radicand = perpendicular * perpendicular - height * height
        if radicand <= GEOMETRIC_TOLERANCE:
            return None
        transformed_perpendicular = math.copysign(math.sqrt(radicand), perpendicular)
        result[vertex] = (a[0] + parallel * ux + transformed_perpendicular * nx,
                          a[1] + parallel * uy + transformed_perpendicular * ny)
    return result


def extents_and_distances(positions, edges):
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    width, height = max(xs) - min(xs), max(ys) - min(ys)
    vertices = sorted(positions)
    pairwise = [math.dist(positions[first], positions[second])
                for index, first in enumerate(vertices) for second in vertices[index + 1:]]
    return {"width": width, "height": height, "area": width * height,
            "edge_length": sum(math.dist(positions[first], positions[second]) for first, second in edges),
            "maximum_pairwise": max(pairwise),
            "minimum_nonzero_pairwise": min(value for value in pairwise if value > 0)}


def failure_condition(validation):
    if validation is None:
        return "two-circle geometric degeneracy"
    reasons = []
    if validation["unrelated_edge_crossings"]:
        reasons.append("unrelated edge crossing")
    if validation["coincident_vertex_pairs"]:
        reasons.append("coincident vertices")
    if validation["vertices_on_unrelated_edges"]:
        reasons.append("vertex on unrelated edge")
    if validation["cyclic_order_violations"]:
        reasons.append("embedding/cyclic-order violation")
    return ", ".join(reasons) or "unknown abstract-graph failure"


def evaluate_height(graph, embedding, initial, edges, model, height):
    positions = positions_at_height(initial, model, height)
    if positions is None:
        return False, None, None
    validation = abstract_graph.validate_abstract(graph, embedding, positions, initial, edges)
    return validation["abstract_graph_valid"], positions, validation


def search_pair(graph, embedding, initial, edges, baseline, pole_a, pole_b):
    model = pair_model(initial, pole_a, pole_b)
    bound = model["theoretical_bound"]
    if bound <= HEIGHT_BACKOFF:
        selected_height = 0.0
        positions = positions_at_height(initial, model, 0.0)
        validation = abstract_graph.validate_abstract(graph, embedding, positions, initial, edges)
        limiting_height, limiting_condition, limiting_vertex = 0.0, "two-circle geometric degeneracy at h=0", model["bound_owner"]
        bisection_iterations = 0
    else:
        low, low_positions, low_validation = 0.0, initial, abstract_graph.validate_abstract(graph, embedding, initial, initial, edges)
        high = bound
        limiting_condition = "two-circle geometric degeneracy"
        limiting_vertex = model["bound_owner"]
        for step in range(1, SWEEP_STEPS + 1):
            candidate = bound * step / SWEEP_STEPS
            valid, candidate_positions, candidate_validation = evaluate_height(graph, embedding, initial, edges, model, candidate)
            if not valid:
                high = candidate
                limiting_condition = failure_condition(candidate_validation)
                if candidate_validation is not None:
                    limiting_vertex = None
                break
            low, low_positions, low_validation = candidate, candidate_positions, candidate_validation
        bisection_iterations = 0
        while high - low > BISECTION_TOLERANCE and bisection_iterations < MAX_BISECTION_ITERATIONS:
            middle = (low + high) / 2
            valid, candidate_positions, candidate_validation = evaluate_height(graph, embedding, initial, edges, model, middle)
            bisection_iterations += 1
            if valid:
                low, low_positions, low_validation = middle, candidate_positions, candidate_validation
            else:
                high = middle
                limiting_condition = failure_condition(candidate_validation)
                if candidate_validation is not None:
                    limiting_vertex = None
        limiting_height = high
        selected_height = min(low, max(0.0, high - HEIGHT_BACKOFF))
        valid, positions, validation = evaluate_height(graph, embedding, initial, edges, model, selected_height)
        if not valid:
            selected_height, positions, validation = low, low_positions, low_validation
    measurements = extents_and_distances(positions, edges)
    record = {
        "pole_a": pole_a, "pole_b": pole_b, "pole_distance": gr.rounded(model["distance"]),
        "theoretical_upper_height": gr.rounded(bound),
        "greatest_valid_equal_height": gr.rounded(selected_height),
        "limiting_height": gr.rounded(limiting_height), "limiting_vertex": limiting_vertex,
        "limiting_condition": limiting_condition, "bisection_iterations": bisection_iterations,
        "width": gr.rounded(measurements["width"]), "height": gr.rounded(measurements["height"]),
        "area": gr.rounded(measurements["area"]), "total_edge_length": gr.rounded(measurements["edge_length"]),
        "maximum_pairwise_vertex_distance": gr.rounded(measurements["maximum_pairwise"]),
        "minimum_nonzero_pairwise_vertex_distance": gr.rounded(measurements["minimum_nonzero_pairwise"]),
        "unrelated_edge_crossings": validation["unrelated_edge_crossings"],
        "coincident_vertex_pairs": validation["coincident_vertex_pairs"],
        "vertices_on_unrelated_edges": validation["vertices_on_unrelated_edges"],
        "cyclic_order_violations": validation["cyclic_order_violations"],
        "cyclic_order_orientation": validation["cyclic_order_orientation"],
        "abstract_graph_valid": validation["abstract_graph_valid"],
        "width_ratio": gr.rounded(measurements["width"] / baseline["width"]),
        "height_ratio": gr.rounded(measurements["height"] / baseline["height"]),
        "area_ratio": gr.rounded(measurements["area"] / baseline["area"]),
        "edge_length_ratio": gr.rounded(measurements["edge_length"] / baseline["edge_length"]),
        "maximum_diameter_ratio": gr.rounded(measurements["maximum_pairwise"] / baseline["maximum_pairwise"]),
    }
    return record, model, positions, validation


def viewport(initial):
    xs, ys = [p[0] for p in initial.values()], [p[1] for p in initial.values()]
    return min(xs) - 80, min(ys) - 80, max(xs) - min(xs) + 160, max(ys) - min(ys) + 160


def graph_markup(positions, edges, model):
    lines = ['<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        if vertex not in (model["pole_a"], model["pole_b"]):
            lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.append('</g>')
    for label, vertex, color in (("A", model["pole_a"], "#d32f2f"), ("B", model["pole_b"], "#7b1fa2")):
        point = positions[vertex]
        lines.append(f'<circle cx="{point[0]:.6f}" cy="{point[1]:.6f}" r="8" fill="{color}" stroke="white" stroke-width="2"><title>{vertex}</title></circle>')
        lines.append(f'<text x="{point[0]+11:.6f}" y="{point[1]-11:.6f}" font-family="monospace" font-size="14" fill="{color}">{label}</text>')
    return lines


def write_svg(path, positions, edges, model, viewbox, title):
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>{title}</title>'] + graph_markup(positions, edges, model) + ['</svg>']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_animation(path, frames, edges, model, viewbox):
    sequence = frames + list(reversed(frames[:-1]))
    count = len(sequence)
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        lines.append(f'@keyframes u{index} {{0%,{start:.6f}%{{opacity:0}} {start:.6f}%,{end:.6f}%{{opacity:1}} {end:.6f}%,100%{{opacity:0}}}}')
        lines.append(f'#u{index}{{opacity:0;animation:u{index} {count*0.5:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, (label, height_value, positions) in enumerate(sequence):
        lines.append(f'<g id="u{index}">')
        lines.extend(graph_markup(positions, edges, model))
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="20">{label}: h={height_value:.6f}</text></g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    baseline = extents_and_distances(initial, edges)
    results = []
    payloads = {}
    for pole_a, pole_b in combinations(sorted(initial), 2):
        record, model, positions, validation = search_pair(graph, embedding, initial, edges, baseline, pole_a, pole_b)
        results.append(record)
        payloads[(pole_a, pole_b)] = (model, positions, validation)
    ranked = sorted(results, key=lambda r: (r["area_ratio"], r["edge_length_ratio"],
                                            r["maximum_diameter_ratio"], -r["greatest_valid_equal_height"],
                                            r["pole_a"], r["pole_b"]))
    best = ranked[0]
    model = payloads[(best["pole_a"], best["pole_b"])][0]
    selected_height = best["greatest_valid_equal_height"]
    viewbox = viewport(initial)
    frames = []
    frame_records = []
    frame_paths = []
    for index, fraction in enumerate(VISUAL_FRACTIONS):
        height_value = selected_height if fraction == 1.0 else fraction * selected_height
        valid, positions, validation = evaluate_height(graph, embedding, initial, edges, model, height_value)
        if not valid:
            raise RuntimeError("best-pair visual frame failed validation")
        label = "final" if fraction == 1.0 else f'{fraction:.2f}'
        filename = "two_pole_final.svg" if fraction == 1.0 else f'two_pole_fraction_{fraction:.2f}'.replace('.', '_') + '.svg'
        path = output_dir / filename
        write_svg(path, positions, edges, model, viewbox, f'Two-pole contraction {label}')
        frames.append((label, height_value, positions)); frame_paths.append(path)
        frame_records.append({"fraction": fraction, "height": gr.rounded(height_value),
                              **{key: gr.rounded(value) for key, value in extents_and_distances(positions, edges).items()},
                              "valid": True})
    animation_path = output_dir / "two_pole_contraction.svg"
    write_animation(animation_path, frames, edges, model, viewbox)
    results_path = output_dir / "two_pole_all_pairs.json"
    gr.write_json(results_path, {"pair_count": len(results), "sweep_steps": SWEEP_STEPS,
                                 "height_backoff": HEIGHT_BACKOFF, "ranking": ranked})
    report_path = output_dir / "two_pole_best_report.json"
    gr.write_json(report_path, {"model": "abstract exact equal-height two-circle contraction",
                                "physical_geometry_participated": False, "best_pair": best,
                                "shared_viewport": {"x": viewbox[0], "y": viewbox[1], "width": viewbox[2], "height": viewbox[3]},
                                "frames": frame_records})
    csv_path = output_dir / "two_pole_all_pairs.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(results)
    return {"runtime_seconds": gr.rounded(time.perf_counter() - started), "pair_count": len(results),
            "best_pair": best, "frame_records": frame_records,
            "artifacts": [str(path) for path in frame_paths + [animation_path, results_path, report_path, csv_path]]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
