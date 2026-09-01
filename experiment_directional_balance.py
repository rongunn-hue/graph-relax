#!/usr/bin/env python3
"""One directional single-centroid balancing sweep of the abstract 4K graph."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


DIRECTIONS = tuple(range(0, 360, 10))
PATH_SCAN_STEPS = 256
BISECTION_TOLERANCE = 1e-8
MAX_BISECTION_ITERATIONS = 40
ANIMATION_STEPS_PER_PULL = 8


def centroid(positions):
    vertices = sorted(positions)
    return (sum(positions[v][0] for v in vertices) / len(vertices),
            sum(positions[v][1] for v in vertices) / len(vertices))


def measures(positions, edges):
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    width, height = max(xs) - min(xs), max(ys) - min(ys)
    return {"width": gr.rounded(width), "height": gr.rounded(height),
            "area": gr.rounded(width * height),
            "total_edge_length": gr.rounded(sum(math.dist(positions[a], positions[b]) for a, b in edges)),
            "centroid": abstract_graph.point_doc(centroid(positions))}


def pulled_positions(current, fixed_center, attached, scale):
    result = {vertex: current[vertex] for vertex in sorted(current)}
    for vertex in attached:
        result[vertex] = (fixed_center[0] + scale * (current[vertex][0] - fixed_center[0]),
                          fixed_center[1] + scale * (current[vertex][1] - fixed_center[1]))
    return result


def failure_reason(validation):
    reasons = []
    if validation["unrelated_edge_crossings"]:
        reasons.append("unrelated edge crossing")
    if validation["coincident_vertex_pairs"]:
        reasons.append("coincident vertices")
    if validation["vertices_on_unrelated_edges"]:
        reasons.append("vertex on unrelated edge")
    if validation["cyclic_order_violations"]:
        reasons.append("cyclic-order/embedding violation")
    return ", ".join(reasons) or "unknown abstract validity event"


def validate(graph, embedding, positions, original, edges):
    return abstract_graph.validate_abstract(graph, embedding, positions, original, edges)


def continuous_pull(graph, embedding, original, edges, current, fixed_center, attached, target):
    last_t = 0.0
    last_positions = current
    first_invalid_t = None
    first_invalid_validation = None
    for step in range(1, PATH_SCAN_STEPS + 1):
        t = step / PATH_SCAN_STEPS
        scale = 1.0 + t * (target - 1.0)
        candidate = pulled_positions(current, fixed_center, attached, scale)
        candidate_validation = validate(graph, embedding, candidate, original, edges)
        if not candidate_validation["abstract_graph_valid"]:
            first_invalid_t = t
            first_invalid_validation = candidate_validation
            break
        last_t, last_positions = t, candidate
    if first_invalid_t is None:
        return target, last_positions, True, None
    low, high = last_t, first_invalid_t
    low_positions = last_positions
    iterations = 0
    while high - low > BISECTION_TOLERANCE and iterations < MAX_BISECTION_ITERATIONS:
        middle = (low + high) / 2
        scale = 1.0 + middle * (target - 1.0)
        candidate = pulled_positions(current, fixed_center, attached, scale)
        candidate_validation = validate(graph, embedding, candidate, original, edges)
        iterations += 1
        if candidate_validation["abstract_graph_valid"]:
            low, low_positions = middle, candidate
        else:
            high, first_invalid_validation = middle, candidate_validation
    accepted_scale = 1.0 + low * (target - 1.0)
    return accepted_scale, low_positions, False, failure_reason(first_invalid_validation)


def viewport(original):
    xs = [point[0] for point in original.values()]
    ys = [point[1] for point in original.values()]
    return min(xs) - 100, min(ys) - 100, max(xs) - min(xs) + 200, max(ys) - min(ys) + 200


def graph_markup(positions, edges, fixed_center, theta=None):
    lines = ['<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.extend(['</g>', f'<circle cx="{fixed_center[0]:.6f}" cy="{fixed_center[1]:.6f}" r="7" fill="#d32f2f" stroke="white" stroke-width="2"/>',
                  f'<text x="{fixed_center[0]+10:.6f}" y="{fixed_center[1]-10:.6f}" font-family="monospace" font-size="14" fill="#d32f2f">C</text>'])
    if theta is not None:
        radians = math.radians(theta)
        end = (fixed_center[0] + 180 * math.cos(radians), fixed_center[1] + 180 * math.sin(radians))
        lines.extend([f'<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#f57c00"/></marker></defs>',
                      f'<line x1="{fixed_center[0]:.6f}" y1="{fixed_center[1]:.6f}" x2="{end[0]:.6f}" y2="{end[1]:.6f}" stroke="#f57c00" stroke-width="3" marker-end="url(#arrow)"/>'])
    return lines


def write_static(path, positions, edges, fixed_center, viewbox, title):
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>{title}</title>'] + graph_markup(positions, edges, fixed_center) + ['</svg>']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_animation(path, frames, edges, fixed_center, viewbox):
    count = len(frames)
    duration = count * 0.10
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Directional centroid balancing sweep</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        lines.append(f'@keyframes d{index} {{0%,{start:.7f}%{{opacity:0}} {start:.7f}%,{end:.7f}%{{opacity:1}} {end:.7f}%,100%{{opacity:0}}}}')
        lines.append(f'#d{index}{{opacity:0;animation:d{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="d{index}">')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">theta={frame["theta"]}°, s={frame["scale"]:.6f}</text>')
        lines.extend(graph_markup(frame["positions"], edges, fixed_center, frame["theta"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, original, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    fixed_center = centroid(original)
    current = {vertex: original[vertex] for vertex in sorted(original)}
    initial_validation = validate(graph, embedding, current, original, edges)
    if not initial_validation["abstract_graph_valid"]:
        raise RuntimeError("original 4K abstract graph is invalid")
    start_measures = measures(current, edges)
    records = []
    animation_frames = [{"theta": 0, "scale": 1.0, "positions": current}]
    for theta in DIRECTIONS:
        radians = math.radians(theta)
        direction = (math.cos(radians), math.sin(radians))
        projections = {vertex: (current[vertex][0] - fixed_center[0]) * direction[0] +
                                   (current[vertex][1] - fixed_center[1]) * direction[1]
                       for vertex in sorted(current)}
        rplus, rminus = max(projections.values()), -min(projections.values())
        imbalance = rplus - rminus
        before = current
        if rplus <= rminus:
            attached = []
            target = 1.0
            accepted = 1.0
            full = True
            limiting = None
        else:
            attached = [vertex for vertex in sorted(current) if projections[vertex] > rminus]
            target = rminus / rplus
            accepted, current, full, limiting = continuous_pull(
                graph, embedding, original, edges, current, fixed_center, attached, target)
            if accepted < 1.0:
                for step in range(1, ANIMATION_STEPS_PER_PULL + 1):
                    scale = 1.0 + step / ANIMATION_STEPS_PER_PULL * (accepted - 1.0)
                    animation_frames.append({"theta": theta, "scale": scale,
                                             "positions": pulled_positions(before, fixed_center, attached, scale)})
        records.append({"theta": theta, "Rplus_before": gr.rounded(rplus),
                        "Rminus_before": gr.rounded(rminus), "imbalance": gr.rounded(imbalance),
                        "attached_vertices": len(attached), "requested_s_target": gr.rounded(target),
                        "accepted_s": gr.rounded(accepted), "full_balance_reached": full,
                        "limiting_validity_event": limiting})
    final_validation = validate(graph, embedding, current, original, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("final directional-balance graph is invalid")
    end_measures = measures(current, edges)
    final_centroid = centroid(current)
    viewbox = viewport(original)
    start_path = output_dir / "directional_balance_start.svg"
    end_path = output_dir / "directional_balance_end.svg"
    motion_path = output_dir / "directional_balance_motion.svg"
    write_static(start_path, original, edges, fixed_center, viewbox, "Directional balance start")
    write_static(end_path, current, edges, fixed_center, viewbox, "Directional balance end")
    write_animation(motion_path, animation_frames, edges, fixed_center, viewbox)
    report_path = output_dir / "directional_balance_report.json"
    gr.write_json(report_path, {"model": "abstract point-and-segment graph", "physical_geometry_participated": False,
                                "frozen_centroid": abstract_graph.point_doc(fixed_center),
                                "final_centroid": abstract_graph.point_doc(final_centroid),
                                "start": {**start_measures, "unrelated_crossings": initial_validation["unrelated_edge_crossings"],
                                          "coincident_vertex_pairs": initial_validation["coincident_vertex_pairs"],
                                          "vertices_on_unrelated_edges": initial_validation["vertices_on_unrelated_edges"]},
                                "end": {**end_measures, "unrelated_crossings": final_validation["unrelated_edge_crossings"],
                                        "coincident_vertex_pairs": final_validation["coincident_vertex_pairs"],
                                        "vertices_on_unrelated_edges": final_validation["vertices_on_unrelated_edges"]},
                                "directions": records, "animation_frame_count": len(animation_frames),
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1], "width": viewbox[2], "height": viewbox[3]}})
    print(json.dumps({"frozen_centroid": abstract_graph.point_doc(fixed_center),
                      "final_centroid": abstract_graph.point_doc(final_centroid),
                      "start": start_measures, "end": end_measures,
                      "accepted_directional_pulls": sum(record["accepted_s"] < 1.0 for record in records),
                      "animation_frames": len(animation_frames),
                      "outputs": [str(start_path), str(end_path), str(motion_path), str(report_path)]},
                     sort_keys=True, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.input, args.output_dir)


if __name__ == "__main__":
    main()
