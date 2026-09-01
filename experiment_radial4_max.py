#!/usr/bin/env python3
"""Maximum legal radial contraction of four cumulative selected half-planes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_directional_balance as balance
import experiment_directional_balance_radial4 as radial4
import experiment_directional_balance_repeated as repeated
import experiment_radial4_fixed as fixed
import graph_relax as gr


DIRECTIONS = (0, 90, 180, 270)
SCAN_RESOLUTION = 4096
BISECTION_TOLERANCE = 1e-9
MAX_BISECTION_ITERATIONS = 48
SAFETY_EPSILON = 1e-9
MOTION_FRAMES = 16
PAUSE_FRAMES = 6


def limiting_summary(validation):
    reasons = []
    details = {}
    if validation["unrelated_edge_crossings"]:
        reasons.append("unrelated edge crossing")
        details["crossing_edges"] = validation["crossing_details"]
    if validation["coincident_vertex_pairs"]:
        reasons.append("coincident vertices")
        details["coincident_vertices"] = validation["coincident_details"]
    if validation["vertices_on_unrelated_edges"]:
        reasons.append("vertex on unrelated edge")
        details["vertex_on_edge"] = validation["vertex_on_edge_details"]
    if validation["cyclic_order_violations"]:
        reasons.append("cyclic-order/embedding violation")
        details["cyclic_violation_vertices"] = validation["cyclic_violation_vertices"]
    return ", ".join(reasons) or "unknown abstract validity event", details


def maximum_pull(graph, embedding, original, edges, current, center, selected):
    last_valid_scale = 1.0
    last_valid_positions = current
    invalid_scale = None
    invalid_validation = None
    for step in range(1, SCAN_RESOLUTION + 1):
        scale = 1.0 - step / SCAN_RESOLUTION
        candidate = fixed.radial_positions(current, center, selected, scale)
        candidate_validation = fixed.valid(graph, embedding, candidate, original, edges)
        if not candidate_validation["abstract_graph_valid"]:
            invalid_scale, invalid_validation = scale, candidate_validation
            break
        last_valid_scale, last_valid_positions = scale, candidate
    if invalid_scale is None:
        accepted = SAFETY_EPSILON
        positions = fixed.radial_positions(current, center, selected, accepted)
        validation = fixed.valid(graph, embedding, positions, original, edges)
        if not validation["abstract_graph_valid"]:
            raise RuntimeError("safety-epsilon result unexpectedly invalid")
        return accepted, positions, "no earlier graph-validity event; numerical safety epsilon", {}, 0
    valid_scale, invalid = last_valid_scale, invalid_scale
    valid_positions = last_valid_positions
    iterations = 0
    while valid_scale - invalid > BISECTION_TOLERANCE and iterations < MAX_BISECTION_ITERATIONS:
        middle = (valid_scale + invalid) / 2
        candidate = fixed.radial_positions(current, center, selected, middle)
        candidate_validation = fixed.valid(graph, embedding, candidate, original, edges)
        iterations += 1
        if candidate_validation["abstract_graph_valid"]:
            valid_scale, valid_positions = middle, candidate
        else:
            invalid, invalid_validation = middle, candidate_validation
    reason, details = limiting_summary(invalid_validation)
    return valid_scale, valid_positions, reason, details, iterations


def write_motion(path, frames, edges, center, viewbox):
    count, epsilon = len(frames), 0.000001
    duration = count * 0.14
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Maximum four-direction radial compaction</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        lines.append(f'@keyframes m{index} {{{rule}}}')
        lines.append(f'#m{index}{{opacity:0;animation:m{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="m{index}"><text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">{frame["label"]}</text>')
        lines.extend(radial4.markup(frame["positions"], edges, center, frame["selected"], frame["theta"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, original, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    current = repeated.load_svg_vertices(output_dir / "directional_balance_start.svg", original)
    starting = {vertex: current[vertex] for vertex in sorted(current)}
    center = balance.centroid(starting)
    if not fixed.valid(graph, embedding, starting, original, edges)["abstract_graph_valid"]:
        raise RuntimeError("starting graph is invalid")
    records, stages = [], [{"stage": "START", **fixed.metrics(starting, edges)}]
    frames = [{"label": "START", "positions": starting, "selected": (), "theta": None} for _ in range(PAUSE_FRAMES)]
    previous_end = {vertex: current[vertex] for vertex in sorted(current)}
    for theta in DIRECTIONS:
        operation_start = {vertex: current[vertex] for vertex in sorted(current)}
        if operation_start != previous_end:
            raise RuntimeError(f"cumulative continuity failure before theta={theta}")
        selected = fixed.selected_vertices(operation_start, center, theta)
        accepted, current, reason, details, iterations = maximum_pull(
            graph, embedding, original, edges, operation_start, center, selected)
        for step in range(1, MOTION_FRAMES + 1):
            scale = 1.0 + step / MOTION_FRAMES * (accepted - 1.0)
            frames.append({"label": f"theta={theta}°, s={scale:.9f}",
                           "positions": fixed.radial_positions(operation_start, center, selected, scale),
                           "selected": selected, "theta": theta})
        previous_end = {vertex: current[vertex] for vertex in sorted(current)}
        if frames[-1]["positions"] != previous_end:
            raise RuntimeError(f"animation endpoint mismatch after theta={theta}")
        for _ in range(PAUSE_FRAMES):
            frames.append({"label": f"theta={theta}° maximum accepted — pause",
                           "positions": previous_end, "selected": (), "theta": theta})
        after = fixed.metrics(current, edges)
        stages.append({"stage": f"after {theta}°", **after})
        records.append({"theta": theta, "selected_nodes": len(selected), "starting_s": 1.0,
                        "accepted_minimum_valid_s": gr.rounded(accepted),
                        "radial_contraction": gr.rounded(1.0 - accepted),
                        "limiting_validity_event": reason, "limiting_details": details,
                        "scan_resolution": SCAN_RESOLUTION, "bisection_iterations": iterations,
                        "width_after": after["width"], "height_after": after["height"],
                        "area_after": after["area"], "total_edge_length_after": after["total_edge_length"],
                        "cumulative_input_matches_previous_output": True})
    final_validation = fixed.valid(graph, embedding, current, original, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("final graph is invalid")
    viewbox = fixed.viewport(starting)
    start_path = output_dir / "radial4_max_start.svg"
    end_path = output_dir / "radial4_max_end.svg"
    motion_path = output_dir / "radial4_max_motion.svg"
    report_path = output_dir / "radial4_max_report.json"
    fixed.write_static(start_path, starting, edges, center, viewbox, "Maximum radial4 start")
    fixed.write_static(end_path, current, edges, center, viewbox, "Maximum radial4 end")
    write_motion(motion_path, frames, edges, center, viewbox)
    start_metrics, final_metrics = stages[0], stages[-1]
    percent_changes = {key: gr.rounded(100 * (final_metrics[key] - start_metrics[key]) / start_metrics[key])
                       for key in ("width", "height", "area", "total_edge_length")}
    gr.write_json(report_path, {"model": "maximum cumulative four-direction radial compaction",
                                "physical_geometry_participated": False,
                                "frozen_centroid": abstract_graph.point_doc(center),
                                "scan_resolution": SCAN_RESOLUTION, "bisection_tolerance": BISECTION_TOLERANCE,
                                "directions": records, "stages": stages, "percent_change_start_to_final": percent_changes,
                                "final_validation": final_validation,
                                "all_cumulative_continuity_assertions_passed": True,
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1],
                                                   "width": viewbox[2], "height": viewbox[3]}})
    print("theta | selected | accepted s | contraction | limiting event | width | height | area | edge length")
    for r in records:
        print(f'{r["theta"]:>5} | {r["selected_nodes"]:>8} | {r["accepted_minimum_valid_s"]:.9f} | '
              f'{r["radial_contraction"]:.9f} | {r["limiting_validity_event"]} | '
              f'{r["width_after"]:.6f} | {r["height_after"]:.6f} | {r["area_after"]:.6f} | '
              f'{r["total_edge_length_after"]:.6f}')
    print(json.dumps({"stages": stages, "percent_change_start_to_final": percent_changes,
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
