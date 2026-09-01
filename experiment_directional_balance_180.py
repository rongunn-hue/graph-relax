#!/usr/bin/env python3
"""One coarse four-direction balancing pass on the abstract graph."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_directional_balance as balance
import experiment_directional_balance_repeated as repeated
import graph_relax as gr


DIRECTIONS = (0, 90, 180, 270)
ANIMATION_STEPS = 10


def directional_positions(current, center, attached, scale, theta):
    result = {vertex: current[vertex] for vertex in sorted(current)}
    horizontal = theta in (0, 180)
    for vertex in attached:
        x, y = current[vertex]
        if horizontal:
            result[vertex] = (center[0] + scale * (x - center[0]), y)
        else:
            result[vertex] = (x, center[1] + scale * (y - center[1]))
    return result


def directional_continuous_pull(graph, embedding, original, edges, current, center, attached, target, theta):
    last_t = 0.0
    last_positions = current
    first_invalid_t = None
    first_invalid_validation = None
    for step in range(1, balance.PATH_SCAN_STEPS + 1):
        t = step / balance.PATH_SCAN_STEPS
        scale = 1.0 + t * (target - 1.0)
        candidate = directional_positions(current, center, attached, scale, theta)
        candidate_validation = balance.validate(graph, embedding, candidate, original, edges)
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
    while high - low > balance.BISECTION_TOLERANCE and iterations < balance.MAX_BISECTION_ITERATIONS:
        middle = (low + high) / 2
        scale = 1.0 + middle * (target - 1.0)
        candidate = directional_positions(current, center, attached, scale, theta)
        candidate_validation = balance.validate(graph, embedding, candidate, original, edges)
        iterations += 1
        if candidate_validation["abstract_graph_valid"]:
            low, low_positions = middle, candidate
        else:
            high, first_invalid_validation = middle, candidate_validation
    return 1.0 + low * (target - 1.0), low_positions, False, balance.failure_reason(first_invalid_validation)


def directional_extents(positions, center, theta):
    radians = math.radians(theta)
    direction = (math.cos(radians), math.sin(radians))
    projections = {vertex: (positions[vertex][0] - center[0]) * direction[0] +
                               (positions[vertex][1] - center[1]) * direction[1]
                   for vertex in sorted(positions)}
    return projections, max(projections.values()), -min(projections.values())


def half_plane_markup(center, theta, viewbox):
    left, top, width, height = viewbox
    right, bottom = left + width, top + height
    if theta == 0:
        x, y, w, h = center[0], top, right - center[0], height
    elif theta == 90:
        x, y, w, h = left, center[1], width, bottom - center[1]
    elif theta == 180:
        x, y, w, h = left, top, center[0] - left, height
    else:
        x, y, w, h = left, top, width, center[1] - top
    return f'<rect x="{x:.6f}" y="{y:.6f}" width="{w:.6f}" height="{h:.6f}" fill="#ffb300" opacity="0.12"/>'


def write_motion(path, frames, edges, center, viewbox):
    count = len(frames)
    duration = count * 0.16
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Coarse 180-degree directional balancing</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        epsilon = 0.000001
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = (f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} '
                    f'{end+epsilon:.7f}%,100%{{opacity:0}}')
        lines.append(f'@keyframes c{index} {{{rule}}}')
        lines.append(f'#c{index}{{opacity:0;animation:c{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="c{index}">{half_plane_markup(center, frame["theta"], viewbox)}')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">theta={frame["theta"]}°, s={frame["scale"]:.6f}</text>')
        lines.extend(balance.graph_markup(frame["positions"], edges, center, frame["theta"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, original_4k, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    source = output_dir / "directional_balance_start.svg"
    current = repeated.load_svg_vertices(source, original_4k)
    starting = {vertex: current[vertex] for vertex in sorted(current)}
    center = balance.centroid(starting)
    initial_validation = balance.validate(graph, embedding, starting, original_4k, edges)
    if not initial_validation["abstract_graph_valid"]:
        raise RuntimeError("directional_balance_start.svg is not a valid abstract graph")
    start_metrics = balance.measures(starting, edges)
    records = []
    stage_records = [{"stage": "START", **start_metrics, "attached_vertices": 0,
                      "requested_s": 1.0, "accepted_s": 1.0}]
    animation_frames = [{"theta": 0, "scale": 1.0, "positions": {v: current[v] for v in sorted(current)}}]
    previous_accepted = {vertex: current[vertex] for vertex in sorted(current)}
    for theta in DIRECTIONS:
        if current != previous_accepted:
            raise RuntimeError(f"state propagation regression before theta={theta}: current coordinates differ from preceding accepted state")
        before = {vertex: current[vertex] for vertex in sorted(current)}
        projections, forward, opposite = directional_extents(before, center, theta)
        imbalance = forward - opposite
        forward_vertices = [vertex for vertex in sorted(before) if projections[vertex] > 0]
        if forward <= opposite:
            attached, target, accepted, full, limiting = [], 1.0, 1.0, True, None
            current = {vertex: before[vertex] for vertex in sorted(before)}
        else:
            attached = forward_vertices
            target = opposite / forward
            accepted, current, full, limiting = directional_continuous_pull(
                graph, embedding, original_4k, edges, before, center, attached, target, theta)
        for step in range(1, ANIMATION_STEPS + 1):
            interpolated_scale = 1.0 + step / ANIMATION_STEPS * (accepted - 1.0)
            frame_positions = directional_positions(before, center, attached, interpolated_scale, theta)
            if theta in (0, 180):
                if any(frame_positions[vertex][1] != before[vertex][1] for vertex in before):
                    raise RuntimeError(f"horizontal operation theta={theta} changed a Y coordinate")
            else:
                if any(frame_positions[vertex][0] != before[vertex][0] for vertex in before):
                    raise RuntimeError(f"vertical operation theta={theta} changed an X coordinate")
            animation_frames.append({"theta": theta, "scale": interpolated_scale,
                                     "positions": frame_positions})
        previous_accepted = {vertex: current[vertex] for vertex in sorted(current)}
        _, forward_after, opposite_after = directional_extents(current, center, theta)
        records.append({"theta": theta, "R_forward_before": gr.rounded(forward),
                        "R_opposite_before": gr.rounded(opposite), "imbalance": gr.rounded(imbalance),
                        "forward_half_plane_vertices": len(forward_vertices),
                        "attached_vertices": len(attached), "requested_s_target": gr.rounded(target),
                        "accepted_s": gr.rounded(accepted), "validity_limited": not full,
                        "limiting_validity_event": limiting, "R_forward_after": gr.rounded(forward_after),
                        "R_opposite_after": gr.rounded(opposite_after)})
        stage_records.append({"stage": f"{theta}°", **balance.measures(current, edges),
                              "attached_vertices": len(attached), "requested_s": gr.rounded(target),
                              "accepted_s": gr.rounded(accepted)})
    final_validation = balance.validate(graph, embedding, current, original_4k, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("coarse directional result is invalid")
    end_metrics = balance.measures(current, edges)
    final_centroid = balance.centroid(current)
    viewbox = balance.viewport(starting)
    start_path = output_dir / "directional_balance_180_start.svg"
    end_path = output_dir / "directional_balance_180_end.svg"
    motion_path = output_dir / "directional_balance_180_motion.svg"
    report_path = output_dir / "directional_balance_180_report.json"
    balance.write_static(start_path, starting, edges, center, viewbox, "Coarse directional balance start")
    balance.write_static(end_path, current, edges, center, viewbox, "Coarse directional balance end")
    write_motion(motion_path, animation_frames, edges, center, viewbox)
    serialized_final = repeated.load_svg_vertices(end_path, original_4k)
    expected_serialized_final = {vertex: (float(f'{current[vertex][0]:.6f}'), float(f'{current[vertex][1]:.6f}'))
                                 for vertex in sorted(current)}
    if serialized_final != expected_serialized_final:
        raise RuntimeError("final SVG coordinates do not equal the fourth cumulative accepted state")
    gr.write_json(report_path, {"source": "output/directional_balance_start.svg",
                                "frozen_centroid": abstract_graph.point_doc(center),
                                "final_centroid": abstract_graph.point_doc(final_centroid),
                                "directions": records, "stages": stage_records,
                                "state_propagation_regression_check": "passed",
                                "final_svg_coordinate_check": "passed",
                                "start": start_metrics, "end": end_metrics,
                                "final_validation": final_validation, "animation_frames": len(animation_frames),
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1],
                                                   "width": viewbox[2], "height": viewbox[3]}})
    print("theta | forward nodes | attached | R_forward | R_opposite | requested s | accepted s | limited")
    for record in records:
        print(f'{record["theta"]:>5} | {record["forward_half_plane_vertices"]:>13} | '
              f'{record["attached_vertices"]:>8} | {record["R_forward_before"]:.6f} | '
              f'{record["R_opposite_before"]:.6f} | {record["requested_s_target"]:.6f} | '
              f'{record["accepted_s"]:.6f} | {str(record["validity_limited"]):>7}')
    print("stage | width | height | area | edge length | attached | requested s | accepted s")
    for stage in stage_records:
        print(f'{stage["stage"]:>5} | {stage["width"]:.6f} | {stage["height"]:.6f} | '
              f'{stage["area"]:.6f} | {stage["total_edge_length"]:.6f} | '
              f'{stage["attached_vertices"]:>8} | {stage["requested_s"]:.6f} | {stage["accepted_s"]:.6f}')
    print(json.dumps({"frozen_centroid": abstract_graph.point_doc(center),
                      "final_centroid": abstract_graph.point_doc(final_centroid),
                      "start": start_metrics, "end": end_metrics,
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
