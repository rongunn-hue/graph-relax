#!/usr/bin/env python3
"""One cumulative four-direction radial envelope-balancing pass."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_directional_balance as balance
import experiment_directional_balance_radial4 as radial4
import experiment_directional_balance_repeated as repeated
import graph_relax as gr


DIRECTIONS = (0, 90, 180, 270)
MOTION_FRAMES = 12
PAUSE_FRAMES = 6


def write_motion(path, frames, edges, center, viewbox):
    count, epsilon = len(frames), 0.000001
    duration = count * 0.14
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Four-direction radial envelope balancing</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        lines.append(f'@keyframes e{index} {{{rule}}}')
        lines.append(f'#e{index}{{opacity:0;animation:e{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="e{index}">')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="20">{frame["label"]}</text>')
        lines.append(f'<text x="{left+25:.6f}" y="{top+62:.6f}" font-family="monospace" font-size="17">R_forward={frame["forward"]:.6f}, R_opposite={frame["opposite"]:.6f}</text>')
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
    if not balance.validate(graph, embedding, starting, original, edges)["abstract_graph_valid"]:
        raise RuntimeError("starting graph is invalid")
    records, frames = [], []
    previous_end = {vertex: current[vertex] for vertex in sorted(current)}
    for theta in DIRECTIONS:
        operation_start = {vertex: current[vertex] for vertex in sorted(current)}
        if operation_start != previous_end:
            raise RuntimeError(f"cumulative continuity failure before theta={theta}")
        projections, forward, opposite = radial4.extents(operation_start, center, theta)
        selected = [vertex for vertex in sorted(operation_start) if projections[vertex] > 0]
        imbalance = forward - opposite
        if forward > opposite:
            correction_required = True
            target = opposite / forward
            accepted, current, full, limiting = balance.continuous_pull(
                graph, embedding, original, edges, operation_start, center, selected, target)
            limited = not full
            for step in range(1, MOTION_FRAMES + 1):
                scale = 1.0 + step / MOTION_FRAMES * (accepted - 1.0)
                frames.append({"label": f"theta={theta}° radial correction, s={scale:.6f}",
                               "positions": balance.pulled_positions(operation_start, center, selected, scale),
                               "selected": selected, "theta": theta, "forward": forward, "opposite": opposite})
        else:
            correction_required = False
            target, accepted, limited, limiting = None, None, False, None
            current = operation_start
            for _ in range(MOTION_FRAMES):
                frames.append({"label": f"theta={theta}° — NO CORRECTION REQUIRED",
                               "positions": operation_start, "selected": selected, "theta": theta,
                               "forward": forward, "opposite": opposite})
        previous_end = {vertex: current[vertex] for vertex in sorted(current)}
        if frames[-1]["positions"] != previous_end:
            raise RuntimeError(f"animation endpoint mismatch at theta={theta}")
        _, forward_after, opposite_after = radial4.extents(current, center, theta)
        for _ in range(PAUSE_FRAMES):
            frames.append({"label": f"theta={theta}° accepted — pause", "positions": previous_end,
                           "selected": (), "theta": theta, "forward": forward_after, "opposite": opposite_after})
        stage_metrics = balance.measures(current, edges)
        records.append({"theta": theta, "selected_node_count": len(selected),
                        "R_forward_before": gr.rounded(forward), "R_opposite_before": gr.rounded(opposite),
                        "imbalance": gr.rounded(imbalance), "correction_required": correction_required,
                        "requested_s_target": None if target is None else gr.rounded(target),
                        "accepted_s": None if accepted is None else gr.rounded(accepted),
                        "validity_limited": limited, "limiting_event": limiting,
                        "R_forward_after": gr.rounded(forward_after), "R_opposite_after": gr.rounded(opposite_after),
                        **stage_metrics, "cumulative_input_matches_previous_output": True})
    final_validation = balance.validate(graph, embedding, current, original, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("final graph is invalid")
    viewbox = balance.viewport(starting)
    start_path = output_dir / "radial4_envelope_balance_start.svg"
    end_path = output_dir / "radial4_envelope_balance_end.svg"
    motion_path = output_dir / "radial4_envelope_balance_motion.svg"
    report_path = output_dir / "radial4_envelope_balance_report.json"
    radial4.write_static(start_path, starting, edges, center, viewbox, "Envelope balance start")
    radial4.write_static(end_path, current, edges, center, viewbox, "Envelope balance end")
    write_motion(motion_path, frames, edges, center, viewbox)
    start_metrics, final_metrics = balance.measures(starting, edges), balance.measures(current, edges)
    gr.write_json(report_path, {"model": "cumulative four-direction radial envelope balance",
                                "physical_geometry_participated": False,
                                "frozen_centroid": abstract_graph.point_doc(center),
                                "final_centroid": abstract_graph.point_doc(balance.centroid(current)),
                                "directions": records, "start": start_metrics, "final": final_metrics,
                                "final_validation": final_validation,
                                "all_cumulative_continuity_assertions_passed": True,
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1],
                                                   "width": viewbox[2], "height": viewbox[3]}})
    print("theta | selected | R_forward | R_opposite | imbalance | correction | requested s | accepted s | limited | R_forward after | R_opposite after | width | height | area | edge length")
    for r in records:
        print(f'{r["theta"]:>5} | {r["selected_node_count"]:>8} | {r["R_forward_before"]:.6f} | '
              f'{r["R_opposite_before"]:.6f} | {r["imbalance"]:.6f} | {str(r["correction_required"]):>10} | '
              f'{str(r["requested_s_target"]):>11} | {str(r["accepted_s"]):>10} | {str(r["validity_limited"]):>7} | '
              f'{r["R_forward_after"]:.6f} | {r["R_opposite_after"]:.6f} | {r["width"]:.6f} | '
              f'{r["height"]:.6f} | {r["area"]:.6f} | {r["total_edge_length"]:.6f}')
    print(json.dumps({"start": start_metrics, "final": final_metrics,
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
