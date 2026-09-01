#!/usr/bin/env python3
"""Four cumulative directional selections with radial centroid pulls."""

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
MOTION_FRAMES = 12
PAUSE_FRAMES = 5


def extents(positions, center, theta):
    radians = math.radians(theta)
    unit = (math.cos(radians), math.sin(radians))
    projections = {vertex: (positions[vertex][0] - center[0]) * unit[0] +
                               (positions[vertex][1] - center[1]) * unit[1]
                   for vertex in sorted(positions)}
    return projections, max(projections.values()), -min(projections.values())


def markup(positions, edges, center, attached, theta):
    attached = set(attached)
    lines = ['<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g>')
    for vertex in sorted(positions):
        color = "#f57c00" if vertex in attached else "#1565c0"
        radius = 5 if vertex in attached else 3
        lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="{radius}" fill="{color}"><title>{vertex}</title></circle>')
    lines.extend([f'<circle cx="{center[0]:.6f}" cy="{center[1]:.6f}" r="7" fill="#d32f2f" stroke="white" stroke-width="2"/>',
                  f'<text x="{center[0]+10:.6f}" y="{center[1]-10:.6f}" font-family="monospace" font-size="14" fill="#d32f2f">C</text>'])
    if theta is not None:
        radians = math.radians(theta)
        endpoint = (center[0] + 180 * math.cos(radians), center[1] + 180 * math.sin(radians))
        lines.extend(['<defs><marker id="radial-arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#f57c00"/></marker></defs>',
                      f'<line x1="{center[0]:.6f}" y1="{center[1]:.6f}" x2="{endpoint[0]:.6f}" y2="{endpoint[1]:.6f}" stroke="#f57c00" stroke-width="3" marker-end="url(#radial-arrow)"/>'])
    return lines


def write_static(path, positions, edges, center, viewbox, title):
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>{title}</title>'] + markup(positions, edges, center, (), None) + ['</svg>']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_motion(path, frames, edges, center, viewbox):
    count = len(frames)
    duration = count * 0.13
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Cumulative four-direction radial balancing</title>', '<style>']
    epsilon = 0.000001
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = (f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} '
                    f'{end+epsilon:.7f}%,100%{{opacity:0}}')
        lines.append(f'@keyframes q{index} {{{rule}}}')
        lines.append(f'#q{index}{{opacity:0;animation:q{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="q{index}">')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">{frame["label"]}</text>')
        lines.extend(markup(frame["positions"], edges, center, frame["attached"], frame["theta"])); lines.append('</g>')
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
    validation = balance.validate(graph, embedding, current, original_4k, edges)
    if not validation["abstract_graph_valid"]:
        raise RuntimeError("starting abstract graph is invalid")
    records = []
    frames = [{"label": "START", "positions": starting, "attached": (), "theta": None} for _ in range(PAUSE_FRAMES)]
    previous_end = {vertex: current[vertex] for vertex in sorted(current)}
    for theta in DIRECTIONS:
        operation_start = {vertex: current[vertex] for vertex in sorted(current)}
        if operation_start != previous_end:
            raise RuntimeError(f"cumulative continuity failure before theta={theta}")
        projections, forward, opposite = extents(operation_start, center, theta)
        selected = [vertex for vertex in sorted(operation_start) if projections[vertex] > 0]
        if forward <= opposite:
            attached, target, accepted, limited, limiting = [], 1.0, 1.0, False, None
            current = operation_start
        else:
            attached = selected
            target = opposite / forward
            accepted, current, full, limiting = balance.continuous_pull(
                graph, embedding, original_4k, edges, operation_start, center, attached, target)
            limited = not full
        for step in range(1, MOTION_FRAMES + 1):
            scale = 1.0 + step / MOTION_FRAMES * (accepted - 1.0)
            frames.append({"label": f"theta={theta}° radial pull, s={scale:.6f}",
                           "positions": balance.pulled_positions(operation_start, center, attached, scale),
                           "attached": attached, "theta": theta})
        previous_end = {vertex: current[vertex] for vertex in sorted(current)}
        if frames[-1]["positions"] != previous_end:
            raise RuntimeError(f"animation endpoint differs from accepted theta={theta} coordinates")
        for _ in range(PAUSE_FRAMES):
            frames.append({"label": f"theta={theta}° accepted — pause", "positions": previous_end,
                           "attached": (), "theta": theta})
        _, forward_after, opposite_after = extents(current, center, theta)
        records.append({"theta": theta, "selected_half_plane_nodes": len(selected),
                        "attached_nodes": len(attached), "R_forward_before": gr.rounded(forward),
                        "R_opposite_before": gr.rounded(opposite), "requested_s": gr.rounded(target),
                        "accepted_s": gr.rounded(accepted), "validity_limited": limited,
                        "limiting_event": limiting, "R_forward_after": gr.rounded(forward_after),
                        "R_opposite_after": gr.rounded(opposite_after),
                        "cumulative_input_matches_previous_output": True})
    final_validation = balance.validate(graph, embedding, current, original_4k, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("final radial4 graph is invalid")
    viewbox = balance.viewport(starting)
    start_path = output_dir / "directional_balance_radial4_start.svg"
    end_path = output_dir / "directional_balance_radial4_end.svg"
    motion_path = output_dir / "directional_balance_radial4_motion.svg"
    report_path = output_dir / "directional_balance_radial4_report.json"
    write_static(start_path, starting, edges, center, viewbox, "Radial4 start")
    write_static(end_path, current, edges, center, viewbox, "Radial4 end")
    write_motion(motion_path, frames, edges, center, viewbox)
    start_metrics, end_metrics = balance.measures(starting, edges), balance.measures(current, edges)
    gr.write_json(report_path, {"source": "output/directional_balance_start.svg",
                                "frozen_centroid": abstract_graph.point_doc(center),
                                "final_centroid": abstract_graph.point_doc(balance.centroid(current)),
                                "directions": records, "start": start_metrics, "end": end_metrics,
                                "final_validation": final_validation,
                                "all_cumulative_continuity_assertions_passed": True,
                                "animation_frame_count": len(frames),
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1],
                                                   "width": viewbox[2], "height": viewbox[3]}})
    print("theta | selected | attached | R_forward | R_opposite | requested s | accepted s | limited | R_forward after | R_opposite after")
    for record in records:
        print(f'{record["theta"]:>5} | {record["selected_half_plane_nodes"]:>8} | {record["attached_nodes"]:>8} | '
              f'{record["R_forward_before"]:.6f} | {record["R_opposite_before"]:.6f} | '
              f'{record["requested_s"]:.6f} | {record["accepted_s"]:.6f} | '
              f'{str(record["validity_limited"]):>7} | {record["R_forward_after"]:.6f} | {record["R_opposite_after"]:.6f}')
    print(json.dumps({"start": start_metrics, "end": end_metrics,
                      "frozen_centroid": abstract_graph.point_doc(center),
                      "final_centroid": abstract_graph.point_doc(balance.centroid(current)),
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
