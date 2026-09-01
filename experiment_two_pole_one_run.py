#!/usr/bin/env python3
"""One fixed-pair, exact two-pole motion demonstration."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_4u as geometry
import graph_relax as gr


POLE_A = "component:RMIN1"
POLE_B = "net:GND"
FRAME_COUNT = 61
HEIGHT_BACKOFF = 1e-6


def markup(positions, edges, model):
    lines = ['<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        if vertex not in (POLE_A, POLE_B):
            lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.append('</g>')
    for label, vertex, color, offset in (("POLE A", POLE_A, "#d32f2f", 1), ("POLE B", POLE_B, "#7b1fa2", -1)):
        point = positions[vertex]
        lines.append(f'<circle cx="{point[0]:.6f}" cy="{point[1]:.6f}" r="9" fill="{color}" stroke="white" stroke-width="2"><title>{vertex}</title></circle>')
        lines.append(f'<text x="{point[0] + 14*offset:.6f}" y="{point[1] + 25:.6f}" text-anchor="{"start" if offset > 0 else "end"}" font-family="monospace" font-size="18" font-weight="bold" fill="{color}">{label}</text>')
    return lines


def write_static(path, positions, edges, model, viewbox, height, title):
    left, top, width, drawing_height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {drawing_height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{drawing_height:.6f}" fill="white"/>',
             f'<title>{title}</title>',
             f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="22" fill="#111">h = {height:.6f}</text>']
    lines.extend(markup(positions, edges, model)); lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_animation(path, frames, edges, model, viewbox):
    left, top, width, drawing_height = viewbox
    duration = FRAME_COUNT * 0.12
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {drawing_height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{drawing_height:.6f}" fill="white"/>',
             '<title>Two simultaneously rising poles</title>', '<style>']
    for index in range(FRAME_COUNT):
        start, end = 100 * index / FRAME_COUNT, 100 * (index + 1) / FRAME_COUNT
        lines.append(f'@keyframes f{index} {{0%,{start:.7f}%{{opacity:0}} {start:.7f}%,{end:.7f}%{{opacity:1}} {end:.7f}%,100%{{opacity:0}}}}')
        lines.append(f'#f{index}{{opacity:0;animation:f{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, (height, positions) in enumerate(frames):
        lines.append(f'<g id="f{index}">')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="22" fill="#111">common pole height h = {height:.6f}</text>')
        lines.extend(markup(positions, edges, model)); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    if POLE_A == POLE_B or initial[POLE_A] == initial[POLE_B]:
        raise RuntimeError("pole bases must be distinct existing graph nodes at distinct coordinates")
    model = geometry.pair_model(initial, POLE_A, POLE_B)
    first_degeneracy = model["theoretical_bound"]
    final_height = first_degeneracy - HEIGHT_BACKOFF
    frames = []
    for index in range(FRAME_COUNT):
        height = final_height * index / (FRAME_COUNT - 1)
        positions = geometry.positions_at_height(initial, model, height)
        validation = abstract_graph.validate_abstract(graph, embedding, positions, initial, edges)
        if not validation["abstract_graph_valid"]:
            raise RuntimeError(f"abstract graph became invalid before cable degeneracy at h={height}")
        frames.append((height, positions))
    viewbox = geometry.viewport(initial)
    start_path = output_dir / "two_pole_start.svg"
    end_path = output_dir / "two_pole_end.svg"
    motion_path = output_dir / "two_pole_motion.svg"
    write_static(start_path, frames[0][1], edges, model, viewbox, frames[0][0], "Two-pole start")
    write_static(end_path, frames[-1][1], edges, model, viewbox, frames[-1][0], "Two-pole end")
    write_animation(motion_path, frames, edges, model, viewbox)
    result = {"pole_a": POLE_A, "pole_a_coordinate": abstract_graph.point_doc(initial[POLE_A]),
              "pole_b": POLE_B, "pole_b_coordinate": abstract_graph.point_doc(initial[POLE_B]),
              "pole_separation": gr.rounded(model["distance"]), "frame_count": FRAME_COUNT,
              "first_degeneracy_height": gr.rounded(first_degeneracy),
              "limiting_vertex": model["bound_owner"], "final_rendered_height": gr.rounded(final_height),
              "fixed_viewport": {"x": viewbox[0], "y": viewbox[1], "width": viewbox[2], "height": viewbox[3]},
              "outputs": [str(start_path), str(end_path), str(motion_path)]}
    print(json.dumps(result, sort_keys=True, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.input, args.output_dir)


if __name__ == "__main__":
    main()
