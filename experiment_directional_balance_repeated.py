#!/usr/bin/env python3
"""Repeat the unchanged directional balancing sweep to its natural fixed point."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_directional_balance as balance
import graph_relax as gr


MAX_SWEEPS = 20


def load_svg_vertices(path, expected_vertices):
    root = ET.parse(path).getroot()
    positions = {}
    for circle in root.iter("{http://www.w3.org/2000/svg}circle"):
        title = circle.find("{http://www.w3.org/2000/svg}title")
        if title is not None and title.text in expected_vertices:
            positions[title.text] = (float(circle.attrib["cx"]), float(circle.attrib["cy"]))
    if set(positions) != set(expected_vertices):
        missing = sorted(set(expected_vertices) - set(positions))
        raise ValueError(f"directional_balance_end.svg does not contain the exact graph vertices: {missing}")
    return {vertex: positions[vertex] for vertex in sorted(positions)}


def write_progress_animation(path, frames, edges, viewbox):
    count = len(frames)
    duration = max(1.0, count * 0.65)
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Repeated directional balancing sweeps</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        lines.append(f'@keyframes r{index} {{0%,{start:.7f}%{{opacity:0}} {start:.7f}%,{end:.7f}%{{opacity:1}} {end:.7f}%,100%{{opacity:0}}}}')
        lines.append(f'#r{index}{{opacity:0;animation:r{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="r{index}">')
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">{frame["label"]}</text>')
        lines.extend(balance.graph_markup(frame["positions"], edges, frame["centroid"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, original_4k, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    source_svg = output_dir / "directional_balance_end.svg"
    current = load_svg_vertices(source_svg, original_4k)
    initial_validation = balance.validate(graph, embedding, current, original_4k, edges)
    if not initial_validation["abstract_graph_valid"]:
        raise RuntimeError("requested starting SVG is not a valid abstract graph")
    start_positions = {vertex: current[vertex] for vertex in sorted(current)}
    fixed_viewport = balance.viewport(start_positions)
    animation_frames = [{"label": "repeated-sweep start", "positions": start_positions,
                         "centroid": balance.centroid(start_positions)}]
    sweep_records = []
    sweep_paths = []
    termination = None
    for sweep_number in range(1, MAX_SWEEPS + 1):
        sweep_start = {vertex: current[vertex] for vertex in sorted(current)}
        sweep_start_metrics = balance.measures(sweep_start, edges)
        sweep_centroid = balance.centroid(sweep_start)
        accepted_pulls = 0
        direction_records = []
        for theta in balance.DIRECTIONS:
            radians = math.radians(theta)
            direction = (math.cos(radians), math.sin(radians))
            projections = {vertex: (current[vertex][0] - sweep_centroid[0]) * direction[0] +
                                       (current[vertex][1] - sweep_centroid[1]) * direction[1]
                           for vertex in sorted(current)}
            rplus, rminus = max(projections.values()), -min(projections.values())
            if rplus <= rminus:
                attached, target, accepted, full, limiting = [], 1.0, 1.0, True, None
            else:
                attached = [vertex for vertex in sorted(current) if projections[vertex] > rminus]
                target = rminus / rplus
                accepted, candidate, full, limiting = balance.continuous_pull(
                    graph, embedding, original_4k, edges, current, sweep_centroid, attached, target)
                if accepted < 1.0:
                    current = candidate
                    accepted_pulls += 1
                    animation_frames.append({"label": f"sweep {sweep_number}, theta={theta}°",
                                             "positions": current, "centroid": sweep_centroid})
            direction_records.append({"theta": theta, "Rplus_before": gr.rounded(rplus),
                                      "Rminus_before": gr.rounded(rminus),
                                      "attached_vertices": len(attached), "requested_s_target": gr.rounded(target),
                                      "accepted_s": gr.rounded(accepted), "full_balance_reached": full,
                                      "limiting_validity_event": limiting})
        final_validation = balance.validate(graph, embedding, current, original_4k, edges)
        if not final_validation["abstract_graph_valid"]:
            raise RuntimeError(f"sweep {sweep_number} ended in an invalid abstract graph")
        end_metrics = balance.measures(current, edges)
        area_change = 100 * (end_metrics["area"] - sweep_start_metrics["area"]) / sweep_start_metrics["area"]
        record = {"sweep": sweep_number, "centroid_used": abstract_graph.point_doc(sweep_centroid),
                  "final_centroid": abstract_graph.point_doc(balance.centroid(current)),
                  "accepted_pulls": accepted_pulls, "width": end_metrics["width"],
                  "height": end_metrics["height"], "area": end_metrics["area"],
                  "total_edge_length": end_metrics["total_edge_length"],
                  "area_change_percent": gr.rounded(area_change), "directions": direction_records}
        sweep_records.append(record)
        sweep_path = output_dir / f"directional_balance_repeated_sweep_{sweep_number:02d}.svg"
        balance.write_static(sweep_path, current, edges, sweep_centroid, fixed_viewport,
                             f"Repeated directional balance after sweep {sweep_number}")
        sweep_paths.append(str(sweep_path))
        animation_frames.append({"label": f"sweep {sweep_number} complete", "positions": current,
                                 "centroid": sweep_centroid})
        if accepted_pulls == 0:
            termination = "zero accepted pulls"
            break
    if termination is None:
        termination = "20-sweep safety limit"
    end_path = output_dir / "directional_balance_repeated_end.svg"
    motion_path = output_dir / "directional_balance_repeated_motion.svg"
    report_path = output_dir / "directional_balance_repeated_report.json"
    balance.write_static(end_path, current, edges, balance.centroid(current), fixed_viewport,
                         "Repeated directional balance final")
    write_progress_animation(motion_path, animation_frames, edges, fixed_viewport)
    gr.write_json(report_path, {"source": "output/directional_balance_end.svg",
                                "algorithm": "unchanged directional balancing sweep",
                                "maximum_sweeps": MAX_SWEEPS, "termination": termination,
                                "starting_centroid": abstract_graph.point_doc(balance.centroid(start_positions)),
                                "final_centroid": abstract_graph.point_doc(balance.centroid(current)),
                                "sweeps": sweep_records, "sweep_svgs": sweep_paths,
                                "animation_frame_count": len(animation_frames),
                                "fixed_viewport": {"x": fixed_viewport[0], "y": fixed_viewport[1],
                                                   "width": fixed_viewport[2], "height": fixed_viewport[3]}})
    print("sweep | centroid used | accepted pulls | width | height | area | edge length | area change %")
    for record in sweep_records:
        center = record["centroid_used"]
        print(f'{record["sweep"]:>5} | ({center["x"]:.6f}, {center["y"]:.6f}) | '
              f'{record["accepted_pulls"]:>14} | {record["width"]:.6f} | {record["height"]:.6f} | '
              f'{record["area"]:.6f} | {record["total_edge_length"]:.6f} | {record["area_change_percent"]:.6f}')
    print(f"termination: {termination}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.input, args.output_dir)


if __name__ == "__main__":
    main()
