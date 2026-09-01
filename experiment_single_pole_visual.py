#!/usr/bin/env python3
"""Visual uniform contraction of the abstract 4K graph toward a pole vertex."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


FIXED_SCALES = (1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
TOLERANCE = 2e-8


def select_pole(initial):
    vertices = sorted(initial)
    nearest = {
        vertex: min(math.dist(initial[vertex], initial[other]) for other in vertices if other != vertex)
        for vertex in vertices
    }
    pole_vertex = sorted(vertices, key=lambda vertex: (-nearest[vertex], vertex))[0]
    return pole_vertex, initial[pole_vertex], nearest[pole_vertex], nearest


def positions_at_scale(initial, pole_vertex, pole, scale):
    positions = {}
    for vertex in sorted(initial):
        if vertex == pole_vertex:
            positions[vertex] = pole
        else:
            positions[vertex] = (pole[0] + scale * (initial[vertex][0] - pole[0]),
                                  pole[1] + scale * (initial[vertex][1] - pole[1]))
    return positions


def extents(positions):
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    return max(xs) - min(xs), max(ys) - min(ys)


def evaluate(graph, embedding, initial, edges, pole_vertex, pole, scale, baseline):
    positions = positions_at_scale(initial, pole_vertex, pole, scale)
    validation = abstract_graph.validate_abstract(graph, embedding, positions, initial, edges)
    width, height = extents(positions)
    area = width * height
    edge_length = sum(math.dist(positions[first], positions[second]) for first, second in edges)
    nearest_to_pole = min(math.dist(pole, positions[vertex]) for vertex in positions if vertex != pole_vertex)
    ratios = {"width": width / baseline["width"], "height": height / baseline["height"],
              "area": area / baseline["area"], "edge_length": edge_length / baseline["edge_length"]}
    errors = (abs(ratios["width"] - scale), abs(ratios["height"] - scale),
              abs(ratios["area"] - scale * scale), abs(ratios["edge_length"] - scale))
    if max(errors) > TOLERANCE:
        raise RuntimeError(f"similarity identity failed at s={scale}: {errors}")
    record = {
        "s": gr.rounded(scale), "width": gr.rounded(width), "height": gr.rounded(height),
        "area": gr.rounded(area), "total_edge_length": gr.rounded(edge_length),
        "nearest_node_distance_to_pole": gr.rounded(nearest_to_pole),
        "width_ratio": gr.rounded(ratios["width"]), "height_ratio": gr.rounded(ratios["height"]),
        "area_ratio": gr.rounded(ratios["area"]), "edge_length_ratio": gr.rounded(ratios["edge_length"]),
        "unrelated_edge_crossings": validation["unrelated_edge_crossings"],
        "coincident_vertex_pairs": validation["coincident_vertex_pairs"],
        "vertices_on_unrelated_edges": validation["vertices_on_unrelated_edges"],
        "cyclic_order_violations": validation["cyclic_order_violations"],
        "abstract_graph_valid": validation["abstract_graph_valid"],
        "similarity_identities_verified": True,
    }
    return record, positions, validation


def shared_viewport(initial):
    xs = [point[0] for point in initial.values()]
    ys = [point[1] for point in initial.values()]
    margin = 80.0
    return min(xs) - margin, min(ys) - margin, max(xs) - min(xs) + 2 * margin, max(ys) - min(ys) + 2 * margin


def graph_group(positions, edges, pole_vertex, pole, opacity=None):
    opacity_attribute = "" if opacity is None else f' opacity="{opacity}"'
    lines = [f'<g{opacity_attribute}><g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        if vertex != pole_vertex:
            lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.extend(['</g>', f'<circle cx="{pole[0]:.6f}" cy="{pole[1]:.6f}" r="8" fill="#d32f2f" stroke="white" stroke-width="2"><title>Pole: {pole_vertex}</title></circle>',
                  f'<text x="{pole[0]+11:.6f}" y="{pole[1]-11:.6f}" font-family="monospace" font-size="14" fill="#d32f2f">P</text>', '</g>'])
    return lines


def write_svg(path, positions, edges, pole_vertex, pole, viewport, scale, label):
    left, top, width, height = viewport
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>Single-pole contraction {label}, s={scale:.9f}</title>']
    lines.extend(graph_group(positions, edges, pole_vertex, pole))
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_animation(path, frames, edges, pole_vertex, pole, viewport):
    sequence = list(frames) + list(reversed(frames[:-1]))
    duration = len(sequence) * 0.55
    left, top, width, height = viewport
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Animated single-pole contraction</title>']
    for index, (label, scale, positions) in enumerate(sequence):
        lines.append(f'<g visibility="hidden">')
        lines.extend(graph_group(positions, edges, pole_vertex, pole))
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="22" fill="#111">{label}: s={scale:.6f}</text>')
        lines.append(f'<animate attributeName="visibility" values="visible;hidden" keyTimes="0;1" begin="{index*0.55:.2f}s" dur="0.55s" repeatCount="indefinite" restart="always"/>')
        lines.append('</g>')
    # A single discrete visibility timeline provides broad SVG-viewer compatibility.
    lines.append(f'<set attributeName="data-duration" to="{duration:.2f}"/>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_animation_css(path, frames, edges, pole_vertex, pole, viewport):
    sequence = list(frames) + list(reversed(frames[:-1]))
    count = len(sequence)
    left, top, width, height = viewport
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Animated single-pole contraction</title>', '<style>']
    for index in range(count):
        start = 100 * index / count
        end = 100 * (index + 1) / count
        lines.append(f'@keyframes frame{index} {{ 0%, {start:.6f}% {{opacity:0}} {start:.6f}%, {end:.6f}% {{opacity:1}} {end:.6f}%, 100% {{opacity:0}} }}')
        lines.append(f'#frame{index} {{opacity:0; animation:frame{index} {count*0.55:.2f}s steps(1,end) infinite;}}')
    lines.append('</style>')
    for index, (label, scale, positions) in enumerate(sequence):
        lines.append(f'<g id="frame{index}">')
        lines.extend(graph_group(positions, edges, pole_vertex, pole))
        lines.append(f'<text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="22" fill="#111">{label}: s={scale:.6f}</text></g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    pole_vertex, pole, minimum_distance, nearest = select_pole(initial)
    stop_scale = 1.0 / minimum_distance
    scales = FIXED_SCALES + (stop_scale,)
    base_width, base_height = extents(initial)
    baseline = {"width": base_width, "height": base_height, "area": base_width * base_height,
                "edge_length": sum(math.dist(initial[a], initial[b]) for a, b in edges)}
    records, payloads = [], []
    viewport = shared_viewport(initial)
    frame_paths = []
    for index, scale in enumerate(scales):
        record, positions, validation = evaluate(graph, embedding, initial, edges, pole_vertex, pole, scale, baseline)
        if not record["abstract_graph_valid"]:
            raise RuntimeError(f"abstract validation failed at scale {scale}")
        records.append(record)
        label = "stop" if index == len(scales) - 1 else f'{scale:.2f}'
        payloads.append((label, scale, positions))
        filename = "single_pole_s_stop.svg" if label == "stop" else f'single_pole_s_{scale:.2f}'.replace('.', '_') + '.svg'
        path = output_dir / filename
        write_svg(path, positions, edges, pole_vertex, pole, viewport, scale, label)
        frame_paths.append(path)
    animation_path = output_dir / "single_pole_contraction.svg"
    write_animation_css(animation_path, payloads, edges, pole_vertex, pole, viewport)
    report_path = output_dir / "single_pole_visual_report.json"
    gr.write_json(report_path, {"model": "abstract point-and-segment graph", "physical_geometry_participated": False,
                                "pole_vertex": pole_vertex, "pole": abstract_graph.point_doc(pole),
                                "Lmin": gr.rounded(minimum_distance), "s_stop": gr.rounded(stop_scale),
                                "nearest_distance_by_vertex": {v: gr.rounded(nearest[v]) for v in sorted(nearest)},
                                "shared_viewport": {"x": viewport[0], "y": viewport[1], "width": viewport[2], "height": viewport[3]},
                                "frames": records})
    csv_path = output_dir / "single_pole_visual_metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    return {"runtime_seconds": gr.rounded(time.perf_counter() - started), "pole_vertex": pole_vertex,
            "pole": abstract_graph.point_doc(pole), "Lmin": gr.rounded(minimum_distance),
            "s_stop": gr.rounded(stop_scale), "frames": records,
            "artifacts": [str(path) for path in frame_paths + [animation_path, report_path, csv_path]]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
