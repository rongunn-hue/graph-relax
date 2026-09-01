#!/usr/bin/env python3
"""Four cumulative fixed-strength radial half-plane pulls on an abstract graph."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


DIRECTIONS = (0, 90, 180, 270)
REQUESTED_SCALE = 0.80
PATH_SCAN_STEPS = 256
BISECTION_TOLERANCE = 1e-8
MAX_BISECTION_ITERATIONS = 40
MOTION_FRAMES = 12
PAUSE_FRAMES = 5


def load_vertices(path, expected):
    root = ET.parse(path).getroot()
    positions = {}
    for circle in root.iter("{http://www.w3.org/2000/svg}circle"):
        title = circle.find("{http://www.w3.org/2000/svg}title")
        if title is not None and title.text in expected:
            positions[title.text] = (float(circle.attrib["cx"]), float(circle.attrib["cy"]))
    if set(positions) != set(expected):
        raise ValueError("starting SVG does not contain the exact abstract graph vertex set")
    return {vertex: positions[vertex] for vertex in sorted(positions)}


def centroid(positions):
    vertices = sorted(positions)
    return (sum(positions[v][0] for v in vertices) / len(vertices),
            sum(positions[v][1] for v in vertices) / len(vertices))


def selected_vertices(positions, center, theta):
    radians = math.radians(theta)
    unit = (math.cos(radians), math.sin(radians))
    return [vertex for vertex in sorted(positions)
            if (positions[vertex][0] - center[0]) * unit[0] +
               (positions[vertex][1] - center[1]) * unit[1] > 0]


def radial_positions(current, center, selected, scale):
    result = {vertex: current[vertex] for vertex in sorted(current)}
    for vertex in selected:
        result[vertex] = (center[0] + scale * (current[vertex][0] - center[0]),
                          center[1] + scale * (current[vertex][1] - center[1]))
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


def valid(graph, embedding, positions, original, edges):
    return abstract_graph.validate_abstract(graph, embedding, positions, original, edges)


def limited_pull(graph, embedding, original, edges, current, center, selected):
    last_t, last_positions = 0.0, current
    first_invalid_t = None
    invalid_validation = None
    for step in range(1, PATH_SCAN_STEPS + 1):
        t = step / PATH_SCAN_STEPS
        scale = 1.0 + t * (REQUESTED_SCALE - 1.0)
        candidate = radial_positions(current, center, selected, scale)
        candidate_validation = valid(graph, embedding, candidate, original, edges)
        if not candidate_validation["abstract_graph_valid"]:
            first_invalid_t, invalid_validation = t, candidate_validation
            break
        last_t, last_positions = t, candidate
    if first_invalid_t is None:
        return REQUESTED_SCALE, last_positions, False, None
    low, high = last_t, first_invalid_t
    low_positions = last_positions
    iterations = 0
    while high - low > BISECTION_TOLERANCE and iterations < MAX_BISECTION_ITERATIONS:
        middle = (low + high) / 2
        scale = 1.0 + middle * (REQUESTED_SCALE - 1.0)
        candidate = radial_positions(current, center, selected, scale)
        candidate_validation = valid(graph, embedding, candidate, original, edges)
        iterations += 1
        if candidate_validation["abstract_graph_valid"]:
            low, low_positions = middle, candidate
        else:
            high, invalid_validation = middle, candidate_validation
    return 1.0 + low * (REQUESTED_SCALE - 1.0), low_positions, True, failure_reason(invalid_validation)


def metrics(positions, edges):
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    width, height = max(xs) - min(xs), max(ys) - min(ys)
    return {"width": gr.rounded(width), "height": gr.rounded(height), "area": gr.rounded(width * height),
            "total_edge_length": gr.rounded(sum(math.dist(positions[a], positions[b]) for a, b in edges)),
            "centroid": abstract_graph.point_doc(centroid(positions))}


def viewport(starting):
    xs = [p[0] for p in starting.values()]
    ys = [p[1] for p in starting.values()]
    return min(xs) - 100, min(ys) - 100, max(xs) - min(xs) + 200, max(ys) - min(ys) + 200


def markup(positions, edges, center, selected, theta):
    selected = set(selected)
    lines = ['<g stroke="#37474f" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g>')
    for vertex in sorted(positions):
        color, radius = ("#f57c00", 5) if vertex in selected else ("#1565c0", 3)
        lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="{radius}" fill="{color}"><title>{vertex}</title></circle>')
    lines.extend([f'<circle cx="{center[0]:.6f}" cy="{center[1]:.6f}" r="7" fill="#d32f2f" stroke="white" stroke-width="2"/>',
                  f'<text x="{center[0]+10:.6f}" y="{center[1]-10:.6f}" font-family="monospace" font-size="14" fill="#d32f2f">C</text>'])
    if theta is not None:
        radians = math.radians(theta)
        end = (center[0] + 180 * math.cos(radians), center[1] + 180 * math.sin(radians))
        lines.extend(['<defs><marker id="fixed-arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#f57c00"/></marker></defs>',
                      f'<line x1="{center[0]:.6f}" y1="{center[1]:.6f}" x2="{end[0]:.6f}" y2="{end[1]:.6f}" stroke="#f57c00" stroke-width="3" marker-end="url(#fixed-arrow)"/>'])
    return lines


def write_static(path, positions, edges, center, viewbox, title):
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>{title}</title>'] + markup(positions, edges, center, (), None) + ['</svg>']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_motion(path, frames, edges, center, viewbox):
    count, epsilon = len(frames), 0.000001
    duration = count * 0.13
    left, top, width, height = viewbox
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Four cumulative fixed radial pulls</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        lines.append(f'@keyframes z{index} {{{rule}}}')
        lines.append(f'#z{index}{{opacity:0;animation:z{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="z{index}"><text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="21">{frame["label"]}</text>')
        lines.extend(markup(frame["positions"], edges, center, frame["selected"], frame["theta"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, original, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    current = load_vertices(output_dir / "directional_balance_start.svg", original)
    starting = {vertex: current[vertex] for vertex in sorted(current)}
    center = centroid(starting)
    if not valid(graph, embedding, starting, original, edges)["abstract_graph_valid"]:
        raise RuntimeError("starting abstract graph is invalid")
    frames = [{"label": "START", "positions": starting, "selected": (), "theta": None} for _ in range(PAUSE_FRAMES)]
    records = []
    previous_end = {vertex: current[vertex] for vertex in sorted(current)}
    for theta in DIRECTIONS:
        operation_start = {vertex: current[vertex] for vertex in sorted(current)}
        if operation_start != previous_end:
            raise RuntimeError(f"cumulative continuity failure before theta={theta}")
        selected = selected_vertices(operation_start, center, theta)
        accepted, current, limited, limiting = limited_pull(
            graph, embedding, original, edges, operation_start, center, selected)
        for step in range(1, MOTION_FRAMES + 1):
            scale = 1.0 + step / MOTION_FRAMES * (accepted - 1.0)
            frames.append({"label": f"theta={theta}°, requested 0.80, s={scale:.6f}",
                           "positions": radial_positions(operation_start, center, selected, scale),
                           "selected": selected, "theta": theta})
        previous_end = {vertex: current[vertex] for vertex in sorted(current)}
        if frames[-1]["positions"] != previous_end:
            raise RuntimeError(f"animation endpoint mismatch after theta={theta}")
        for _ in range(PAUSE_FRAMES):
            frames.append({"label": f"theta={theta}° accepted — pause", "positions": previous_end,
                           "selected": (), "theta": theta})
        after = metrics(current, edges)
        records.append({"theta": theta, "selected_nodes": len(selected), "requested_s": REQUESTED_SCALE,
                        "accepted_s": gr.rounded(accepted), "validity_limited": limited,
                        "limiting_event": limiting, "width_after": after["width"],
                        "height_after": after["height"], "area_after": after["area"],
                        "total_edge_length_after": after["total_edge_length"],
                        "cumulative_input_matches_previous_output": True})
    final_validation = valid(graph, embedding, current, original, edges)
    if not final_validation["abstract_graph_valid"]:
        raise RuntimeError("final graph is invalid")
    viewbox = viewport(starting)
    start_path = output_dir / "radial4_fixed_start.svg"
    end_path = output_dir / "radial4_fixed_end.svg"
    motion_path = output_dir / "radial4_fixed_motion.svg"
    report_path = output_dir / "radial4_fixed_report.json"
    write_static(start_path, starting, edges, center, viewbox, "Fixed radial4 start")
    write_static(end_path, current, edges, center, viewbox, "Fixed radial4 end")
    write_motion(motion_path, frames, edges, center, viewbox)
    start_metrics, final_metrics = metrics(starting, edges), metrics(current, edges)
    gr.write_json(report_path, {"model": "four cumulative fixed-strength radial half-plane pulls",
                                "physical_geometry_participated": False,
                                "frozen_centroid": abstract_graph.point_doc(center),
                                "directions": records, "start": start_metrics, "final": final_metrics,
                                "final_validation": final_validation,
                                "all_cumulative_continuity_assertions_passed": True,
                                "animation_frames": len(frames),
                                "fixed_viewport": {"x": viewbox[0], "y": viewbox[1],
                                                   "width": viewbox[2], "height": viewbox[3]}})
    print("theta | selected | requested s | accepted s | limited | width | height | area | edge length")
    for record in records:
        print(f'{record["theta"]:>5} | {record["selected_nodes"]:>8} | {record["requested_s"]:.6f} | '
              f'{record["accepted_s"]:.6f} | {str(record["validity_limited"]):>7} | '
              f'{record["width_after"]:.6f} | {record["height_after"]:.6f} | '
              f'{record["area_after"]:.6f} | {record["total_edge_length_after"]:.6f}')
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
