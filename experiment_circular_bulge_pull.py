#!/usr/bin/env python3
"""Part 1: cumulatively pull measured positive support-envelope bulges inward."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_circular_envelope as envelope
import experiment_radial4_fixed as radial
import graph_relax as gr


EMERGENCY_ITERATION_CAP = 500
MEANINGFUL_SCALE_CHANGE = 1e-10
FALLBACK_SCAN_STEPS = 256
FALLBACK_BISECTION_TOLERANCE = 1e-9
MOTION_FRAMES = 12
PAUSE_FRAMES = 5


def measure(points):
    samples = envelope.support_samples(points)
    cx, cy, radius = envelope.fit_circle(samples)
    errors = [sample["e"] for sample in samples]
    rms = math.sqrt(sum(error * error for error in errors) / len(errors))
    lobes = envelope.positive_lobes(samples)
    tolerance = 1e-9 * max(1.0, abs(radius))
    return {"samples": samples, "center": (cx, cy), "radius": radius,
            "rms": rms, "normalized_rms": rms / radius, "lobes": lobes,
            "positive_tolerance": tolerance, "maximum_positive": max(errors)}


def selected_for_lobe(points, measurement, lobe):
    cx, cy = measurement["center"]
    radius = measurement["radius"]
    tolerance = 1e-9 * max(1.0, abs(radius))
    selected = []
    contributions = {}
    for vertex in sorted(points):
        dx, dy = points[vertex][0] - cx, points[vertex][1] - cy
        values = [dx * measurement["samples"][index]["cos"] +
                  dy * measurement["samples"][index]["sin"] for index in lobe["sample_indices"]]
        contributions[vertex] = max(values)
        if any(value > radius + tolerance for value in values):
            selected.append(vertex)
    if not selected:
        # The strict defect exists; tolerances may only collapse at machine scale.
        return [], None
    maximum = max(contributions[vertex] for vertex in selected)
    return selected, maximum


def moved(points, center, selected, scale):
    result = {vertex: points[vertex] for vertex in sorted(points)}
    for vertex in selected:
        result[vertex] = (center[0] + scale * (points[vertex][0] - center[0]),
                          center[1] + scale * (points[vertex][1] - center[1]))
    return result


def validity(graph, embedding, points, original, edges):
    return abstract_graph.validate_abstract(graph, embedding, points, original, edges)


def validity_event(validation):
    reasons, details = [], {}
    if validation["unrelated_edge_crossings"]:
        reasons.append("unrelated edge crossing"); details["crossing_edges"] = validation["crossing_details"]
    if validation["coincident_vertex_pairs"]:
        reasons.append("coincident vertices"); details["coincident_vertices"] = validation["coincident_details"]
    if validation["vertices_on_unrelated_edges"]:
        reasons.append("vertex on unrelated edge"); details["vertex_on_edge"] = validation["vertex_on_edge_details"]
    if validation["cyclic_order_violations"]:
        reasons.append("cyclic-order/embedding violation"); details["cyclic_vertices"] = validation["cyclic_violation_vertices"]
    return ", ".join(reasons) or "unknown abstract validity event", details


def fallback_pull(graph, embedding, original, edges, points, center, selected, target):
    last_t, last_scale, last_points = 0.0, 1.0, points
    invalid_t = None
    invalid_validation = None
    for step in range(1, FALLBACK_SCAN_STEPS + 1):
        t = step / FALLBACK_SCAN_STEPS
        scale = 1.0 + t * (target - 1.0)
        candidate = moved(points, center, selected, scale)
        check = validity(graph, embedding, candidate, original, edges)
        if not check["abstract_graph_valid"]:
            invalid_t, invalid_validation = t, check
            break
        last_t, last_scale, last_points = t, scale, candidate
    if invalid_t is None:
        raise RuntimeError("direct target was invalid but fallback scan did not reproduce an invalid state")
    low, high = last_t, invalid_t
    iterations = 0
    while high - low > FALLBACK_BISECTION_TOLERANCE and iterations < 48:
        middle = (low + high) / 2
        scale = 1.0 + middle * (target - 1.0)
        candidate = moved(points, center, selected, scale)
        check = validity(graph, embedding, candidate, original, edges)
        iterations += 1
        if check["abstract_graph_valid"]:
            low, last_scale, last_points = middle, scale, candidate
        else:
            high, invalid_validation = middle, check
    reason, details = validity_event(invalid_validation)
    return last_scale, last_points, reason, details


def geometry_metrics(points, edges):
    xs, ys = [p[0] for p in points.values()], [p[1] for p in points.values()]
    width, height = max(xs) - min(xs), max(ys) - min(ys)
    return {"width": gr.rounded(width), "height": gr.rounded(height),
            "area": gr.rounded(width * height),
            "total_edge_length": gr.rounded(sum(math.dist(points[a], points[b]) for a, b in edges))}


def viewbox(starting, initial_measurement):
    cx, cy = initial_measurement["center"]
    radius = initial_measurement["radius"]
    xs = [p[0] for p in starting.values()] + [cx - radius, cx + radius]
    ys = [p[1] for p in starting.values()] + [cy - radius, cy + radius]
    margin = 100
    return min(xs) - margin, min(ys) - margin, max(xs) - min(xs) + 2 * margin, max(ys) - min(ys) + 2 * margin


def markup(points, edges, measurement, selected, lobe):
    cx, cy = measurement["center"]
    radius = measurement["radius"]
    selected = set(selected)
    lines = [f'<circle cx="{cx:.6f}" cy="{cy:.6f}" r="{radius:.6f}" fill="none" stroke="#2e7d32" stroke-width="3" stroke-dasharray="10 7"/>',
             '<g stroke="#90a4ae" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{points[first][0]:.6f}" y1="{points[first][1]:.6f}" x2="{points[second][0]:.6f}" y2="{points[second][1]:.6f}"/>')
    lines.append('</g>')
    if lobe is not None:
        ray = radius * 1.15
        for angle, color, width in ((lobe["start"], "#ef6c00", 2), (lobe["end"], "#ef6c00", 2),
                                    (lobe["center"], "#c62828", 4)):
            radians = math.radians(angle)
            endpoint = (cx + ray * math.cos(radians), cy + ray * math.sin(radians))
            lines.append(f'<line x1="{cx:.6f}" y1="{cy:.6f}" x2="{endpoint[0]:.6f}" y2="{endpoint[1]:.6f}" stroke="{color}" stroke-width="{width}"/>')
    for vertex in sorted(points):
        color, size = ("#f57c00", 7) if vertex in selected else ("#1565c0", 4)
        lines.append(f'<circle cx="{points[vertex][0]:.6f}" cy="{points[vertex][1]:.6f}" r="{size}" fill="{color}"><title>{vertex}</title></circle>')
    lines.extend([f'<circle cx="{cx:.6f}" cy="{cy:.6f}" r="8" fill="#2e7d32" stroke="white" stroke-width="2"/>',
                  f'<text x="{cx+12:.6f}" y="{cy-12:.6f}" font-family="monospace" font-size="15" fill="#2e7d32">fitted C</text>'])
    return lines


def write_static(path, points, edges, measurement, view, title):
    left, top, width, height = view
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             f'<title>{title}</title>'] + markup(points, edges, measurement, (), None) + ['</svg>']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_motion(path, frames, edges, view):
    count, epsilon = len(frames), 0.000001
    duration = count * 0.14
    left, top, width, height = view
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {width:.6f} {height:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{width:.6f}" height="{height:.6f}" fill="white"/>',
             '<title>Circular positive-bulge pulls</title>', '<style>']
    for index in range(count):
        start, end = 100 * index / count, 100 * (index + 1) / count
        if index == 0:
            rule = f'0%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        elif index == count - 1:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,100%{{opacity:1}}'
        else:
            rule = f'0%,{start:.7f}%{{opacity:0}} {start+epsilon:.7f}%,{end:.7f}%{{opacity:1}} {end+epsilon:.7f}%,100%{{opacity:0}}'
        lines.append(f'@keyframes b{index} {{{rule}}}')
        lines.append(f'#b{index}{{opacity:0;animation:b{index} {duration:.2f}s steps(1,end) infinite}}')
    lines.append('</style>')
    for index, frame in enumerate(frames):
        lines.append(f'<g id="b{index}"><text x="{left+25:.6f}" y="{top+35:.6f}" font-family="monospace" font-size="20">{frame["label"]}</text>')
        lines.extend(markup(frame["points"], edges, frame["measurement"], frame["selected"], frame["lobe"])); lines.append('</g>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, authoritative, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    starting = {vertex: authoritative[vertex] for vertex in sorted(authoritative)}
    current = {vertex: starting[vertex] for vertex in sorted(starting)}
    initial_measurement = measure(current)
    initial_metrics = geometry_metrics(current, edges)
    view = viewbox(starting, initial_measurement)
    frames = [{"label": "START — measured fitted circle", "points": current,
               "measurement": initial_measurement, "selected": (), "lobe": None} for _ in range(PAUSE_FRAMES)]
    iteration_records = []
    terminal_blocked_lobes = []
    total_pulls = 0
    termination = None
    for iteration in range(1, EMERGENCY_ITERATION_CAP + 1):
        before_measurement = measure(current)
        if before_measurement["maximum_positive"] <= before_measurement["positive_tolerance"]:
            termination = "no positive lobe above numerical fitting tolerance"
            break
        accepted_this_measurement = False
        attempt_records = []
        for lobe_rank, lobe in enumerate(before_measurement["lobes"], 1):
            selected, maximum = selected_for_lobe(current, before_measurement, lobe)
            if not selected or maximum is None:
                attempt_records.append({"lobe_rank": lobe_rank, "blocked": True,
                                        "reason": "no vertex exceeds fitted support beyond deterministic tolerance"})
                continue
            radius = before_measurement["radius"]
            target = radius / maximum
            if not 0 < target < 1:
                attempt_records.append({"lobe_rank": lobe_rank, "blocked": True,
                                        "reason": "derived target has no meaningful inward contraction"})
                continue
            candidate = moved(current, before_measurement["center"], selected, target)
            direct_check = validity(graph, embedding, candidate, starting, edges)
            direct = direct_check["abstract_graph_valid"]
            if direct:
                accepted, accepted_points, limiting, limiting_details = target, candidate, None, {}
            else:
                accepted, accepted_points, limiting, limiting_details = fallback_pull(
                    graph, embedding, starting, edges, current, before_measurement["center"], selected, target)
            if 1.0 - accepted <= MEANINGFUL_SCALE_CHANGE:
                attempt_records.append({"lobe_rank": lobe_rank, "blocked": True,
                                        "reason": limiting or "no meaningful valid inward correction",
                                        "limiting_details": limiting_details})
                continue
            operation_start = {vertex: current[vertex] for vertex in sorted(current)}
            for step in range(1, MOTION_FRAMES + 1):
                scale = 1.0 + step / MOTION_FRAMES * (accepted - 1.0)
                frames.append({"label": f"iteration {iteration}, lobe {lobe_rank}, s={scale:.9f}",
                               "points": moved(operation_start, before_measurement["center"], selected, scale),
                               "measurement": before_measurement, "selected": selected, "lobe": lobe})
            current = accepted_points
            total_pulls += 1
            after_measurement = measure(current)
            for _ in range(PAUSE_FRAMES):
                frames.append({"label": f"iteration {iteration} accepted — circle refitted",
                               "points": current, "measurement": after_measurement,
                               "selected": (), "lobe": None})
            iteration_records.append({"iteration": iteration,
                "fitted_before": {"Cx": gr.rounded(before_measurement["center"][0]),
                                  "Cy": gr.rounded(before_measurement["center"][1]),
                                  "R": gr.rounded(before_measurement["radius"])},
                "normalized_RMS_before": gr.rounded(before_measurement["normalized_rms"]),
                "positive_lobe_count": len(before_measurement["lobes"]),
                "chosen_lobe": {"start": gr.rounded(lobe["start"]), "end": gr.rounded(lobe["end"]),
                                "width": gr.rounded(lobe["width"]), "center": gr.rounded(lobe["center"]),
                                "integrated_positive_error": gr.rounded(lobe["integrated_positive_error"])},
                "selected_vertices": selected, "M": gr.rounded(maximum), "s_target": gr.rounded(target),
                "accepted_s": gr.rounded(accepted), "direct_target_reached": direct,
                "limiting_graph_validity_event": limiting,
                "limiting_details": limiting_details, "earlier_blocked_lobe_attempts": attempt_records,
                "fitted_after": {"Cx": gr.rounded(after_measurement["center"][0]),
                                 "Cy": gr.rounded(after_measurement["center"][1]),
                                 "R": gr.rounded(after_measurement["radius"])},
                "normalized_RMS_after": gr.rounded(after_measurement["normalized_rms"])})
            accepted_this_measurement = True
            break
        if not accepted_this_measurement:
            terminal_blocked_lobes = attempt_records
            termination = "every remaining positive lobe blocked from meaningful inward correction"
            break
    if termination is None:
        termination = "emergency iteration cap reached"
    final_measurement = measure(current)
    final_check = validity(graph, embedding, current, starting, edges)
    if not final_check["abstract_graph_valid"]:
        raise RuntimeError("final graph is invalid")
    final_metrics = geometry_metrics(current, edges)
    start_path = output_dir / "circular_bulge_pull_start.svg"
    end_path = output_dir / "circular_bulge_pull_end.svg"
    motion_path = output_dir / "circular_bulge_pull_motion.svg"
    report_path = output_dir / "circular_bulge_pull_report.json"
    write_static(start_path, starting, edges, initial_measurement, view, "Circular bulge pull start")
    write_static(end_path, current, edges, final_measurement, view, "Circular bulge pull end")
    write_motion(motion_path, frames, edges, view)
    report = {"part": 1, "source": "authoritative 4K coordinates used by circular_envelope_measurement.svg",
              "outward_or_inverse_push_implemented": False, "emergency_iteration_cap": EMERGENCY_ITERATION_CAP,
              "termination": termination, "total_pulls": total_pulls,
              "initial_normalized_RMS": gr.rounded(initial_measurement["normalized_rms"]),
              "final_normalized_RMS": gr.rounded(final_measurement["normalized_rms"]),
              "initial_metrics": initial_metrics, "final_metrics": final_metrics,
              "terminal_blocked_lobes": terminal_blocked_lobes,
              "final_graph_validity": final_check, "iterations": iteration_records,
              "fixed_viewport": {"x": view[0], "y": view[1], "width": view[2], "height": view[3]}}
    gr.write_json(report_path, report)
    print(json.dumps({"termination": termination, "total_pulls": total_pulls,
                      "initial_normalized_RMS": report["initial_normalized_RMS"],
                      "final_normalized_RMS": report["final_normalized_RMS"],
                      "initial_metrics": initial_metrics, "final_metrics": final_metrics,
                      "outputs": [str(start_path), str(end_path), str(motion_path), str(report_path)]},
                     sort_keys=True, indent=2))
    for record in iteration_records:
        lobe = record["chosen_lobe"]
        print(f'iteration {record["iteration"]}: C=({record["fitted_before"]["Cx"]:.6f},'
              f'{record["fitted_before"]["Cy"]:.6f}), R={record["fitted_before"]["R"]:.6f}, '
              f'NRMS={record["normalized_RMS_before"]:.6f}, lobes={record["positive_lobe_count"]}, '
              f'lobe={lobe["start"]:.6f}°->{lobe["end"]:.6f}° center={lobe["center"]:.6f}° '
              f'integrated={lobe["integrated_positive_error"]:.6f}, selected={record["selected_vertices"]}, '
              f'M={record["M"]:.6f}, target={record["s_target"]:.9f}, accepted={record["accepted_s"]:.9f}, '
              f'direct={record["direct_target_reached"]}, after_NRMS={record["normalized_RMS_after"]:.6f}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.input, args.output_dir)


if __name__ == "__main__":
    main()
