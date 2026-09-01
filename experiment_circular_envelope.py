#!/usr/bin/env python3
"""Measure the circular support-envelope residual of the unchanged 4K points."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import graph_relax as gr


SAMPLE_COUNT = 360
SUPPORT_TOLERANCE = 1e-9
RESIDUAL_CHECK_TOLERANCE = 1e-7


def arithmetic_centroid(points):
    values = list(points.values())
    return (sum(p[0] for p in values) / len(values), sum(p[1] for p in values) / len(values))


def support_samples(points):
    samples = []
    for theta in range(SAMPLE_COUNT):
        radians = math.radians(theta)
        cosine, sine = math.cos(radians), math.sin(radians)
        projections = {vertex: point[0] * cosine + point[1] * sine for vertex, point in points.items()}
        support = max(projections.values())
        tolerance = SUPPORT_TOLERANCE * max(1.0, abs(support))
        supporters = sorted(vertex for vertex, value in projections.items() if abs(value - support) <= tolerance)
        samples.append({"theta": theta, "cos": cosine, "sin": sine, "H": support,
                        "support_vertices": supporters})
    return samples


def fit_circle(samples):
    # For a complete equally spaced circle, [cos(theta), sin(theta), 1]
    # are mutually orthogonal, so these are the ordinary least-squares coefficients.
    count = len(samples)
    cx = 2.0 / count * sum(sample["H"] * sample["cos"] for sample in samples)
    cy = 2.0 / count * sum(sample["H"] * sample["sin"] for sample in samples)
    radius = sum(sample["H"] for sample in samples) / count
    for sample in samples:
        fitted = cx * sample["cos"] + cy * sample["sin"] + radius
        sample["H_circle"] = fitted
        sample["e"] = sample["H"] - fitted
    return cx, cy, radius


def interpolated_zero(theta_a, error_a, error_b):
    return theta_a + (-error_a) / (error_b - error_a)


def positive_lobes(samples):
    positive = [sample["e"] > 0.0 for sample in samples]
    starts = [index for index in range(SAMPLE_COUNT)
              if positive[index] and not positive[(index - 1) % SAMPLE_COUNT]]
    lobes = []
    for start_index in starts:
        indices = []
        index = start_index
        while positive[index]:
            indices.append(index)
            index = (index + 1) % SAMPLE_COUNT
            if index == start_index:
                raise RuntimeError("all residual samples are positive, inconsistent with least squares")
        previous = (start_index - 1) % SAMPLE_COUNT
        last = indices[-1]
        following = (last + 1) % SAMPLE_COUNT
        start_base = start_index - 1
        if start_index == 0:
            start_base = -1
        start_angle = interpolated_zero(start_base, samples[previous]["e"], samples[start_index]["e"])
        last_unwrapped = last
        if last < start_index:
            last_unwrapped += 360
        end_angle = interpolated_zero(last_unwrapped, samples[last]["e"], samples[following]["e"])
        unwrapped_samples = []
        for sample_index in indices:
            angle = sample_index + (360 if sample_index < start_index else 0)
            unwrapped_samples.append((angle, samples[sample_index]["e"], sample_index))
        integration_points = [(start_angle, 0.0)] + [(a, e) for a, e, _ in unwrapped_samples] + [(end_angle, 0.0)]
        integrated = sum((integration_points[i + 1][0] - integration_points[i][0]) *
                         (integration_points[i][1] + integration_points[i + 1][1]) / 2
                         for i in range(len(integration_points) - 1))
        vector_x = sum(error * math.cos(math.radians(index)) for _, error, index in unwrapped_samples)
        vector_y = sum(error * math.sin(math.radians(index)) for _, error, index in unwrapped_samples)
        center = math.degrees(math.atan2(vector_y, vector_x)) % 360
        peak_index = max(indices, key=lambda i: (samples[i]["e"], -i))
        supporters = sorted({vertex for i in indices for vertex in samples[i]["support_vertices"]})
        lobes.append({"start": start_angle % 360, "end": end_angle % 360,
                      "start_unwrapped": start_angle, "end_unwrapped": end_angle,
                      "width": end_angle - start_angle, "center": center,
                      "peak_angle": peak_index, "peak_error": samples[peak_index]["e"],
                      "integrated_positive_error": integrated,
                      "support_vertices": supporters, "sample_indices": indices})
    return sorted(lobes, key=lambda lobe: (-lobe["integrated_positive_error"], lobe["center"]))


def convex_hull(points):
    ordered = sorted((point[0], point[1], vertex) for vertex, point in points.items())
    def cross(origin, a, b):
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])
    lower = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def write_measurement_svg(path, points, edges, center, radius, hull, largest_lobe):
    circle_extremes = [(center[0] - radius, center[1] - radius), (center[0] + radius, center[1] + radius)]
    all_points = list(points.values()) + circle_extremes
    margin = 100
    left, right = min(p[0] for p in all_points) - margin, max(p[0] for p in all_points) + margin
    top, bottom = min(p[1] for p in all_points) - margin, max(p[1] for p in all_points) + margin
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {right-left:.6f} {bottom-top:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{right-left:.6f}" height="{bottom-top:.6f}" fill="white"/>',
             '<title>Circular support-envelope measurement</title>',
             f'<circle cx="{center[0]:.6f}" cy="{center[1]:.6f}" r="{radius:.6f}" fill="none" stroke="#2ca02c" stroke-width="3" stroke-dasharray="10 7"/>',
             '<g stroke="#b0bec5" stroke-width="1.2" fill="none">']
    for first, second in edges:
        lines.append(f'<line x1="{points[first][0]:.6f}" y1="{points[first][1]:.6f}" x2="{points[second][0]:.6f}" y2="{points[second][1]:.6f}"/>')
    lines.append('</g>')
    hull_points = " ".join(f'{x:.6f},{y:.6f}' for x, y, _ in hull)
    lines.append(f'<polygon points="{hull_points}" fill="none" stroke="#6a1b9a" stroke-width="3"/>')
    ray_length = radius * 1.12
    for angle, color, width in ((largest_lobe["start"], "#ef6c00", 2),
                                (largest_lobe["end"], "#ef6c00", 2),
                                (largest_lobe["center"], "#d32f2f", 4)):
        radians = math.radians(angle)
        endpoint = (center[0] + ray_length * math.cos(radians), center[1] + ray_length * math.sin(radians))
        lines.append(f'<line x1="{center[0]:.6f}" y1="{center[1]:.6f}" x2="{endpoint[0]:.6f}" y2="{endpoint[1]:.6f}" stroke="{color}" stroke-width="{width}"/>')
    highlighted = set(largest_lobe["support_vertices"])
    for vertex in sorted(points):
        color, size = ("#e53935", 8) if vertex in highlighted else ("#1565c0", 4)
        lines.append(f'<circle cx="{points[vertex][0]:.6f}" cy="{points[vertex][1]:.6f}" r="{size}" fill="{color}"><title>{vertex}</title></circle>')
    lines.extend([f'<circle cx="{center[0]:.6f}" cy="{center[1]:.6f}" r="8" fill="#2ca02c" stroke="white" stroke-width="2"/>',
                  f'<text x="{center[0]+12:.6f}" y="{center[1]-12:.6f}" font-family="monospace" font-size="16" fill="#2ca02c">fitted C</text>', '</svg>'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_error_svg(path, samples, lobes):
    width, height = 1100, 520
    left, right, top, bottom = 75, 35, 45, 65
    plot_width, plot_height = width - left - right, height - top - bottom
    maximum = max(abs(sample["e"]) for sample in samples) * 1.10
    def sx(theta): return left + theta / 360 * plot_width
    def sy(error): return top + (maximum - error) / (2 * maximum) * plot_height
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
             f'<rect width="{width}" height="{height}" fill="white"/>', '<title>Circular envelope residual by support direction</title>']
    for lobe in lobes:
        intervals = [(lobe["start_unwrapped"], lobe["end_unwrapped"])]
        for start, end in intervals:
            if start < 0:
                lines.append(f'<rect x="{sx(0):.3f}" y="{top}" width="{sx(end)-sx(0):.3f}" height="{plot_height}" fill="#ffccbc" opacity="0.55"/>')
                lines.append(f'<rect x="{sx(start+360):.3f}" y="{top}" width="{sx(360)-sx(start+360):.3f}" height="{plot_height}" fill="#ffccbc" opacity="0.55"/>')
            elif end > 360:
                lines.append(f'<rect x="{sx(start):.3f}" y="{top}" width="{sx(360)-sx(start):.3f}" height="{plot_height}" fill="#ffccbc" opacity="0.55"/>')
                lines.append(f'<rect x="{sx(0):.3f}" y="{top}" width="{sx(end-360)-sx(0):.3f}" height="{plot_height}" fill="#ffccbc" opacity="0.55"/>')
            else:
                lines.append(f'<rect x="{sx(start):.3f}" y="{top}" width="{sx(end)-sx(start):.3f}" height="{plot_height}" fill="#ffccbc" opacity="0.55"/>')
    zero_y = sy(0)
    lines.append(f'<line x1="{left}" y1="{zero_y:.3f}" x2="{left+plot_width}" y2="{zero_y:.3f}" stroke="#555" stroke-width="1.5"/>')
    for angle in (0, 90, 180, 270, 360):
        x = sx(angle)
        lines.append(f'<line x1="{x:.3f}" y1="{top}" x2="{x:.3f}" y2="{top+plot_height}" stroke="#ddd"/>')
        lines.append(f'<text x="{x:.3f}" y="{height-28}" text-anchor="middle" font-family="sans-serif" font-size="14">{angle}°</text>')
    plot_samples = samples + [{**samples[0], "theta": 360}]
    points = " ".join(f'{sx(sample["theta"]):.3f},{sy(sample["e"]):.3f}' for sample in plot_samples)
    lines.append(f'<polyline points="{points}" fill="none" stroke="#1565c0" stroke-width="2"/>')
    for rank, lobe in enumerate(lobes, 1):
        x = sx(lobe["center"])
        lines.append(f'<text x="{x:.3f}" y="{top+18+(rank%3)*17}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#b71c1c">L{rank}: {lobe["center"]:.1f}°</text>')
    lines.extend([f'<text x="{width/2}" y="{height-5}" text-anchor="middle" font-family="sans-serif" font-size="15">support direction θ</text>',
                  f'<text transform="translate(18 {height/2}) rotate(-90)" text-anchor="middle" font-family="sans-serif" font-size="15">e(θ)</text>', '</svg>'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    _, _, points, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    original = {vertex: points[vertex] for vertex in sorted(points)}
    samples = support_samples(points)
    cx, cy, radius = fit_circle(samples)
    residuals = [sample["e"] for sample in samples]
    rms = math.sqrt(sum(error * error for error in residuals) / SAMPLE_COUNT)
    checks = {"mean_e": sum(residuals) / SAMPLE_COUNT,
              "sum_e_cos": sum(sample["e"] * sample["cos"] for sample in samples),
              "sum_e_sin": sum(sample["e"] * sample["sin"] for sample in samples)}
    if max(abs(value) for value in checks.values()) > RESIDUAL_CHECK_TOLERANCE:
        raise RuntimeError(f"least-squares normal-equation residual checks failed: {checks}")
    lobes = positive_lobes(samples)
    if not lobes:
        raise RuntimeError("no positive residual lobes found")
    hull = convex_hull(points)
    xs, ys = [p[0] for p in points.values()], [p[1] for p in points.values()]
    arithmetic = arithmetic_centroid(points)
    bbox_center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    maximum_index = max(range(SAMPLE_COUNT), key=lambda i: samples[i]["e"])
    minimum_index = min(range(SAMPLE_COUNT), key=lambda i: samples[i]["e"])
    report_samples = [{"theta": sample["theta"], "H": gr.rounded(sample["H"]),
                       "H_circle": gr.rounded(sample["H_circle"]), "e": gr.rounded(sample["e"]),
                       "support_vertices": sample["support_vertices"]} for sample in samples]
    report_lobes = [{key: (gr.rounded(value) if isinstance(value, float) else value)
                     for key, value in lobe.items() if key not in ("sample_indices", "start_unwrapped", "end_unwrapped")}
                    for lobe in lobes]
    report = {"source": "authoritative output/planar_xy_compact_coordinates.json coordinates represented by radial4_max_start.svg",
              "coordinates_unchanged": original == points, "vertex_count": len(points), "sample_count": SAMPLE_COUNT,
              "fitted_circle": {"Cx": gr.rounded(cx), "Cy": gr.rounded(cy), "R": gr.rounded(radius)},
              "arithmetic_centroid": abstract_graph.point_doc(arithmetic),
              "bounding_box_center": abstract_graph.point_doc(bbox_center),
              "RMS_error": gr.rounded(rms), "normalized_RMS_error": gr.rounded(rms / radius),
              "maximum_positive_residual": {"error": gr.rounded(samples[maximum_index]["e"]), "theta": maximum_index},
              "maximum_negative_residual": {"error": gr.rounded(samples[minimum_index]["e"]), "theta": minimum_index},
              "peak_to_peak_residual": gr.rounded(samples[maximum_index]["e"] - samples[minimum_index]["e"]),
              "least_squares_residual_checks": {key: gr.rounded(value) for key, value in checks.items()},
              "samples": report_samples, "positive_lobes": report_lobes}
    gr.write_json(output_dir / "circular_envelope_report.json", report)
    write_measurement_svg(output_dir / "circular_envelope_measurement.svg", points, edges,
                          (cx, cy), radius, hull, lobes[0])
    write_error_svg(output_dir / "circular_envelope_error.svg", samples, lobes)
    print(f"fitted center: ({cx:.6f}, {cy:.6f})")
    print(f"fitted radius: {radius:.6f}")
    print(f"normalized RMS envelope error: {rms/radius:.9f}")
    print(f"positive lobes: {len(lobes)}")
    for rank, lobe in enumerate(lobes, 1):
        print(f'{rank}: center={lobe["center"]:.6f}°, width={lobe["width"]:.6f}°, '
              f'start={lobe["start"]:.6f}°, end={lobe["end"]:.6f}°, peak={lobe["peak_error"]:.6f}, '
              f'integrated={lobe["integrated_positive_error"]:.6f}, supporters={",".join(lobe["support_vertices"])}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.input, args.output_dir)


if __name__ == "__main__":
    main()
