#!/usr/bin/env python3
"""One deterministic, genuinely separated two-pole contraction demonstration."""

from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_4u as two_pole
import graph_relax as gr


def select_pair(initial):
    candidates = []
    for pole_a, pole_b in combinations(sorted(initial), 2):
        model = two_pole.pair_model(initial, pole_a, pole_b)
        if model["theoretical_bound"] > two_pole.HEIGHT_BACKOFF:
            candidates.append(model)
    if not candidates:
        raise RuntimeError("no distinct existing-node pair has a positive two-circle height range")
    return sorted(candidates, key=lambda model: (-model["distance"], model["pole_a"], model["pole_b"]))[0]


def run(input_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = abstract_graph.load_abstract_start(circuit, output_dir)
    selected_model = select_pair(initial)
    baseline = two_pole.extents_and_distances(initial, edges)
    record, model, final_positions, validation = two_pole.search_pair(
        graph, embedding, initial, edges, baseline,
        selected_model["pole_a"], selected_model["pole_b"])
    if model["pole_a"] == model["pole_b"] or model["a"] == model["b"]:
        raise RuntimeError("two poles are not distinct")
    frames = []
    frame_records = []
    for fraction in two_pole.VISUAL_FRACTIONS:
        height = record["greatest_valid_equal_height"] if fraction == 1.0 else fraction * record["greatest_valid_equal_height"]
        valid, positions, frame_validation = two_pole.evaluate_height(
            graph, embedding, initial, edges, model, height)
        if not valid:
            raise RuntimeError(f"invalid demonstration frame at fraction {fraction}")
        label = "final" if fraction == 1.0 else f"{fraction:.2f}"
        frames.append((label, height, positions))
        measures = two_pole.extents_and_distances(positions, edges)
        frame_records.append({"fraction": fraction, "height": gr.rounded(height),
                              "width": gr.rounded(measures["width"]), "drawing_height": gr.rounded(measures["height"]),
                              "area": gr.rounded(measures["area"]), "edge_length": gr.rounded(measures["edge_length"]),
                              "crossings": frame_validation["unrelated_edge_crossings"],
                              "coincidences": frame_validation["coincident_vertex_pairs"], "valid": True})
    viewbox = two_pole.viewport(initial)
    animation_path = output_dir / "two_pole_distinct_contraction.svg"
    two_pole.write_animation(animation_path, frames, edges, model, viewbox)
    initial_path = output_dir / "two_pole_distinct_initial.svg"
    final_path = output_dir / "two_pole_distinct_final.svg"
    two_pole.write_svg(initial_path, frames[0][2], edges, model, viewbox, "Distinct two-pole initial graph")
    two_pole.write_svg(final_path, final_positions, edges, model, viewbox, "Distinct two-pole final graph")
    report_path = output_dir / "two_pole_distinct_report.json"
    gr.write_json(report_path, {
        "selection_rule": "maximum pole separation among distinct existing-node pairs with positive height range; stable IDs break ties",
        "pole_a": model["pole_a"], "pole_a_coordinate": abstract_graph.point_doc(model["a"]),
        "pole_b": model["pole_b"], "pole_b_coordinate": abstract_graph.point_doc(model["b"]),
        "pole_distance": gr.rounded(model["distance"]), "poles_are_distinct_vertices": True,
        "poles_have_distinct_coordinates": True, "best_height_for_this_pair": record["greatest_valid_equal_height"],
        "limiting_vertex": record["limiting_vertex"], "limiting_condition": record["limiting_condition"],
        "shared_viewport": {"x": viewbox[0], "y": viewbox[1], "width": viewbox[2], "height": viewbox[3]},
        "final": record, "frames": frame_records,
    })
    return {"pole_a": model["pole_a"], "pole_b": model["pole_b"],
            "pole_a_coordinate": abstract_graph.point_doc(model["a"]),
            "pole_b_coordinate": abstract_graph.point_doc(model["b"]),
            "pole_distance": gr.rounded(model["distance"]), "final": record,
            "animation": str(animation_path), "report": str(report_path)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
