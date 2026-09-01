#!/usr/bin/env python3
"""Experiment 4K: exhaustive independent integer X/Y scaling of 4I."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import experiment_4i as e4i
import experiment_4j as e4j
import graph_relax as gr


SCALE_VALUES = tuple(range(1, 53))
SCALE_PAIRS = tuple((sx, sy) for sx in SCALE_VALUES for sy in SCALE_VALUES)
SOURCE_TOLERANCE = 2e-5


def geometry_at_scales(circuit, unscaled_components, unscaled_nets, sx, sy):
    """Reconstruct one pair directly from immutable unscaled 4I points."""
    components = {ref: (point[0] * sx, point[1] * sy)
                  for ref, point in unscaled_components.items()}
    nets = {name: (point[0] * sx, point[1] * sy)
            for name, point in unscaled_nets.items()}
    attachments, edges = {}, []
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            attachment = e4i.ray_rectangle_attachment(components[ref], nets[net.name])
            attachments[terminal] = attachment
            edges.append({"terminal": terminal, "component": ref, "net": net.name,
                          "points": [attachment, nets[net.name]]})
    return components, nets, attachments, edges


def fully_valid(validation):
    return (validation["component_overlaps"] == 0 and
            validation["clearance_violations"] == 0 and
            validation["different_net_crossings"] == 0 and
            validation["unrelated_component_body_intersections"] == 0 and
            validation["electrical_incidence_exact"] and
            validation["expected_incidence_count"] == 44 and
            validation["realized_incidence_count"] == 44 and
            validation["cyclic_order_orientation"] in ("same", "global_reflection"))


def record_for_pair(circuit, graph, embedding, uc, un, sx, sy):
    geometry = geometry_at_scales(circuit, uc, un, sx, sy)
    validation = e4i.validate_geometry(circuit, graph, embedding, *geometry)
    metrics = e4i.measure(circuit, geometry[0], geometry[3])
    extents = e4j.geometry_extents(geometry[0], geometry[1])
    return ({"sx": sx, "sy": sy,
             "drawing_width": extents["width"], "drawing_height": extents["height"],
             "drawing_area": gr.rounded(extents["width"] * extents["height"]),
             "total_connection_length": gr.rounded(metrics["total_connection_length"]),
             "component_overlaps": validation["component_overlaps"],
             "clearance_violations": validation["clearance_violations"],
             "different_net_crossings": validation["different_net_crossings"],
             "unrelated_component_body_intersections": validation["unrelated_component_body_intersections"],
             "electrical_incidence_exact": validation["electrical_incidence_exact"],
             "expected_incidences": validation["expected_incidence_count"],
             "realized_incidences": validation["realized_incidence_count"],
             "added_incidences": 0 if validation["electrical_incidence_exact"] else None,
             "missing_incidences": 0 if validation["electrical_incidence_exact"] else None,
             "crowding": gr.rounded(metrics["energy"]["crowding"]),
             "energy": {key: gr.rounded(value) for key, value in metrics["energy"].items()},
             "objective": gr.rounded(metrics["objective"]),
             "cyclic_order_orientation": validation["cyclic_order_orientation"],
             "fully_valid": fully_valid(validation)}, geometry, validation, metrics, extents)


def sweep(circuit, graph, embedding, uc, un):
    records = []
    selected_payload = None
    for sx, sy in SCALE_PAIRS:
        record, geometry, validation, metrics, extents = record_for_pair(
            circuit, graph, embedding, uc, un, sx, sy)
        records.append(record)
        if record["fully_valid"]:
            key = (record["drawing_area"], record["total_connection_length"],
                   max(sx, sy), sx, sy)
            if selected_payload is None or key < selected_payload[0]:
                selected_payload = (key, record, geometry, validation, metrics, extents)
    if selected_payload is None:
        raise ValueError("no fully valid X/Y scale pair")
    return records, selected_payload


def failure_reasons(record):
    return e4j.failure_reasons(record)


def verify_baseline(circuit, graph, embedding, uc, un, output_dir):
    record, geometry, validation, metrics, extents = record_for_pair(
        circuit, graph, embedding, uc, un, 52, 52)
    if not record["fully_valid"]:
        raise ValueError("stored Experiment 4J (52,52) baseline is not fully valid")
    expected = json.loads((output_dir / "planar_scale_compact_metrics.json").read_text())
    expected_metrics = expected["experiment_4j"]
    if abs(metrics["total_connection_length"] - expected_metrics["total_connection_length"]) > SOURCE_TOLERANCE:
        raise ValueError("stored Experiment 4J connection-length baseline discrepancy")
    for key in ("component_overlaps", "clearance_violations", "net_crossings"):
        if metrics[key] != expected_metrics[key]:
            raise ValueError(f"stored Experiment 4J metric discrepancy: {key}")
    return record, geometry, validation, metrics, extents


def coordinate_document(circuit, sx, sy, geometry, validation):
    document = e4j.coordinate_document(circuit, {"sx": sx, "sy": sy}, geometry, validation)
    document.pop("scale")
    document["x_scale"] = sx
    document["y_scale"] = sy
    return document


def grid_text(records):
    by_pair = {(r["sx"], r["sy"]): r["fully_valid"] for r in records}
    lines = ["Experiment 4K fully-valid grid: # valid, . invalid",
             "Rows: SY=52 down to 1; columns: SX=1 through 52",
             "    " + "".join(str(sx // 10 or " ")[-1] for sx in SCALE_VALUES),
             "    " + "".join(str(sx % 10) for sx in SCALE_VALUES)]
    for sy in reversed(SCALE_VALUES):
        lines.append(f"{sy:2d}  " + "".join("#" if by_pair[sx, sy] else "." for sx in SCALE_VALUES))
    return "\n".join(lines) + "\n"


def percent_change(new, old):
    return gr.rounded(100.0 * (new - old) / old)


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    _, _, _, graph, embedding, uc, un = e4j.load_stored(circuit, output_dir)
    baseline_record, baseline_geometry, _, baseline_metrics, baseline_extents = verify_baseline(
        circuit, graph, embedding, uc, un, output_dir)
    records, selected = sweep(circuit, graph, embedding, uc, un)
    _, selected_record, geometry, validation, metrics, extents = selected
    sx, sy = selected_record["sx"], selected_record["sy"]
    valid_records = [record for record in records if record["fully_valid"]]
    record_map = {(record["sx"], record["sy"]): record for record in records}
    neighbors = {}
    for label, pair in (("lower_x", (sx - 1, sy)), ("lower_y", (sx, sy - 1))):
        record = record_map.get(pair)
        neighbors[label] = None if record is None else {
            "sx": pair[0], "sy": pair[1], "fully_valid": record["fully_valid"],
            "failure_reasons": [] if record["fully_valid"] else failure_reasons(record)}

    sweep_path = output_dir / "planar_xy_scale_sweep.json"
    sweep_document = {
        "scale_values": list(SCALE_VALUES), "pair_enumeration_order": "SX ascending, then SY ascending",
        "expected_pair_count": 2704, "tested_pair_count": len(records), "records": records,
        "fully_valid_pair_count": len(valid_records), "invalid_pair_count": len(records)-len(valid_records),
        "smallest_sx_in_any_valid_pair": min(r["sx"] for r in valid_records),
        "smallest_sy_in_any_valid_pair": min(r["sy"] for r in valid_records),
        "selection_rule": ["minimum drawing area", "minimum connection length", "minimum max(SX,SY)", "minimum SX", "minimum SY"],
        "selected_pair": {"sx": sx, "sy": sy}, "immediate_lower_neighbors": neighbors,
        "validity_grid": grid_text(records).splitlines()}
    gr.write_json(sweep_path, sweep_document)

    coordinate_path = output_dir / "planar_xy_compact_coordinates.json"
    gr.write_json(coordinate_path, coordinate_document(circuit, sx, sy, geometry, validation))
    original_coordinates = json.loads((output_dir / "planar_node_initial_coordinates.json").read_text())
    original_geometry = geometry_at_scales(circuit, uc, un, 72, 72)
    original_extents = e4j.geometry_extents(original_geometry[0], original_geometry[1])
    original_metrics = json.loads((output_dir / "planar_node_initial_metrics.json").read_text())["metrics"]
    comparisons = {}
    for name, old_extents, old_metrics in (("experiment_4i", original_extents, original_metrics),
                                           ("experiment_4j", baseline_extents, baseline_metrics)):
        comparisons[name] = {
            "width_percent_change": percent_change(extents["width"], old_extents["width"]),
            "height_percent_change": percent_change(extents["height"], old_extents["height"]),
            "area_percent_change": percent_change(extents["width"]*extents["height"], old_extents["width"]*old_extents["height"]),
            "connection_length_percent_change": percent_change(metrics["total_connection_length"], old_metrics["total_connection_length"])}
    metrics_document = {
        "selected_sx": sx, "selected_sy": sy, "geometry_extents": extents,
        "drawing_area": gr.rounded(extents["width"]*extents["height"]),
        "metrics": gr.rounded_metrics(metrics), "comparisons": comparisons,
        "experiment_4i_scale": {"sx": 72, "sy": 72},
        "experiment_4j_scale": {"sx": 52, "sy": 52}, "measurement_only": True}
    metrics_path = output_dir / "planar_xy_compact_metrics.json"
    gr.write_json(metrics_path, metrics_document)
    validation_document = {
        "baseline_4j_52_52_reproduced": True, "embedding_unchanged": True,
        "unscaled_coordinates_unchanged": True, "incidence_graph_unchanged": True,
        "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
        "attachment_algorithm": "Experiment 4I ray/rectangle boundary intersection",
        "selected_pair": {"sx": sx, "sy": sy}, "fully_valid": True,
        "validation": validation, "immediate_lower_neighbors": neighbors}
    validation_path = output_dir / "planar_xy_compact_validation.json"
    gr.write_json(validation_path, validation_document)
    svg_path = output_dir / "planar_xy_compact.svg"
    e4i.write_svg(svg_path, circuit, *geometry)
    grid_path = output_dir / "planar_xy_scale_grid.txt"
    grid_path.write_text(grid_text(records), encoding="utf-8")
    paths = (sweep_path, coordinate_path, metrics_path, validation_path, svg_path, grid_path)
    return {"selected_pair": {"sx": sx, "sy": sy}, "selected_record": selected_record,
            "valid_pair_count": len(valid_records), "invalid_pair_count": len(records)-len(valid_records),
            "smallest_sx": sweep_document["smallest_sx_in_any_valid_pair"],
            "smallest_sy": sweep_document["smallest_sy_in_any_valid_pair"],
            "neighbors": neighbors, "metrics": metrics_document,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
