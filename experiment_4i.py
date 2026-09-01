#!/usr/bin/env python3
"""Experiment 4I: deterministic planar incidence-node initialization."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import networkx as nx

import experiment_4p as e4p
import graph_relax as gr


INITIAL_SCALE = max(gr.WIDTH, gr.HEIGHT) + gr.BODY_CLEARANCE  # 72
SCALE_SEQUENCE = tuple(INITIAL_SCALE * (2 ** exponent) for exponent in range(20))
EPSILON = 1e-7


def stable_embedding(circuit):
    graph = e4p.incidence_graph(circuit)
    if (graph.number_of_nodes(), graph.number_of_edges()) != (34, 44):
        raise ValueError("4I incidence graph must have 34 vertices and 44 edges")
    planar, source_embedding = nx.check_planarity(graph, counterexample=True)
    if not planar:
        raise ValueError("4I incidence graph unexpectedly nonplanar")
    names = sorted(graph)
    name_to_integer = {name: index for index, name in enumerate(names)}
    integer_to_name = {index: name for name, index in name_to_integer.items()}
    rotations = {name_to_integer[name]: [name_to_integer[neighbor]
                                         for neighbor in source_embedding.neighbors_cw_order(name)]
                 for name in names}
    embedding = nx.PlanarEmbedding()
    embedding.set_data(rotations)
    embedding.check_structure()
    positions = nx.combinatorial_embedding_to_pos(embedding, fully_triangulate=False)
    return graph, source_embedding, embedding, name_to_integer, integer_to_name, positions


def ray_rectangle_attachment(center, target):
    dx, dy = target[0] - center[0], target[1] - center[1]
    if dx == 0 and dy == 0:
        raise ValueError("component and incident net point coincide")
    x_scale = (gr.WIDTH / 2.0) / abs(dx) if dx else math.inf
    y_scale = (gr.HEIGHT / 2.0) / abs(dy) if dy else math.inf
    scale = min(x_scale, y_scale)
    return center[0] + scale * dx, center[1] + scale * dy


def realize(circuit, integer_to_name, point_positions, scale):
    points = {integer_to_name[index]: (float(point[0]) * scale, float(point[1]) * scale)
              for index, point in point_positions.items()}
    components = {ref: points[e4p.component_vertex(ref)] for ref in circuit.refs}
    nets = {net.name: points[e4p.net_vertex(net.name)] for net in circuit.nets}
    attachments, edges = {}, []
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            attachment = ray_rectangle_attachment(components[ref], nets[net.name])
            attachments[terminal] = attachment
            edges.append({"terminal": terminal, "component": ref, "net": net.name,
                          "points": [attachment, nets[net.name]]})
    return components, nets, attachments, edges


def cyclic_equal(first, second):
    if len(first) != len(second):
        return False
    if not first:
        return True
    return any(first == second[offset:] + second[:offset] for offset in range(len(second)))


def perimeter_parameter(point, center):
    x, y = point[0] - center[0], point[1] - center[1]
    half_width, half_height = gr.WIDTH / 2.0, gr.HEIGHT / 2.0
    if abs(y + half_height) <= EPSILON:
        return x + half_width
    if abs(x - half_width) <= EPSILON:
        return gr.WIDTH + y + half_height
    if abs(y - half_height) <= EPSILON:
        return gr.WIDTH + gr.HEIGHT + half_width - x
    if abs(x + half_width) <= EPSILON:
        return 2.0 * gr.WIDTH + gr.HEIGHT + half_height - y
    raise ValueError("attachment is not on component perimeter")


def geometric_component_orders(circuit, components, attachments, embedding):
    terminal_for = {(terminal.rsplit(".", 1)[0], net.name): terminal
                    for net in circuit.nets for terminal in net.terminals}
    result = {}
    for ref in circuit.refs:
        vertex = e4p.component_vertex(ref)
        incident = list(embedding.neighbors_cw_order(vertex))
        ordered = sorted(incident, key=lambda net_vertex: perimeter_parameter(
            attachments[terminal_for[(ref, net_vertex.split(":", 1)[1])]], components[ref]))
        result[ref] = ordered
    return result


def validate_geometry(circuit, graph, embedding, components, nets, attachments, edges):
    expected = {(terminal, net.name) for net in circuit.nets for terminal in net.terminals}
    actual = {(edge["terminal"], edge["net"]) for edge in edges}
    overlaps = clearance = 0
    for index, first in enumerate(circuit.refs):
        for second in circuit.refs[index + 1:]:
            if gr.overlap_area(components[first], components[second]) > 0:
                overlaps += 1
            elif gr.rect_gap(components[first], components[second]) < gr.BODY_CLEARANCE:
                clearance += 1
    crossings = 0
    for index, first in enumerate(edges):
        for second in edges[index + 1:]:
            if first["net"] != second["net"] and gr.proper_intersection(
                    first["points"][0], first["points"][1],
                    second["points"][0], second["points"][1]):
                crossings += 1
    body_intersections = 0
    body_details = []
    for edge in edges:
        start, end = edge["points"]
        for ref in circuit.refs:
            if ref == edge["component"]:
                continue
            if gr.segment_rect_distance(start, end, gr.rect_bounds(components[ref])) <= EPSILON:
                body_intersections += 1
                body_details.append({"terminal": edge["terminal"], "net": edge["net"],
                                     "unrelated_component": ref})
    geometric_orders = geometric_component_orders(circuit, components, attachments, embedding)
    orientation_rows, degree_three_states = [], []
    for ref in circuit.refs:
        required = list(embedding.neighbors_cw_order(e4p.component_vertex(ref)))
        geometric = geometric_orders[ref]
        if len(required) < 3:
            state = "ambiguous_degree_below_3"
        elif cyclic_equal(geometric, required):
            state = "same"
            degree_three_states.append("same")
        elif cyclic_equal(geometric, list(reversed(required))):
            state = "global_reflection"
            degree_three_states.append("global_reflection")
        else:
            state = "inconsistent"
            degree_three_states.append("inconsistent")
        orientation_rows.append({"component": ref, "required_clockwise": required,
                                 "geometric_perimeter_order": geometric, "state": state})
    unique_states = set(degree_three_states)
    orientation = (next(iter(unique_states)) if len(unique_states) == 1 and
                   unique_states <= {"same", "global_reflection"} else "inconsistent_local_reversal")
    embedding_edges = {frozenset((node, neighbor)) for node in embedding
                       for neighbor in embedding.neighbors_cw_order(node)}
    graph_edges = {frozenset(edge) for edge in graph.edges()}
    return {"graph_vertex_count": graph.number_of_nodes(), "graph_edge_count": graph.number_of_edges(),
            "planar": True, "electrical_incidence_exact": actual == expected,
            "expected_incidence_count": len(expected), "realized_incidence_count": len(actual),
            "component_overlaps": overlaps, "clearance_violations": clearance,
            "different_net_crossings": crossings,
            "unrelated_component_body_intersections": body_intersections,
            "body_intersection_details": body_details,
            "embedding_vertex_set_exact": set(embedding) == set(graph),
            "embedding_edge_set_exact": embedding_edges == graph_edges,
            "cyclic_order_orientation": orientation,
            "cyclic_order_rows": orientation_rows,
            "valid": (actual == expected and overlaps == clearance == crossings == body_intersections == 0
                      and set(embedding) == set(graph) and embedding_edges == graph_edges
                      and orientation in ("same", "global_reflection"))}


def measure(circuit, components, edges):
    length = sum(math.dist(*edge["points"]) for edge in edges)
    length_raw = sum((edge["points"][0][0] - edge["points"][1][0]) ** 2 +
                     (edge["points"][0][1] - edge["points"][1][1]) ** 2 for edge in edges)
    overlaps = overlap_raw = clearance = clearance_raw = 0
    for index, first in enumerate(circuit.refs):
        for second in circuit.refs[index + 1:]:
            area = gr.overlap_area(components[first], components[second])
            if area > 0:
                overlaps += 1; overlap_raw += 1.0 + (area / (gr.WIDTH * gr.HEIGHT)) ** 2
            else:
                gap = gr.rect_gap(components[first], components[second])
                if gap < gr.BODY_CLEARANCE:
                    clearance += 1; clearance_raw += ((gr.BODY_CLEARANCE-gap)/gr.BODY_CLEARANCE)**2
    crossings = 0; crowding_raw = 0.0
    boxes = [gr.segment_bounds(*edge["points"]) for edge in edges]
    body_boxes = {ref: gr.rect_bounds(components[ref]) for ref in circuit.refs}
    for index, first in enumerate(edges):
        a, b = first["points"]
        for ref in circuit.refs:
            if ref != first["component"] and gr.bounds_gap(boxes[index], body_boxes[ref]) < gr.NET_BODY_CLEARANCE:
                distance = gr.segment_rect_distance(a, b, body_boxes[ref])
                if distance < gr.NET_BODY_CLEARANCE:
                    crowding_raw += ((gr.NET_BODY_CLEARANCE-distance)/gr.NET_BODY_CLEARANCE)**2
        for offset, second in enumerate(edges[index + 1:], index + 1):
            if first["net"] == second["net"] or gr.bounds_gap(boxes[index], boxes[offset]) >= gr.NET_NET_CLEARANCE:
                continue
            c, d = second["points"]
            if gr.proper_intersection(a, b, c, d): crossings += 1
            distance = gr.segment_distance(a, b, c, d)
            if distance < gr.NET_NET_CLEARANCE:
                crowding_raw += ((gr.NET_NET_CLEARANCE-distance)/gr.NET_NET_CLEARANCE)**2
    raw = {"length": length_raw, "overlap": overlap_raw, "clearance": clearance_raw,
           "crossing": float(crossings), "crowding": crowding_raw}
    energy = {name: raw[name] * gr.WEIGHTS[name] for name in gr.WEIGHTS}
    return {"total_connection_length": length, "component_overlaps": overlaps,
            "clearance_violations": clearance, "net_crossings": crossings,
            "energy": energy, "objective": sum(energy.values())}


def point_doc(point): return {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}


def write_svg(path, circuit, components, nets, attachments, edges):
    all_points = list(components.values()) + list(nets.values())
    left, top = min(p[0] for p in all_points)-60, min(p[1] for p in all_points)-50
    width, height = max(p[0] for p in all_points)-left+60, max(p[1] for p in all_points)-top+50
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.3f} {top:.3f} {width:.3f} {height:.3f}">',
             f'<rect x="{left:.3f}" y="{top:.3f}" width="{width:.3f}" height="{height:.3f}" fill="white"/>',
             '<g stroke="#777" stroke-width="1" fill="none">']
    for edge in edges:
        a,b=edge["points"]; lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    lines.append('</g><g stroke="#111" stroke-width="1.5" fill="#f5f5f5">')
    for ref in circuit.refs:
        x,y=components[ref]; lines.append(f'<rect x="{x-30:.3f}" y="{y-20:.3f}" width="60" height="40"/>')
    lines.append('</g><g font-family="monospace" font-size="10" text-anchor="middle">')
    for ref in circuit.refs:
        x,y=components[ref]; lines.append(f'<text x="{x:.3f}" y="{y:.3f}">{ref}</text>')
    lines.append('</g><g fill="#d22">')
    for point in attachments.values(): lines.append(f'<circle cx="{point[0]:.3f}" cy="{point[1]:.3f}" r="2.5"/>')
    lines.append('</g><g fill="#1769aa">')
    for name in sorted(nets):
        x,y=nets[name]; lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3.5"/>')
    lines.append('</g></svg>'); path.write_text("\n".join(lines)+"\n",encoding="utf-8")


def run(input_path, output_dir):
    circuit=gr.parse_circuit(input_path); graph, source_emb, int_emb, mapping, reverse, point_positions=stable_embedding(circuit)
    embedding_doc={"graph":e4p.graph_document(graph),"vertex_count":34,"edge_count":44,"planar":True,
                   "stable_integer_relabeling":mapping,
                   "cyclic_neighbor_order_clockwise":{node:list(source_emb.neighbors_cw_order(node)) for node in sorted(graph)},
                   "verification":e4p.verify_embedding(graph,source_emb)[0]}
    embedding_path=output_dir/"node_planar_embedding.json"; gr.write_json(embedding_path,embedding_doc)
    scale_results=[]; selected=None
    for scale in SCALE_SEQUENCE:
        realization=realize(circuit,reverse,point_positions,scale)
        validation=validate_geometry(circuit,graph,source_emb,*realization)
        scale_results.append({"scale":scale,"component_overlaps":validation["component_overlaps"],
                              "clearance_violations":validation["clearance_violations"],
                              "different_net_crossings":validation["different_net_crossings"],
                              "unrelated_component_body_intersections":validation["unrelated_component_body_intersections"],
                              "valid":validation["valid"]})
        if validation["valid"]: selected=(scale,realization,validation); break
    if selected is None: raise RuntimeError("no valid scale in predetermined sequence")
    scale,(components,nets,attachments,edges),validation=selected
    if scale != 72 or not validation["valid"]: raise RuntimeError("scale-72 prototype result did not reproduce")
    validation.update({"selected_scale":scale,"tested_scales":scale_results,
                       "unscaled_point_coordinate_extents":{"min_x":0,"max_x":64,"min_y":0,"max_y":21}})
    validation_path=output_dir/"planar_node_initial_validation.json"; gr.write_json(validation_path,validation)
    coordinates={"circuit":circuit.name,"scale":scale,"component_dimensions":{"width":gr.WIDTH,"height":gr.HEIGHT},
                 "components":{ref:point_doc(components[ref]) for ref in circuit.refs},
                 "net_junctions":{name:point_doc(nets[name]) for name in sorted(nets)},
                 "attachments":{terminal:{**point_doc(attachments[terminal]),"net":next(n.name for n in circuit.nets if terminal in n.terminals)} for terminal in sorted(attachments)},
                 "edges":[{"terminal":e["terminal"],"component":e["component"],"net":e["net"],"points":[point_doc(p) for p in e["points"]]} for e in edges],
                 "component_cyclic_incidence_order":{row["component"]:row["geometric_perimeter_order"] for row in validation["cyclic_order_rows"]}}
    coordinate_path=output_dir/"planar_node_initial_coordinates.json"; gr.write_json(coordinate_path,coordinates)
    metrics=measure(circuit,components,edges)
    comparisons={"experiment_4c_initial":json.loads((output_dir/"closeness_metrics.json").read_text())["closeness_initial"],
                 "experiment_4f_final":json.loads((output_dir/"closeness_hv_annealed_metrics.json").read_text())["states"]["experiment_4f_annealed_and_quenched"],
                 "experiment_4h_junction_only":json.loads((output_dir/"junction_untangled_metrics.json").read_text())["junction_only_final"],
                 "experiment_4i_planar_node_initial":gr.rounded_metrics(metrics)}
    metrics_doc={"metrics":gr.rounded_metrics(metrics),"measurement_only":True,"comparison":comparisons,
                 "construction":{"method":"NetworkX combinatorial_embedding_to_pos after stable integer relabeling; straight incidence rays truncated at rectangles","selected_scale":scale,"scale_sequence":list(SCALE_SEQUENCE),"tested_scales":scale_results}}
    metrics_path=output_dir/"planar_node_initial_metrics.json"; gr.write_json(metrics_path,metrics_doc)
    svg_path=output_dir/"planar_node_initial.svg"; write_svg(svg_path,circuit,components,nets,attachments,edges)
    return {"metrics":metrics_doc,"validation":validation,"hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (embedding_path,coordinate_path,metrics_path,validation_path,svg_path)}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("output_dir",type=Path);args=parser.parse_args()
    print(json.dumps(run(args.input,args.output_dir),sort_keys=True,indent=2))


if __name__=="__main__":main()
