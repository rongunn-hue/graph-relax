#!/usr/bin/env python3
"""Experiment 4P: exact ordinary planarity certification for IAMP."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import networkx as nx

import graph_relax as gr


def component_vertex(ref):
    return f"component:{ref}"


def net_vertex(name):
    return f"net:{name}"


def terminal_vertex(name):
    return f"terminal:{name}"


def incidence_graph(circuit):
    graph = nx.Graph()
    graph.add_nodes_from((component_vertex(ref), {"kind": "component", "name": ref})
                         for ref in circuit.refs)
    graph.add_nodes_from((net_vertex(net.name), {"kind": "net", "name": net.name})
                         for net in circuit.nets)
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            graph.add_edge(component_vertex(ref), net_vertex(net.name), terminal=terminal)
    return graph


def terminal_expanded_graph(circuit):
    graph = nx.Graph()
    graph.add_nodes_from((component_vertex(ref), {"kind": "component_core", "name": ref})
                         for ref in circuit.refs)
    graph.add_nodes_from((net_vertex(net.name), {"kind": "net_junction", "name": net.name})
                         for net in circuit.nets)
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            vertex = terminal_vertex(terminal)
            graph.add_node(vertex, kind="terminal", name=terminal, component=ref, net=net.name)
            graph.add_edge(component_vertex(ref), vertex, role="component_terminal_association")
            graph.add_edge(vertex, net_vertex(net.name), role="electrical_incidence")
    return graph


def graph_document(graph):
    return {
        "vertices": [{"id": node, **dict(sorted(graph.nodes[node].items()))}
                     for node in sorted(graph)],
        "edges": [{"u": first, "v": second, **dict(sorted(data.items()))}
                  for first, second, data in sorted(
                      ((min(a, b), max(a, b), data) for a, b, data in graph.edges(data=True)),
                      key=lambda item: (item[0], item[1]))],
    }


def cyclic_terminal_orders(circuit):
    # local_terminal traverses top, right, bottom, left: clockwise in SVG coordinates.
    return {component_vertex(ref): [terminal_vertex(f"{ref}.{pin}") for pin in circuit.pins[ref]]
            for ref in circuit.refs}


def terminal_order_wheel_graph(circuit):
    """Necessary wheel-gadget expansion for disk-terminal embeddings.

    Rim cycles are invariant under reversal, so this graph relaxes C by allowing
    every component order to reverse independently. Nonplanarity of the relaxed
    graph therefore proves C nonplanar and also rules out all D orientations.
    """
    graph = terminal_expanded_graph(circuit)
    for ref in circuit.refs:
        terminals = [terminal_vertex(f"{ref}.{pin}") for pin in circuit.pins[ref]]
        if len(terminals) < 3:
            continue
        for index, first in enumerate(terminals):
            second = terminals[(index + 1) % len(terminals)]
            graph.add_edge(first, second, role="fixed_component_perimeter_rim")
    return graph


def embedding_faces(embedding):
    visited = set()
    faces = []
    for first in sorted(embedding):
        for second in embedding.neighbors_cw_order(first):
            if (first, second) in visited:
                continue
            face = embedding.traverse_face(first, second, visited)
            faces.append(face)
    return faces


def verify_embedding(graph, embedding):
    embedding.check_structure()
    graph_edges = {frozenset(edge) for edge in graph.edges()}
    embedding_edges = {frozenset((node, neighbor)) for node in embedding
                       for neighbor in embedding.neighbors_cw_order(node)}
    same_vertices = set(graph) == set(embedding)
    same_edges = graph_edges == embedding_edges
    faces = embedding_faces(embedding)
    connected = nx.number_connected_components(graph)
    euler_left = graph.number_of_nodes() - graph.number_of_edges() + len(faces)
    euler_right = 1 + connected
    return {"embedding_structure_valid": True, "same_vertex_set": same_vertices,
            "same_edge_set": same_edges, "face_count": len(faces),
            "euler_characteristic": euler_left, "required_euler_characteristic": euler_right,
            "euler_verified": euler_left == euler_right,
            "verified": same_vertices and same_edges and euler_left == euler_right}, faces


def exact_planarity(graph):
    planar, certificate = nx.check_planarity(graph, counterexample=True)
    if planar:
        verification, faces = verify_embedding(graph, certificate)
        if not verification["verified"]:
            raise AssertionError("NetworkX embedding failed independent verification")
        return {"result": "PLANAR", "embedding": certificate,
                "verification": verification, "faces": faces}
    # NetworkX returns a Kuratowski subgraph when counterexample=True.
    planar_again, _ = nx.check_planarity(certificate)
    if planar_again:
        raise AssertionError("NetworkX nonplanarity counterexample is planar")
    return {"result": "NONPLANAR", "certificate": certificate,
            "verification": {"certificate_independently_retested_nonplanar": True}}


def embedding_document(name, graph, result):
    embedding = result["embedding"]
    return {"representation": name,
            "cyclic_edge_order_clockwise": {
                node: list(embedding.neighbors_cw_order(node)) for node in sorted(embedding)},
            "faces": result["faces"], "verification": result["verification"],
            "graph": graph_document(graph)}


def representation_summary(graph, result):
    return {"vertex_count": graph.number_of_nodes(), "edge_count": graph.number_of_edges(),
            "connected_components": [sorted(group) for group in
                                     sorted(nx.connected_components(graph), key=lambda group: min(group))],
            "connected_component_count": nx.number_connected_components(graph),
            "result": result["result"], "graph": graph_document(graph),
            "verification": result["verification"]}


def certificate_participants(certificate):
    components, nets = set(), set()
    for node in certificate:
        prefix, name = node.split(":", 1)
        if prefix == "component":
            components.add(name)
        elif prefix == "terminal":
            components.add(name.rsplit(".", 1)[0])
        elif prefix == "net":
            nets.add(name)
    return sorted(components), sorted(nets)


def incidence_preserved(circuit, graph):
    expected = {(terminal, net.name) for net in circuit.nets for terminal in net.terminals}
    actual = set()
    for terminal, net_name in expected:
        vertex = terminal_vertex(terminal)
        if vertex in graph and graph.has_edge(vertex, net_vertex(net_name)):
            ref = terminal.rsplit(".", 1)[0]
            if graph.has_edge(component_vertex(ref), vertex):
                actual.add((terminal, net_name))
    return expected == actual


def text_report(report):
    a, b, c, d = (report["representations"][key] for key in ("A", "B", "C", "D"))
    lines = ["Experiment 4P: exact planarity certification", "",
             "Representation | Constraint level | Result",
             f"A | component/net incidence | {a['result']}",
             f"B | explicit terminals | {b['result']}",
             f"C | fixed terminal cyclic order | {c['result']}",
             f"D | permitted orientation freedom | {d['result']}", "",
             "Answers", "",
             f"1. Abstract incidence graph planar: {'yes' if a['result'] == 'PLANAR' else 'no'}.",
             "2. Zero crossings exists at abstract topology level: yes, using a planar embedding with curved/polyline edges.",
             f"3. Terminal expansion changed planarity: {'no' if a['result'] == b['result'] else 'yes'}.",
             "4. Fixed terminal cyclic-order result: NONPLANAR. A necessary wheel-gadget relaxation is already nonplanar.",
             "5. Orientation freedom does not restore planarity: the wheel relaxation permits independent reversal and remains nonplanar.",
             "6. Planarization: not applicable because A is planar.",
             "7. Next experiment must acknowledge a terminal-order obstruction; zero crossings requires changing the representation or special graphical/layer treatment.", "",
             "Important distinctions", "",
             "Graph planarity does not imply that a straight-line drawing with fixed coordinates is crossing-free.",
             "A planar embedding may require curved or polyline edges.",
             "Fixed terminal cyclic order is an additional constraint not decided by ordinary graph planarity.",
             "Minimum-crossing drawing and planarization by edge deletion are separate problems.", "",
             f"Library: NetworkX {report['implementation']['library_version']}",
             "Function: networkx.check_planarity(..., counterexample=True), Left-Right Planarity Test."]
    return "\n".join(lines) + "\n"


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    graph_a = incidence_graph(circuit)
    graph_b = terminal_expanded_graph(circuit)
    graph_c = terminal_order_wheel_graph(circuit)
    result_a = exact_planarity(graph_a)
    result_b = exact_planarity(graph_b)
    result_c = exact_planarity(graph_c)
    if result_a["result"] != "PLANAR" or result_b["result"] != "PLANAR":
        raise AssertionError("IAMP ordinary representations unexpectedly changed expected result")
    embedding_a_path = output_dir / "planarity_A_embedding.json"
    embedding_b_path = output_dir / "planarity_B_embedding.json"
    gr.write_json(embedding_a_path, embedding_document("A", graph_a, result_a))
    gr.write_json(embedding_b_path, embedding_document("B", graph_b, result_b))
    if result_c["result"] != "NONPLANAR":
        raise AssertionError("terminal-order wheel relaxation unexpectedly planar")
    certificate_path = output_dir / "planarity_C_nonplanarity_certificate.json"
    certificate = result_c["certificate"]
    certificate_components, certificate_nets = certificate_participants(certificate)
    gr.write_json(certificate_path, {
        "representation": "C wheel-gadget relaxation",
        "certificate_type": "Kuratowski subgraph returned by NetworkX",
        "graph": graph_document(certificate),
        "independent_verification": result_c["verification"],
        "participating_components": certificate_components,
        "participating_nets": certificate_nets})
    orientation_components = sorted(ref for ref in circuit.refs if len(circuit.pins[ref]) >= 3)
    orientation_combinations = 2 ** len(orientation_components)
    d_path = output_dir / "planarity_D_result.json"
    gr.write_json(d_path, {
        "result": "NONPLANAR",
        "unique_topological_states_per_component": {
            ref: (["current_or_180", "horizontal_or_vertical_mirror"]
                  if ref in orientation_components else ["all_transformations_equivalent"])
            for ref in circuit.refs},
        "components_with_two_states": orientation_components,
        "unique_orientation_combinations": orientation_combinations,
        "individually_enumerated": 0,
        "collectively_ruled_out_by_nonplanar_reversal_relaxation": orientation_combinations,
        "proof": ("The wheel rim is unchanged when its cyclic order is reversed. Thus the ordinary "
                  "wheel-augmented graph permits each component to choose either topological state "
                  "independently. Its verified nonplanarity rules out every assignment at once."),
        "certificate_artifact": certificate_path.name})
    constraints = cyclic_terminal_orders(circuit)
    report = {
        "input": str(input_path),
        "implementation": {"library": "NetworkX", "library_version": nx.__version__,
                           "function": "networkx.check_planarity(graph, counterexample=True)",
                           "algorithm": "Left-Right Planarity Test",
                           "ordinary_planarity_exact": True},
        "representations": {
            "A": {"constraint_level": "component/net bipartite incidence",
                  **representation_summary(graph_a, result_a),
                  "singleton_net_incidence": {net.name: list(net.terminals) for net in circuit.nets
                                               if len(net.terminals) == 1},
                  "embedding_artifact": embedding_a_path.name},
            "B": {"constraint_level": "explicit terminals; unconstrained component rotation",
                  **representation_summary(graph_b, result_b),
                  "two_terminal_rule": "terminal -- net junction -- terminal",
                  "electrical_incidence_preserved": incidence_preserved(circuit, graph_b),
                  "embedding_artifact": embedding_b_path.name},
            "C": {"constraint_level": "fixed clockwise component terminal cyclic orders",
                  **representation_summary(graph_c, result_c),
                  "fixed_terminal_orders": constraints,
                  "method": ("Add a rim cycle through each component's 3+ terminals in fixed cyclic order, "
                             "retaining the component-core spokes. Any disk embedding induces a planar "
                             "embedding of this graph. The gadget allows independent rim reversal, making "
                             "it a relaxation; nonplanarity of the relaxation is a proof for C."),
                  "certificate_artifact": certificate_path.name},
            "D": {"constraint_level": "current/180-degree/horizontal-mirror/vertical-mirror states",
                  "result": "NONPLANAR",
                  "topological_state_note": ("For cyclic order, current and 180-degree rotation are equivalent; "
                                             "horizontal and vertical mirrors are the reversed order."),
                  "components_with_two_unique_states": orientation_components,
                  "unique_orientation_combinations": orientation_combinations,
                  "combinations_collectively_pruned_by_certificate": orientation_combinations,
                  "result_artifact": d_path.name},
        },
        "comparison": {"terminal_expansion_changed_planarity": False,
                       "reason": "B subdivides each A incidence edge and adds the corresponding component-core star without imposing rotations"},
        "planarization": {"performed": False, "reason": "Representation A is planar"},
        "conclusions": {"abstract_incidence_graph_planar": True,
                        "zero_crossing_abstract_embedding_exists": True,
                        "fixed_terminal_order_decided": True,
                        "representation_level_making_zero_crossing_impossible": "C: component terminal cyclic-order constraints"}}
    report_path = output_dir / "planarity_report.json"
    gr.write_json(report_path, report)
    (output_dir / "planarity_report.txt").write_text(text_report(report), encoding="utf-8")
    return {"results": {key: value["result"] for key, value in report["representations"].items()},
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in (report_path, embedding_a_path, embedding_b_path,
                                    certificate_path, d_path)}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
