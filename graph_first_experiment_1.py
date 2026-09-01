#!/usr/bin/env python3
"""Graph-first Experiment 1: authoritative netlist to abstract graph."""

from __future__ import annotations

import argparse
import json
import math
import shlex
from dataclasses import dataclass
from pathlib import Path

import networkx as nx


@dataclass(frozen=True)
class Component:
    name: str
    device_id: str
    source_line: int


@dataclass(frozen=True)
class Incidence:
    component: str
    terminal: str
    pin: str
    net: str
    source_line: int
    source_position: int


@dataclass(frozen=True)
class ParsedCircuit:
    name: str
    source: str
    components: tuple[Component, ...]
    nets: tuple[str, ...]
    incidences: tuple[Incidence, ...]


def component_id(name: str) -> str:
    return f"component:{name}"


def net_id(name: str) -> str:
    return f"net:{name}"


def parse_circuit(path: Path) -> ParsedCircuit:
    """Parse electrical declarations without importing placement semantics."""
    title = None
    components: dict[str, Component] = {}
    nets: list[str] = []
    incidences: list[Incidence] = []
    used_terminals: set[str] = set()
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        words = shlex.split(raw, comments=True)
        if not words:
            continue
        if words[0] == "circuit" and len(words) == 2:
            if title is not None:
                raise ValueError(f"line {line_number}: duplicate circuit declaration")
            title = words[1]
        elif words[0] == "component" and len(words) == 3:
            name, device_id = words[1:]
            if name in components:
                raise ValueError(f"line {line_number}: duplicate component {name}")
            components[name] = Component(name, device_id, line_number)
        elif words[0] == "net" and len(words) >= 3:
            name = words[1]
            if name in nets:
                raise ValueError(f"line {line_number}: duplicate net {name}")
            nets.append(name)
            for position, terminal in enumerate(words[2:], 1):
                if "." not in terminal:
                    raise ValueError(f"line {line_number}: malformed terminal {terminal}")
                ref, pin = terminal.rsplit(".", 1)
                if ref not in components:
                    raise ValueError(f"line {line_number}: unknown component {ref}")
                if terminal in used_terminals:
                    raise ValueError(f"line {line_number}: terminal used more than once: {terminal}")
                used_terminals.add(terminal)
                incidences.append(Incidence(ref, terminal, pin, name, line_number, position))
        else:
            raise ValueError(f"line {line_number}: unsupported statement")
    if title is None:
        raise ValueError("missing circuit declaration")
    return ParsedCircuit(
        title,
        path.as_posix(),
        tuple(sorted(components.values(), key=lambda item: item.name)),
        tuple(sorted(nets)),
        tuple(sorted(incidences, key=lambda item: (item.component, item.net, item.terminal))),
    )


def build_graph(circuit: ParsedCircuit) -> nx.Graph:
    graph = nx.Graph()
    for component in circuit.components:
        graph.add_node(component_id(component.name), type="COMPONENT", name=component.name)
    for name in circuit.nets:
        graph.add_node(net_id(name), type="NET", name=name)
    for item in circuit.incidences:
        left, right = component_id(item.component), net_id(item.net)
        if graph.has_edge(left, right):
            raise ValueError(
                f"multiple terminal incidences collapse onto graph edge {left} -- {right}"
            )
        graph.add_edge(
            left,
            right,
            component=item.component,
            net=item.net,
            pin=item.pin,
            terminal=item.terminal,
            source_line=item.source_line,
            source_position=item.source_position,
        )
    return graph


def sorted_edges(graph: nx.Graph):
    return sorted((min(a, b), max(a, b), data) for a, b, data in graph.edges(data=True))


def graph_document(circuit: ParsedCircuit, graph: nx.Graph) -> dict:
    return {
        "schema": "graph-relax.abstract-electrical-graph.v1",
        "source_circuit": circuit.source,
        "circuit_name": circuit.name,
        "device_id_policy": "parsed from component declarations; excluded from graph topology",
        "source_components": [
            {
                "name": item.name,
                "device_id": item.device_id,
                "source_line": item.source_line,
            }
            for item in circuit.components
        ],
        "vertices": [
            {"id": node, "name": graph.nodes[node]["name"], "type": graph.nodes[node]["type"]}
            for node in sorted(graph)
        ],
        "edges": [
            {"source": left, "target": right, **dict(sorted(data.items()))}
            for left, right, data in sorted_edges(graph)
        ],
    }


def verify_embedding(graph: nx.Graph, embedding: nx.PlanarEmbedding) -> dict:
    embedding.check_structure()
    expected_edges = {frozenset(edge) for edge in graph.edges}
    actual_edges = {
        frozenset((node, neighbor))
        for node in embedding
        for neighbor in embedding.neighbors_cw_order(node)
    }
    return {
        "same_vertices": set(graph) == set(embedding),
        "same_edges": expected_edges == actual_edges,
        "embedding_structure_valid": True,
        "valid": set(graph) == set(embedding) and expected_edges == actual_edges,
    }


def analyze(circuit: ParsedCircuit, graph: nx.Graph) -> dict:
    expected = {(item.component, item.net, item.terminal) for item in circuit.incidences}
    actual = {
        (data["component"], data["net"], data["terminal"])
        for _, _, data in graph.edges(data=True)
    }
    duplicate_source = len(expected) != len(circuit.incidences)
    component_nodes = {component_id(item.name) for item in circuit.components}
    net_nodes = {net_id(name) for name in circuit.nets}
    partitions_valid = nx.is_bipartite(graph) and all(
        (a in component_nodes and b in net_nodes) or (b in component_nodes and a in net_nodes)
        for a, b in graph.edges
    )
    connected = [sorted(group) for group in nx.connected_components(graph)]
    connected.sort(key=lambda group: (group[0], len(group), group))
    distances = {
        source: {target: distance for target, distance in sorted(lengths.items())}
        for source, lengths in sorted(nx.all_pairs_shortest_path_length(graph))
    }
    closeness = nx.closeness_centrality(graph)
    degrees = dict(graph.degree())
    max_degree = max(degrees.values(), default=0)
    max_closeness = max(closeness.values(), default=0.0)
    highest_degree = sorted(node for node, value in degrees.items() if value == max_degree)
    highest_closeness = sorted(
        node for node, value in closeness.items() if math.isclose(value, max_closeness, abs_tol=1e-15)
    )

    # Insertions are deterministic, so NetworkX's exact Left-Right test returns
    # a reproducible rotation system for this fixed graph and library version.
    planar, certificate = nx.check_planarity(graph, counterexample=True)
    planarity: dict = {"result": "PLANAR" if planar else "NONPLANAR"}
    if planar:
        verification = verify_embedding(graph, certificate)
        if not verification["valid"]:
            raise AssertionError("planarity embedding verification failed")
        planarity["implementation"] = {
            "library": "NetworkX",
            "version": nx.__version__,
            "function": "networkx.check_planarity(counterexample=True)",
            "algorithm": "Left-Right Planarity Test",
        }
        planarity["embedding_cyclic_neighbor_order"] = {
            node: list(certificate.neighbors_cw_order(node)) for node in sorted(certificate)
        }
        planarity["embedding_validation"] = verification
    else:
        if nx.check_planarity(certificate)[0]:
            raise AssertionError("nonplanarity certificate independently tested planar")
        planarity["certificate"] = {
            "vertices": sorted(certificate),
            "edges": [[a, b] for a, b, _ in sorted_edges(certificate)],
            "independently_retested_nonplanar": True,
        }

    bridges = sorted([sorted(edge) for edge in nx.bridges(graph)])
    return {
        "schema": "graph-relax.abstract-electrical-analysis.v1",
        "source_circuit": circuit.source,
        "counts": {
            "component_vertices": len(component_nodes),
            "net_vertices": len(net_nodes),
            "total_vertices": graph.number_of_nodes(),
            "electrical_incidence_edges": graph.number_of_edges(),
            "connected_components": len(connected),
        },
        "validation": {
            "bipartite": nx.is_bipartite(graph),
            "edges_cross_partitions_only": partitions_valid,
            "all_parsed_incidences_appear_exactly_once": expected == actual and not duplicate_source,
            "parsed_incidence_count": len(circuit.incidences),
            "graph_incidence_count": graph.number_of_edges(),
            "duplicate_electrical_incidences": duplicate_source,
            "self_loops": sorted(nx.nodes_with_selfloops(graph)),
            "isolated_vertices": sorted(nx.isolates(graph)),
        },
        "connected_component_members": connected,
        "connected_neighborhoods": {
            node: sorted(graph.neighbors(node)) for node in sorted(graph)
        },
        "degree": {node: degrees[node] for node in sorted(graph)},
        "shortest_path_distances": distances,
        "closeness_centrality": {node: closeness[node] for node in sorted(graph)},
        "structural_centers": {
            "highest_degree": {"degree": max_degree, "vertices": highest_degree},
            "highest_closeness": {"value": max_closeness, "vertices": highest_closeness},
        },
        "articulation_vertices": sorted(nx.articulation_points(graph)),
        "bridges": bridges,
        "planarity": planarity,
    }


def write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(input_path: Path, output_dir: Path) -> tuple[dict, dict]:
    circuit = parse_circuit(input_path)
    graph = build_graph(circuit)
    graph_doc = graph_document(circuit, graph)
    analysis_doc = analyze(circuit, graph)
    counts = analysis_doc["counts"]
    expected_cross_check = (19, 15, 34, 44)
    actual = (
        counts["component_vertices"], counts["net_vertices"],
        counts["total_vertices"], counts["electrical_incidence_edges"],
    )
    if actual != expected_cross_check:
        raise RuntimeError(f"independent counts {actual} differ from prior cross-check {expected_cross_check}")
    write_json(output_dir / "iamp_graph.json", graph_doc)
    write_json(output_dir / "iamp_graph_analysis.json", analysis_doc)
    return graph_doc, analysis_doc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("circuit_examples/iamp.circuit"))
    parser.add_argument("--output", type=Path, default=Path("output/graph_first"))
    args = parser.parse_args()
    _, analysis = run(args.input, args.output)
    print(json.dumps({
        "counts": analysis["counts"],
        "planarity": analysis["planarity"]["result"],
        "highest_degree": analysis["structural_centers"]["highest_degree"],
        "highest_closeness": analysis["structural_centers"]["highest_closeness"],
        "articulation_vertices": analysis["articulation_vertices"],
        "bridges": analysis["bridges"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
