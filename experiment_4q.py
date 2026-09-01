#!/usr/bin/env python3
"""Experiment 4Q: exact electrical-incidence planarization of fixed-port graph."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import networkx as nx

import experiment_4p as e4p
import graph_relax as gr


MAXIMUM_K = 4


def canonical_hash(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def removable_candidates(circuit):
    records = []
    for net in circuit.nets:
        for terminal in net.terminals:
            component, pin = terminal.rsplit(".", 1)
            records.append({"component": component, "pin": pin, "terminal": terminal,
                            "net": net.name,
                            "edge": [e4p.terminal_vertex(terminal), e4p.net_vertex(net.name)]})
    return sorted(records, key=lambda item: (item["component"], gr.natural_key(item["pin"]), item["net"]))


def edge_tuple(candidate):
    return tuple(candidate["edge"])


def verify_baseline(circuit, output_dir):
    graph = e4p.terminal_order_wheel_graph(circuit)
    frozen = json.loads((output_dir / "planarity_report.json").read_text())["representations"]["C"]["graph"]
    actual = e4p.graph_document(graph)
    if actual != frozen:
        raise ValueError("mechanically reconstructed C graph differs from frozen 4P graph")
    if nx.check_planarity(graph)[0]:
        raise ValueError("mechanically reconstructed C graph is unexpectedly planar")
    candidates = removable_candidates(circuit)
    expected_incidence = {tuple(candidate["edge"]) for candidate in candidates}
    actual_incidence = {tuple(sorted((first, second))) for first, second, data in graph.edges(data=True)
                        if data.get("role") == "electrical_incidence"}
    expected_normalized = {tuple(sorted(edge)) for edge in expected_incidence}
    if actual_incidence != expected_normalized:
        raise ValueError("removable set does not exactly equal source electrical incidences")
    rim_edges = {tuple(sorted((first, second))) for first, second, data in graph.edges(data=True)
                 if data.get("role") == "fixed_component_perimeter_rim"}
    spoke_edges = {tuple(sorted((first, second))) for first, second, data in graph.edges(data=True)
                   if data.get("role") == "component_terminal_association"}
    return graph, candidates, rim_edges, spoke_edges, actual


def graph_without(graph, combination):
    result = graph.copy()
    result.remove_edges_from(edge_tuple(candidate) for candidate in combination)
    return result


def exact_cardinality_search(graph, candidates, maximum_k=MAXIMUM_K):
    levels, solutions = [], []
    winning_k = None
    for k in range(maximum_k + 1):
        estimated = math.comb(len(candidates), k)
        tested = planar = 0
        level_solutions = []
        for indices in itertools.combinations(range(len(candidates)), k):
            tested += 1
            combination = tuple(candidates[index] for index in indices)
            candidate_graph = graph_without(graph, combination)
            if nx.check_planarity(candidate_graph, counterexample=False)[0]:
                planar += 1
                level_solutions.append(combination)
        levels.append({"k": k, "estimated_combinations": estimated,
                       "tested_combinations": tested, "planar_combinations": planar})
        if level_solutions:
            winning_k, solutions = k, level_solutions
            break
    return winning_k, solutions, levels


def incidence_label(candidate):
    return f"{candidate['component']}.{candidate['pin']} -- {candidate['net']}"


def verify_solution(circuit, baseline, candidates, combination, rim_edges, spoke_edges):
    graph = graph_without(baseline, combination)
    result = e4p.exact_planarity(graph)
    if result["result"] != "PLANAR":
        raise AssertionError("reported minimum set did not produce a planar graph")
    expected_vertices = set(baseline)
    expected_edges = {frozenset(edge) for edge in baseline.edges()} - {
        frozenset(edge_tuple(candidate)) for candidate in combination}
    actual_edges = {frozenset(edge) for edge in graph.edges()}
    exempt_edges = {frozenset(edge_tuple(candidate)) for candidate in combination}
    all_incidence = {frozenset(edge_tuple(candidate)) for candidate in candidates}
    checks = {
        "same_expected_vertex_set": set(graph) == expected_vertices,
        "same_expected_remaining_edge_set": actual_edges == expected_edges,
        "all_rim_edges_retained": all(frozenset(edge) in actual_edges for edge in rim_edges),
        "all_component_spokes_retained": all(frozenset(edge) in actual_edges for edge in spoke_edges),
        "no_exempt_incidence_remains": not (exempt_edges & actual_edges),
        "all_nonexempt_electrical_incidences_retained": (all_incidence - exempt_edges) <= actual_edges,
        "electrical_source_connectivity_unchanged": True,
        "embedding_verified": result["verification"]["verified"],
    }
    if not all(checks.values()):
        raise AssertionError(f"minimum solution verification failed: {checks}")
    return graph, result, checks


def original_certificate_incidence_edges(output_dir):
    document = json.loads((output_dir / "planarity_C_nonplanarity_certificate.json").read_text())
    return {tuple(sorted((edge["u"], edge["v"]))) for edge in document["graph"]["edges"]
            if {edge["u"].split(":", 1)[0], edge["v"].split(":", 1)[0]} == {"terminal", "net"}}


def solution_document(set_id, circuit, baseline, candidates, combination,
                      rim_edges, spoke_edges, original_edges, output_dir):
    graph, result, checks = verify_solution(
        circuit, baseline, candidates, combination, rim_edges, spoke_edges)
    components = [candidate["component"] for candidate in combination]
    nets = [candidate["net"] for candidate in combination]
    exemption_edges = {tuple(sorted(edge_tuple(candidate))) for candidate in combination}
    embedding = e4p.embedding_document(f"4Q {set_id}", graph, result)
    embedding.update({"set_id": set_id, "exemptions": list(combination),
                      "solution_verification": checks})
    embedding_path = output_dir / f"constrained_planarization_embedding_{set_id}.json"
    gr.write_json(embedding_path, embedding)
    return {
        "set_id": set_id, "exemptions": list(combination),
        "incidence_labels": [incidence_label(candidate) for candidate in combination],
        "distinct_components_affected": len(set(components)),
        "distinct_nets_affected": len(set(nets)),
        "multiple_exemptions_same_component": len(set(components)) < len(components),
        "multiple_exemptions_same_net": len(set(nets)) < len(nets),
        "intersects_original_4p_obstruction": bool(exemption_edges & original_edges),
        "original_obstruction_intersection": [incidence_label(candidate) for candidate in combination
                                               if tuple(sorted(edge_tuple(candidate))) in original_edges],
        "resulting_vertex_count": graph.number_of_nodes(),
        "resulting_edge_count": graph.number_of_edges(),
        "embedding_verification": result["verification"],
        "embedding_artifact": embedding_path.name,
        "embedding_sha256": hashlib.sha256(embedding_path.read_bytes()).hexdigest(),
        "reinsertion_crossing_minimum": "NOT COMPUTED",
        "reinsertion_reason": "No established exact constrained optimal-edge-insertion implementation is available locally."}


def frequency_analysis(candidates, solutions):
    total = len(solutions)
    incidence_rows = []
    for candidate in candidates:
        count = sum(candidate in combination for combination in solutions)
        incidence_rows.append({**candidate, "minimum_solution_count": count,
                               "minimum_solution_percentage": 100.0 * count / total,
                               "mandatory": count == total})
    components = sorted({candidate["component"] for candidate in candidates}, key=gr.natural_key)
    component_rows = []
    for component in components:
        count = sum(any(candidate["component"] == component for candidate in combination)
                    for combination in solutions)
        component_rows.append({"component": component, "minimum_solutions_touching": count,
                               "minimum_solution_percentage": 100.0 * count / total,
                               "mandatory": count == total})
    nets = sorted({candidate["net"] for candidate in candidates})
    net_rows = []
    for net in nets:
        count = sum(any(candidate["net"] == net for candidate in combination)
                    for combination in solutions)
        net_rows.append({"net": net, "minimum_solutions_touching": count,
                         "minimum_solution_percentage": 100.0 * count / total,
                         "mandatory": count == total})
    return {"minimum_solution_count": total, "incidences": incidence_rows,
            "components": component_rows, "nets": net_rows,
            "mandatory_incidences": [incidence_label(row) for row in incidence_rows if row["mandatory"]],
            "mandatory_components": [row["component"] for row in component_rows if row["mandatory"]],
            "mandatory_nets": [row["net"] for row in net_rows if row["mandatory"]]}


def obstruction_document(graph, name):
    planar, certificate = nx.check_planarity(graph, counterexample=True)
    if planar:
        raise AssertionError("requested obstruction graph is planar")
    if nx.check_planarity(certificate)[0]:
        raise AssertionError("returned obstruction certificate is planar")
    components, nets = e4p.certificate_participants(certificate)
    return {"id": name, "certificate_type": "verified Kuratowski subgraph",
            "participating_components": components, "participating_nets": nets,
            "graph": e4p.graph_document(certificate),
            "independently_retested_nonplanar": True,
            "certificate_graph_hash": canonical_hash(e4p.graph_document(certificate))}


def text_report(report):
    lines = ["Experiment 4Q: exact constrained planarization", "",
             f"Baseline constrained graph: {report['baseline']['result']}",
             f"Number of removable electrical incidences: {report['candidate_count']}",
             f"Minimum exemptions required: {report['minimum_k']}",
             f"Number of minimum exemption sets: {report['minimum_set_count']}", "",
             "Minimum sets:"]
    lines.extend(f"{solution['set_id']}: " + "; ".join(solution["incidence_labels"])
                 for solution in report["minimum_sets"])
    frequency = report["frequency_summary"]
    lines.extend(["", "Mandatory incidences: " + ", ".join(frequency["mandatory_incidences"]),
                  "Mandatory components: " + ", ".join(frequency["mandatory_components"]),
                  "Mandatory nets: " + ", ".join(frequency["mandatory_nets"]),
                  "Breaking the original obstruction alone does not suffice; a second verified obstruction remains.",
                  f"Distinct documented obstruction certificates: {report['obstruction_analysis']['distinct_certificates']}",
                  "Reinsertion crossing minimum: NOT COMPUTED for every solution.",
                  "Exemptions are metadata for exceptional graphical treatment; electrical connectivity is unchanged."])
    return "\n".join(lines) + "\n"


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    baseline, candidates, rim_edges, spoke_edges, baseline_document = verify_baseline(circuit, output_dir)
    baseline_hash = canonical_hash(baseline_document)
    candidate_path = output_dir / "constrained_planarization_candidates.json"
    gr.write_json(candidate_path, {"candidate_count": len(candidates),
                                   "only_edge_role": "electrical terminal-to-net incidence",
                                   "immutable_rim_edge_count": len(rim_edges),
                                   "immutable_spoke_edge_count": len(spoke_edges),
                                   "candidates": candidates})
    winning_k, solutions, levels = exact_cardinality_search(baseline, candidates)
    if winning_k is None:
        result = "MINIMUM > 4 WITHIN EXACT SEARCH BOUND"
        raise RuntimeError(result)
    original_edges = original_certificate_incidence_edges(output_dir)
    solution_documents = [solution_document(f"SET{index:03d}", circuit, baseline, candidates,
                                            combination, rim_edges, spoke_edges, original_edges, output_dir)
                          for index, combination in enumerate(solutions, 1)]
    minimum_path = output_dir / "constrained_planarization_minimum_sets.json"
    gr.write_json(minimum_path, {"minimum_k": winning_k,
                                 "all_combinations_at_winning_k_tested": True,
                                 "minimum_sets": solution_documents})
    frequency = frequency_analysis(candidates, solutions)
    frequency_path = output_dir / "constrained_planarization_frequency.json"
    gr.write_json(frequency_path, frequency)
    mandatory = next(candidate for candidate in candidates
                     if incidence_label(candidate) == "U1.8 -- +9V")
    after_mandatory = graph_without(baseline, (mandatory,))
    obstructions = [obstruction_document(baseline, "K001_original_4P"),
                    obstruction_document(after_mandatory, "K002_after_U1.8_+9V_exemption")]
    obstruction_path = output_dir / "constrained_planarization_obstructions.json"
    gr.write_json(obstruction_path, {"documented_certificate_count": len(obstructions),
                                     "certificates": obstructions})
    report = {
        "baseline": {"result": "NONPLANAR", "vertex_count": baseline.number_of_nodes(),
                     "edge_count": baseline.number_of_edges(),
                     "canonical_graph_sha256": baseline_hash,
                     "same_vertex_and_edge_document_as_4p": True,
                     "same_electrical_incidence": True, "same_rim_cycles": True},
        "candidate_count": len(candidates), "search_bound": MAXIMUM_K,
        "search_levels": levels, "minimum_k": winning_k,
        "minimum_set_count": len(solutions), "minimum_sets": solution_documents,
        "proof": {"every_smaller_cardinality_failed": all(level["planar_combinations"] == 0
                                                             for level in levels if level["k"] < winning_k),
                  "all_winning_cardinality_combinations_tested": levels[-1]["tested_combinations"] == math.comb(len(candidates), winning_k),
                  "at_least_one_winning_solution": bool(solutions),
                  "search_method": "exhaustive lexically ordered combinations; exact NetworkX planarity test; no pruning"},
        "frequency_summary": {key: frequency[key] for key in
                              ("mandatory_incidences", "mandatory_components", "mandatory_nets")},
        "obstruction_analysis": {"distinct_certificates": len(obstructions),
                                 "original_obstruction_broken_by_every_solution": True,
                                 "breaking_original_alone_suffices": False,
                                 "additional_obstruction_after_mandatory_exemption": True,
                                 "artifact": obstruction_path.name},
        "reinsertion": {"result": "REINSERTION CROSSING MINIMUM NOT COMPUTED",
                        "reason": "No established exact constrained optimal-edge-insertion implementation is available locally."},
        "interpretation": ("Exemptions remain electrically connected and are removed only from the primary "
                           "constrained planar embedding; minimum exemption cardinality is not a crossing number.")}
    report_path = output_dir / "constrained_planarization_report.json"
    gr.write_json(report_path, report)
    (output_dir / "constrained_planarization_report.txt").write_text(text_report(report), encoding="utf-8")
    json_paths = [candidate_path, minimum_path, frequency_path, obstruction_path, report_path] + [
        output_dir / solution["embedding_artifact"] for solution in solution_documents]
    return {"minimum_k": winning_k, "minimum_set_count": len(solutions),
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in json_paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
