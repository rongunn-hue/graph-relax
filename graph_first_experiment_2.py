#!/usr/bin/env python3
"""Graph-first Experiment 2: coordinate-free topological placement blueprint."""

from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path

import networkx as nx


def load_authoritative_graph(graph_path: Path, analysis_path: Path) -> tuple[nx.Graph, dict, dict]:
    graph_doc = json.loads(graph_path.read_text(encoding="utf-8"))
    analysis_doc = json.loads(analysis_path.read_text(encoding="utf-8"))
    graph = nx.Graph()
    for vertex in graph_doc["vertices"]:
        graph.add_node(vertex["id"], type=vertex["type"], name=vertex["name"])
    for edge in graph_doc["edges"]:
        metadata = {key: value for key, value in edge.items() if key not in {"source", "target"}}
        graph.add_edge(edge["source"], edge["target"], **metadata)
    if set(graph) != {item["id"] for item in graph_doc["vertices"]}:
        raise ValueError("authoritative graph vertex mismatch")
    if graph.number_of_edges() != len(graph_doc["edges"]):
        raise ValueError("authoritative graph edge mismatch")
    return graph, graph_doc, analysis_doc


def select_structural_center(graph: nx.Graph) -> tuple[str, dict, dict]:
    degree = dict(graph.degree())
    closeness = nx.closeness_centrality(graph)
    center = min(graph, key=lambda node: (-closeness[node], -degree[node], node))
    return center, degree, closeness


def restore_embedding(graph: nx.Graph, analysis: dict) -> nx.PlanarEmbedding:
    stored = analysis["planarity"]
    if stored["result"] != "PLANAR" or not stored["embedding_validation"]["valid"]:
        raise ValueError("Experiment 1 does not contain a validated planar embedding")
    rotations = stored["embedding_cyclic_neighbor_order"]
    embedding = nx.PlanarEmbedding()
    embedding.set_data({node: list(rotations[node]) for node in sorted(rotations)})
    embedding.check_structure()
    expected = {frozenset(edge) for edge in graph.edges}
    actual = {
        frozenset((node, neighbor))
        for node in embedding
        for neighbor in embedding.neighbors_cw_order(node)
    }
    if set(graph) != set(embedding) or expected != actual:
        raise ValueError("stored rotation system does not match authoritative graph")
    return embedding


def enumerate_faces(embedding: nx.PlanarEmbedding) -> list[list[str]]:
    visited: set[tuple[str, str]] = set()
    faces: list[list[str]] = []
    for first in sorted(embedding):
        for second in embedding.neighbors_cw_order(first):
            if (first, second) not in visited:
                faces.append(list(embedding.traverse_face(first, second, visited)))
    # Face IDs must not depend on the directed half-edge at which traversal began.
    def canonical_cycle(walk: list[str]) -> tuple[str, ...]:
        rotations = [tuple(walk[index:] + walk[:index]) for index in range(len(walk))]
        return min(rotations)
    return [list(face) for face in sorted(faces, key=canonical_cycle)]


def biconnected_structure(graph: nx.Graph, center: str, distances: dict[str, int]) -> dict:
    articulation = sorted(nx.articulation_points(graph))
    articulation_set = set(articulation)
    edge_blocks = []
    for edges in nx.biconnected_component_edges(graph):
        normalized_edges = sorted([sorted(edge) for edge in edges])
        members = sorted({node for edge in normalized_edges for node in edge})
        edge_blocks.append((members, normalized_edges))
    edge_blocks.sort(key=lambda item: (item[0], item[1]))
    blocks = []
    for index, (members, edges) in enumerate(edge_blocks, 1):
        block_id = f"block:{index:03d}"
        blocks.append({
            "id": block_id,
            "member_vertices": members,
            "articulation_vertices": sorted(set(members) & articulation_set),
            "edge_count": len(edges),
            "edges": edges,
            "bridge_block": len(edges) == 1,
            "distance_from_structural_center": min(distances[node] for node in members),
            "contains_structural_center": center in members,
        })

    tree = nx.Graph()
    tree.add_nodes_from(block["id"] for block in blocks)
    tree.add_nodes_from(f"articulation:{node}" for node in articulation)
    for block in blocks:
        for node in block["articulation_vertices"]:
            tree.add_edge(block["id"], f"articulation:{node}")
    root_blocks = sorted(block["id"] for block in blocks if block["contains_structural_center"])
    if not root_blocks:
        raise AssertionError("structural center is absent from block decomposition")
    root = root_blocks[0]
    tree_distance = nx.single_source_shortest_path_length(tree, root)
    for block in blocks:
        block["block_cut_tree_distance"] = tree_distance[block["id"]]
        block["leaf_block"] = tree.degree(block["id"]) <= 1 and block["id"] != root
    leaf_blocks = sorted(block["id"] for block in blocks if block["leaf_block"])

    # Root the block-cut tree to expose complete outward pendant subtrees.
    parents = {root: None}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for neighbor in sorted(tree.neighbors(node)):
            if neighbor not in parents:
                parents[neighbor] = node
                queue.append(neighbor)
    pendant_branches = []
    for articulation_node in sorted(node for node in tree if node.startswith("articulation:")):
        parent = parents[articulation_node]
        for child in sorted(tree.neighbors(articulation_node)):
            if child == parent:
                continue
            subtree = nx.node_connected_component(nx.subgraph_view(
                tree, filter_edge=lambda a, b, cut=frozenset((articulation_node, child)):
                    frozenset((a, b)) != cut
            ), child)
            branch_blocks = sorted(node for node in subtree if node.startswith("block:"))
            members = sorted({vertex for block in blocks if block["id"] in branch_blocks
                              for vertex in block["member_vertices"]} -
                             {articulation_node.removeprefix("articulation:")})
            pendant_branches.append({
                "attachment_articulation": articulation_node.removeprefix("articulation:"),
                "root_block": child,
                "blocks": branch_blocks,
                "vertices_beyond_articulation": members,
            })
    pendant_branches.sort(key=lambda item: (item["attachment_articulation"], item["root_block"]))
    return {
        "blocks": blocks,
        "block_cut_tree": {
            "root_block": root,
            "nodes": ([{"id": block["id"], "type": "BLOCK"} for block in blocks] +
                      [{"id": f"articulation:{node}", "type": "ARTICULATION", "vertex": node}
                       for node in articulation]),
            "edges": [sorted(edge) for edge in sorted(tree.edges())],
        },
        "articulation_vertices": articulation,
        "bridges": sorted([sorted(edge) for edge in nx.bridges(graph)]),
        "central_blocks": root_blocks,
        "leaf_blocks": leaf_blocks,
        "pendant_branches": pendant_branches,
    }


def face_structure(embedding: nx.PlanarEmbedding, center: str, distances: dict[str, int],
                   closeness: dict[str, float]) -> dict:
    walks = enumerate_faces(embedding)
    faces = []
    for index, walk in enumerate(walks, 1):
        faces.append({
            "id": f"face:{index:03d}",
            "boundary_walk": walk,
            "boundary_vertex_count": len(walk),
            "distinct_boundary_vertex_count": len(set(walk)),
            "contains_structural_center": center in walk,
            "average_graph_distance": sum(distances[node] for node in walk) / len(walk),
            "average_closeness": sum(closeness[node] for node in walk) / len(walk),
        })

    eligible = [face for face in faces if not face["contains_structural_center"]] or faces
    distance_ranking = [face["id"] for face in sorted(
        eligible, key=lambda face: (-face["average_graph_distance"], face["id"]))]
    closeness_ranking = [face["id"] for face in sorted(
        eligible, key=lambda face: (face["average_closeness"], face["id"]))]

    # Unweighted Pareto peeling: farther is better and lower closeness is better.
    remaining = list(eligible)
    tiers = []
    while remaining:
        nondominated = []
        for candidate in remaining:
            dominated = any(
                other["average_graph_distance"] >= candidate["average_graph_distance"]
                and other["average_closeness"] <= candidate["average_closeness"]
                and (other["average_graph_distance"] > candidate["average_graph_distance"]
                     or other["average_closeness"] < candidate["average_closeness"])
                for other in remaining
            )
            if not dominated:
                nondominated.append(candidate)
        nondominated.sort(key=lambda face: face["id"])
        tiers.append([face["id"] for face in nondominated])
        ids = {face["id"] for face in nondominated}
        remaining = [face for face in remaining if face["id"] not in ids]
    return {
        "faces": faces,
        "outer_face_selection_method": {
            "first_filter": "prefer faces not containing the structural center when any exist",
            "criteria": ["maximize average graph distance", "minimize average closeness"],
            "combination": "unweighted Pareto dominance; deterministic face ID only orders ties within a tier",
        },
        "candidate_face_ids": [face["id"] for face in eligible],
        "distance_criterion_ranking": distance_ranking,
        "closeness_criterion_ranking": closeness_ranking,
        "criteria_agree_on_first": distance_ranking[0] == closeness_ranking[0],
        "pareto_tiers": tiers,
        "best_outer_face_candidates": tiers[0],
    }


def sector_structure(graph: nx.Graph, embedding: nx.PlanarEmbedding, center: str,
                     distances: dict[str, int], blocks: dict) -> dict:
    neighbors = list(embedding.neighbors_cw_order(center))
    articulation = set(blocks["articulation_vertices"])
    peripheral = {vertex for branch in blocks["pendant_branches"]
                  for vertex in branch["vertices_beyond_articulation"]}
    assignments = {}
    sectors = {neighbor: [] for neighbor in neighbors}
    for vertex in sorted(graph):
        if vertex == center:
            assignments[vertex] = {"classification": "center", "sectors": []}
            continue
        first_hops = [neighbor for neighbor in neighbors
                      if nx.shortest_path_length(graph, neighbor, vertex) + 1 == distances[vertex]]
        for neighbor in first_hops:
            sectors[neighbor].append(vertex)
        if vertex in articulation:
            classification = "articulation/branch structure"
        elif vertex in peripheral:
            classification = "peripheral"
        elif len(first_hops) == 1:
            classification = "sector-local"
        else:
            classification = "shared between sectors"
        assignments[vertex] = {"classification": classification, "sectors": first_hops}
    shared = sorted(vertex for vertex, item in assignments.items()
                    if item["classification"] == "shared between sectors")
    return {
        "center_cyclic_neighbor_order": neighbors,
        "method": "first-hop membership across all shortest paths from the structural center",
        "sectors": [{"via_neighbor": neighbor, "shortest_path_members": sectors[neighbor]}
                    for neighbor in neighbors],
        "vertex_sector_assignments": assignments,
        "vertices_not_honestly_assignable_to_one_sector": shared,
        "coherent_disjoint_sectors_exist": not shared,
    }


def build_blueprint(graph_path: Path, analysis_path: Path) -> dict:
    graph, graph_doc, analysis_doc = load_authoritative_graph(graph_path, analysis_path)
    center, degree, closeness = select_structural_center(graph)
    distances = nx.single_source_shortest_path_length(graph, center)
    layers = {
        str(distance): sorted(node for node, value in distances.items() if value == distance)
        for distance in range(max(distances.values()) + 1)
    }
    embedding = restore_embedding(graph, analysis_doc)
    blocks = biconnected_structure(graph, center, distances)
    faces = face_structure(embedding, center, distances, closeness)
    sectors = sector_structure(graph, embedding, center, distances, blocks)
    block_membership = {node: [] for node in graph}
    for block in blocks["blocks"]:
        for node in block["member_vertices"]:
            block_membership[node].append(block["id"])
    bridge_incident = {node for edge in blocks["bridges"] for node in edge}
    leaf_members = {node for block in blocks["blocks"] if block["leaf_block"]
                    for node in block["member_vertices"]}
    vertex_blueprint = []
    for node in sorted(graph):
        assignment = sectors["vertex_sector_assignments"][node]
        vertex_blueprint.append({
            "id": node,
            "type": graph.nodes[node]["type"],
            "structural_center": node == center,
            "bfs_layer": distances[node],
            "degree": degree[node],
            "closeness": closeness[node],
            "block_membership": block_membership[node],
            "articulation": node in blocks["articulation_vertices"],
            "incident_to_bridge": node in bridge_incident,
            "member_of_leaf_block": node in leaf_members,
            "embedding_cyclic_neighbors": list(embedding.neighbors_cw_order(node)),
            "sector_classification": assignment["classification"],
            "candidate_sectors": assignment["sectors"],
        })
    return {
        "schema": "graph-relax.topological-placement-blueprint.v1",
        "authoritative_input": graph_path.as_posix(),
        "coordinate_policy": "no coordinates assigned",
        "graph_counts": {
            "vertices": graph.number_of_nodes(), "edges": graph.number_of_edges()
        },
        "structural_center_selection": {
            "rule": ["maximum closeness centrality", "maximum degree", "lexical vertex identifier"],
            "chosen_vertex": center,
            "closeness": closeness[center],
            "degree": degree[center],
        },
        "bfs_layers": layers,
        "block_structure": blocks,
        "planar_embedding": {
            "cyclic_neighbor_order": {
                node: list(embedding.neighbors_cw_order(node)) for node in sorted(embedding)
            },
            **faces,
        },
        "embedding_sectors": sectors,
        "vertices": vertex_blueprint,
        "source_graph_schema": graph_doc["schema"],
    }


def write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(graph_path: Path, analysis_path: Path, output_path: Path) -> dict:
    document = build_blueprint(graph_path, analysis_path)
    write_json(output_path, document)
    return document


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=Path, default=Path("output/graph_first/iamp_graph.json"))
    parser.add_argument("--analysis", type=Path,
                        default=Path("output/graph_first/iamp_graph_analysis.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("output/graph_first/iamp_topology_blueprint.json"))
    args = parser.parse_args()
    document = run(args.graph, args.analysis, args.output)
    print(json.dumps({
        "structural_center": document["structural_center_selection"],
        "bfs_layers": document["bfs_layers"],
        "block_count": len(document["block_structure"]["blocks"]),
        "face_count": len(document["planar_embedding"]["faces"]),
        "best_outer_face_candidates": document["planar_embedding"]["best_outer_face_candidates"],
        "center_cyclic_neighbor_order": document["embedding_sectors"]["center_cyclic_neighbor_order"],
        "shared_sector_vertices": document["embedding_sectors"]["vertices_not_honestly_assignable_to_one_sector"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
