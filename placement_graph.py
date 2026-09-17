"""Separate electrical incidence and source-derived device placement views.

SIGNAL means unordered co-membership in a non-rail net, never a physical wire.
Only explicit SOURCE_FACT SOFT_NEARNESS records add non-electrical adjacency.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations
from urllib.parse import quote

import networkx as nx


@dataclass(frozen=True)
class PlacementGraph:
    topology: nx.Graph
    rails: tuple[str, ...]
    rail_classification: dict
    blocks: tuple[dict, ...]
    groups: tuple[dict, ...]
    rail_attachments: tuple[dict, ...]


def classify_rails(circuit):
    """Consume upstream rail semantics; never discover rails from placement data."""
    metadata = circuit.source_metadata or {}
    semantics = metadata.get('semantics', {})
    record = semantics.get('rail_classification')

    # Backward-compatible handling for older .circuit files that predate the
    # semantic contract.  This path still uses only explicit source metadata.
    if record is None:
        normalized = metadata.get('normalized', {})
        rails = normalized.get('rails', [])
        if not isinstance(rails, (list, tuple)) or any(not isinstance(x, str) for x in rails):
            raise ValueError('source_metadata.normalized.rails must be a list of net names')
        if len(set(rails)) != len(rails):
            raise ValueError('duplicate rail classification')
        unknown = set(rails) - {n.name for n in circuit.nets}
        if unknown:
            raise ValueError(f'rail classification references unknown nets: {sorted(unknown)}')
        provenance = {
            'state': 'EXPLICIT' if 'rails' in normalized else 'UNCLASSIFIED',
            'status': 'EXPLICIT_SOURCE_METADATA' if 'rails' in normalized else 'UNCLASSIFIED',
            'pointer': '/source_metadata/normalized/rails' if 'rails' in normalized else None,
            'source': metadata.get('source'),
            'sha256': metadata.get('sha256'),
            'authority': deepcopy(normalized.get('provenance', {})),
            'policy': 'No name, voltage, degree or electronics-convention inference.',
        }
        return tuple(sorted(rails)), provenance

    if not isinstance(record, dict):
        raise ValueError('source_metadata.semantics.rail_classification must be an object')

    state = record.get('state')
    if state not in {'EXPLICIT', 'EVIDENCE_DERIVED', 'UNCLASSIFIED'}:
        raise ValueError(f'invalid rail classification state: {state!r}')

    rails = record.get('rails')
    if not isinstance(rails, (list, tuple)) or any(not isinstance(x, str) for x in rails):
        raise ValueError('rail classification rails must be a list of net names')
    if len(set(rails)) != len(rails):
        raise ValueError('duplicate rail classification')
    if state == 'UNCLASSIFIED' and rails:
        raise ValueError('UNCLASSIFIED rail semantics cannot contain rails')

    unknown = set(rails) - {n.name for n in circuit.nets}
    if unknown:
        raise ValueError(f'rail classification references unknown nets: {sorted(unknown)}')

    provenance = deepcopy(record)
    provenance['status'] = state
    provenance['source'] = metadata.get('source')
    provenance['sha256'] = metadata.get('sha256')
    provenance['pointer'] = '/source_metadata/semantics/rail_classification'

    return tuple(sorted(rails)), provenance


def build_electrical_graph(circuit):
    """Complete authoritative component/net multigraph, one edge per terminal.

    NC declarations have no net edge; retain them and the complete pin inventory
    as graph metadata. Parallel pins on the same component/net are never collapsed.
    """
    graph = nx.MultiGraph(view='ELECTRICAL')
    graph.graph['nc'] = tuple(circuit.nc)
    graph.graph['nets'] = tuple((n.name, tuple(n.terminals)) for n in circuit.nets)
    for ref in circuit.refs:
        graph.add_node('component:' + ref, kind='COMPONENT', ref=ref,
                       pins=tuple(circuit.pins[ref]), device=circuit.devices.get(ref))
    for net in circuit.nets:
        graph.add_node('net:' + net.name, kind='NET', name=net.name)
        for terminal in net.terminals:
            ref, pin = terminal.rsplit('.', 1)
            graph.add_edge('component:' + ref, 'net:' + net.name, key=terminal,
                           terminal=terminal, pin=pin, net=net.name, kind='ELECTRICAL_INCIDENCE')
    return nx.freeze(graph)


def _add_relation(graph, first, second, record):
    if not graph.has_edge(first, second):
        graph.add_edge(first, second, relationships=[])
    graph[first][second]['relationships'].append(record)


def power_port_id(host_ref, host_pin, electrical_net):
    """Escape every identity segment so separators cannot collide with source names."""
    return 'POWER_PORT::' + '::'.join(quote(x, safe='') for x in (host_ref, host_pin, electrical_net))


def device_topology(model):
    """Physical-device projection; POWER_ATTACHMENT never creates device adjacency."""
    return model.topology.subgraph(n for n, data in model.topology.nodes(data=True)
                                  if data['object_type'] == 'DEVICE')


def build_placement_graph(circuit):
    """Physical-device neighborhoods plus one local leaf per rail terminal."""
    rails, classification = classify_rails(circuit)
    graph = nx.Graph(view='PLACEMENT')
    graph.add_nodes_from((ref, {'object_type': 'DEVICE', 'ref': ref}) for ref in sorted(circuit.refs))
    metadata = circuit.source_metadata or {}
    parts = metadata.get('normalized', {}).get('parts', [])
    source_parts = {p['ref']: (i, p) for i, p in enumerate(parts)}
    signals = {}
    attachments = []
    for net in sorted(circuit.nets, key=lambda n: n.name):
        if net.name in rails:
            for terminal in sorted(net.terminals):
                ref, pin = terminal.rsplit('.', 1)
                identity = power_port_id(ref, pin, net.name)
                if identity in graph:
                    raise ValueError(f'duplicate POWER_PORT identity {identity}')
                provenance = {'kind': 'ELECTRICAL_INCIDENCE_PROJECTION',
                              'terminal': terminal, 'electrical_net': net.name,
                              'source': metadata.get('source'), 'sha256': metadata.get('sha256'),
                              'authority': deepcopy(metadata.get('normalized', {}).get('provenance', {}))}
                if ref in source_parts:
                    index, part = source_parts[ref]
                    pointer_pin = pin.replace('~', '~0').replace('/', '~1')
                    provenance['pointer'] = f'/parts/{index}/terminals/{pointer_pin}/net'
                    provenance['component'] = deepcopy(part.get('provenance', {}))
                graph.add_node(identity, object_type='POWER_PORT', host_ref=ref, host_pin=pin,
                               electrical_net=net.name, provenance=provenance)
                _add_relation(graph, ref, identity, {
                    'kind': 'POWER_ATTACHMENT', 'host_ref': ref, 'host_pin': pin,
                    'electrical_net': net.name, 'provenance': deepcopy(provenance),
                    'semantics': 'LOCAL_PLACEMENT_LEAF_NOT_PHYSICAL_DEVICE_ADJACENCY'})
                attachments.append({'terminal': terminal, 'net': net.name, 'power_port': identity})
            continue
        members = sorted({t.rsplit('.', 1)[0] for t in net.terminals})
        for first, second in combinations(members, 2):
            signals.setdefault((first, second), []).append({
                'net': net.name,
                'terminals': sorted(t for t in net.terminals if t.rsplit('.', 1)[0] in (first, second)),
            })
    for (first, second), evidence in sorted(signals.items()):
        _add_relation(graph, first, second, {
            'kind': 'SIGNAL', 'semantics': 'UNORDERED_SHARED_NON_RAIL_NET',
            'evidence': evidence,
        })
    for record in sorted(circuit.soft_nearness, key=lambda r: (r['ref'], r['near'])):
        first, second = record['ref'], record['near']
        if first not in circuit.refs or second not in circuit.refs or first == second:
            raise ValueError('invalid SOFT_NEARNESS endpoints')
        if record.get('kind') != 'SOFT_NEARNESS' or record.get('provenance', {}).get('kind') != 'SOURCE_FACT':
            raise ValueError('SOFT_NEARNESS must be an explicit SOURCE_FACT')
        _add_relation(graph, first, second, deepcopy(record))
    normalized = (circuit.source_metadata or {}).get('normalized', {})
    return PlacementGraph(nx.freeze(graph), rails, classification,
                          tuple(deepcopy(normalized.get('blocks', []))),
                          tuple(deepcopy(normalized.get('groups', []))), tuple(attachments))


def analyze_placement(model: PlacementGraph):
    """Never accept the electrical graph as a placement metric input."""
    if not isinstance(model, PlacementGraph) or model.topology.graph.get('view') != 'PLACEMENT':
        raise TypeError('placement analysis requires a PlacementGraph')
    graph = device_topology(model)
    neighbors = {ref: sorted(graph.neighbors(ref)) for ref in sorted(graph)}
    overlaps = []
    for first, second in combinations(sorted(graph), 2):
        a, b = set(neighbors[first]), set(neighbors[second])
        common = sorted(a & b)
        if common:
            overlaps.append({'first': first, 'second': second, 'common_neighbors': common,
                             'count': len(common), 'jaccard': len(common) / len(a | b)})
    return {
        'view': 'PLACEMENT', 'metric_scope': 'PHYSICAL_DEVICES_ONLY',
        'degree': dict(sorted(graph.degree())),
        'device_only_degree': dict(sorted(graph.degree())),
        'total_device_degree': {r: model.topology.degree(r) for r in sorted(graph)},
        'total_placement_degree': dict(sorted(model.topology.degree())),
        'degree_centrality': nx.degree_centrality(graph),
        'closeness_centrality': nx.closeness_centrality(graph),
        'neighbors': neighbors, 'positive_neighborhood_overlaps': overlaps,
        'connected_components': sorted([sorted(c) for c in nx.connected_components(graph)]),
    }


def placement_document(model):
    graph = model.topology
    return {
        'schema': 'graph-relax.placement-graph.v2', 'view': 'PLACEMENT',
        'semantics': 'Non-rail device co-membership, explicit soft nearness and local POWER_PORT leaves; no physical wire topology.',
        'vertices': sorted(graph),
        'objects': [{'id': n, **deepcopy(graph.nodes[n]), 'placement_degree': graph.degree(n)} for n in sorted(graph)],
        'edges': [{'source': a, 'target': b, 'relationships': deepcopy(graph[a][b]['relationships'])}
                  for a, b in sorted(tuple(sorted(pair)) for pair in graph.edges())],
        'rails': list(model.rails), 'rail_classification': deepcopy(model.rail_classification),
        'blocks': deepcopy(list(model.blocks)), 'groups': deepcopy(list(model.groups)),
        'rail_incidences_projected_to_local_ports': list(model.rail_attachments),
        'analysis': analyze_placement(model),
    }


def electrical_document(graph):
    if graph.graph.get('view') != 'ELECTRICAL':
        raise TypeError('expected electrical graph')
    return {
        'vertices': [{'id': v, **deepcopy(graph.nodes[v])} for v in sorted(graph)],
        'incidences': sorted([{'component': a if a.startswith('component:') else b,
                              'net_vertex': b if b.startswith('net:') else a, **deepcopy(data)}
                             for a, b, _, data in graph.edges(keys=True, data=True)],
                            key=lambda x: x['terminal']),
        'nets': deepcopy(graph.graph['nets']), 'nc': graph.graph['nc'],
    }


def comparison_statistics(circuit, placement):
    """Diagnostic baseline only: all shared nets, with identical soft records.

    The unfiltered graph is never exposed as an initial-placement input and rail
    relationships in this baseline are never mislabeled SIGNAL.
    """
    before = nx.Graph()
    before.add_nodes_from(circuit.refs)
    rail_pairs = set()
    for net in circuit.nets:
        pairs = list(combinations(sorted({t.rsplit('.', 1)[0] for t in net.terminals}), 2))
        before.add_edges_from(pairs)
        if net.name in placement.rails:
            rail_pairs.update(pairs)
    before.add_edges_from((r['ref'], r['near']) for r in circuit.soft_nearness)
    after = device_topology(placement)
    total = placement.topology
    ports = [n for n, d in total.nodes(data=True) if d["object_type"] == "POWER_PORT"]
    electrical = build_electrical_graph(circuit)
    def top(graph):
        return [{'ref': r, 'degree': d} for r, d in sorted(graph.degree(), key=lambda x: (-x[1], x[0]))[:10]]
    def counts(kind):
        return sum(r['kind'] == kind for _, _, edge in total.edges(data=True) for r in edge['relationships'])
    return {
        'schema': 'graph-relax.placement-graph-statistics.v2',
        'baseline': 'All unordered shared-net device pairs plus the same explicit SOFT_NEARNESS; diagnostic only.',
        'electrical_vertex_count': electrical.number_of_nodes(),
        'electrical_incidence_count': electrical.number_of_edges(),
        'electrical_net_count': len(circuit.nets), 'electrical_nc_count': len(circuit.nc),
        'placement_vertex_count': total.number_of_nodes(), 'placement_edge_count': total.number_of_edges(),
        'physical_device_count': after.number_of_nodes(), 'power_port_count': len(ports),
        'power_attachment_relationship_count': counts('POWER_ATTACHMENT'),
        'no_global_rail_placement_hubs': all(n not in total and 'net:' + n not in total for n in placement.rails)
            and all(total.degree(p) == 1 and set(total.neighbors(p)) == {total.nodes[p]['host_ref']} for p in ports),
        'maximum_placement_degree_before_global_rail_suppression': max(dict(before.degree()).values(), default=0),
        'maximum_placement_degree_before_local_power_port_conversion': max(dict(after.degree()).values(), default=0),
        'maximum_device_only_degree_after_conversion': max(dict(after.degree()).values(), default=0),
        'maximum_total_device_degree_after_conversion': max((total.degree(n) for n in after), default=0),
        'maximum_power_port_degree': max((total.degree(n) for n in ports), default=0),
        'local_power_port_conversion': {
            'before': {'view': 'PREVIOUS_RAIL_FILTERED_DEVICE_ONLY',
                       'vertices': after.number_of_nodes(), 'edges': after.number_of_edges(),
                       'maximum_degree': max(dict(after.degree()).values(), default=0)},
            'after': {'view': 'DEVICES_AND_LOCAL_POWER_PORTS',
                      'vertices': total.number_of_nodes(), 'edges': total.number_of_edges(),
                      'device_only_degree': dict(sorted(after.degree())),
                      'total_device_degree': {n: total.degree(n) for n in sorted(after)}}},
        'rail_incidences_excluded_from_device_adjacency': len(placement.rail_attachments),
        'rail_generated_pair_count_before': len(rail_pairs),
        'signal_relationship_count': counts('SIGNAL'), 'soft_nearness_relationship_count': counts('SOFT_NEARNESS'),
        'before_suppression': {'vertex_count': before.number_of_nodes(), 'edge_count': before.number_of_edges(),
                               'highest_degree_devices': top(before), 'degree': dict(sorted(before.degree()))},
        'after_suppression': {'vertex_count': after.number_of_nodes(), 'edge_count': after.number_of_edges(),
                              'highest_degree_devices': top(after), 'degree': dict(sorted(after.degree())),
                              'isolated_devices': sorted(nx.isolates(after))},
    }
