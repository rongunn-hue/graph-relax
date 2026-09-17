"""Versioned B/C interfaces and geometry-free semantic validation.

This module does not import a producer, circuit parser, or placement algorithm.
"""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import argparse
import json
import math
from pathlib import Path
import re
from types import SimpleNamespace

SEMANTIC_SCHEMA = 'graph-relax.semantic-placement.v1'
INITIAL_SCHEMA = 'graph-relax.initial-placement.v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical_bytes(document):
    return (json.dumps(document, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def digest(document):
    return sha256(canonical_bytes(document)).hexdigest()


def load_document(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f'duplicate JSON key: {key}')
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def text(value):
    require(isinstance(value, str) and bool(value), 'expected nonempty string')
    return value


def records(value):
    require(isinstance(value, list), 'expected array')
    return value


def unique_strings(value):
    values = [text(v) for v in records(value)]
    require(len(values) == len(set(values)), 'duplicate identity')
    return set(values)


def validate_soft_nearness(record):
    require(record.get('kind') == 'SOFT_NEARNESS', 'invalid nearness kind')
    for key in ('ref', 'near'):
        require(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', text(record.get(key))) is not None,
                'invalid nearness reference')
    require(record['ref'] != record['near'], 'self nearness')
    if 'purpose' in record:
        require(isinstance(record['purpose'], str), 'invalid purpose')
    p = record.get('provenance', {})
    require(p.get('kind') == 'SOURCE_FACT', 'nearness requires SOURCE_FACT')
    for key in ('source', 'pointer'):
        text(p.get(key))
    require(re.fullmatch('[0-9a-f]{64}', text(p.get('sha256'))) is not None, 'invalid source digest')


def _semantic(doc):
    require(doc['schema'] == SEMANTIC_SCHEMA, 'wrong semantic schema')
    source, e, p = doc['source'], doc['electrical'], doc['placement']
    text(source['path'])
    require(re.fullmatch('[0-9a-f]{64}', text(source['sha256'])) is not None, 'invalid input hash')
    require(isinstance(source['metadata'], dict), 'invalid source metadata')
    components, pins = {}, set()
    for c in records(e['components']):
        ref = text(c['ref'])
        require(ref not in components, 'duplicate component')
        require(c['device'] is None or isinstance(c['device'], str), 'invalid device type')
        components[ref] = c
        for pin in unique_strings(c['pins']):
            terminal = ref + '.' + pin
            require(terminal not in pins, 'ambiguous terminal identity')
            pins.add(terminal)
    nets, assignments = {}, {}
    for n in records(e['nets']):
        name = text(n['name'])
        require(name not in nets, 'duplicate net')
        nets[name] = unique_strings(n['terminals'])
        for terminal in nets[name]:
            require(terminal in pins and terminal not in assignments, 'unknown/conflicting terminal')
            assignments[terminal] = name
    nc = unique_strings(e['nc'])
    require(nc <= pins and not nc.intersection(assignments), 'invalid NC reference/disposition')
    require(set(assignments) | nc == pins, 'lost terminal disposition')
    incidences = {}
    for inc in records(e['incidences']):
        terminal = text(inc['terminal'])
        require(terminal not in incidences, 'duplicate electrical incidence')
        require(inc['ref'] in components and inc['pin'] in components[inc['ref']]['pins'], 'unknown pin')
        require(terminal == inc['ref'] + '.' + inc['pin'], 'inconsistent terminal identity')
        require(terminal in assignments and assignments[terminal] == inc['net'], 'invalid incidence')
        incidences[terminal] = inc['net']
    require(incidences == assignments, 'lost electrical incidence')
    body = {k: v for k, v in e.items() if k != 'sha256'}
    require(e['sha256'] == digest(body), 'electrical identity mismatch')
    text(e['name'])
    rails = unique_strings(p['rails'])
    require(rails <= set(nets), 'unknown rail')
    classification = p['rail_classification']
    state = classification['state']
    for key in ('source', 'sha256', 'pointer', 'status'):
        require(key in classification, 'missing rail provenance/status')
    require(classification['source'] == source['metadata'].get('source') and
            classification['sha256'] == source['metadata'].get('sha256'), 'rail source identity mismatch')
    require(state in ('EXPLICIT', 'EVIDENCE_DERIVED', 'UNCLASSIFIED'), 'invalid rail state')
    require(state != 'UNCLASSIFIED' or not rails, 'unclassified rails must be empty')
    if 'rails' in classification:
        require(unique_strings(classification['rails']) == rails, 'inconsistent rail list')
    # Blocks/groups remain opaque source metadata, not graph edges.
    normalized = source['metadata'].get('normalized', {})
    require(records(p['blocks']) == normalized.get('blocks', []), 'source blocks changed')
    require(records(p['groups']) == normalized.get('groups', []), 'source groups changed')
    upstream = source['metadata'].get('semantics', {}).get('rail_classification')
    if upstream is not None:
        require(upstream['state'] == state and set(upstream['rails']) == rails,
                'classification differs from upstream')
        require(all(classification.get(k) == v for k, v in upstream.items()), 'rail provenance changed')
    else:
        require(state == ('EXPLICIT' if 'rails' in normalized else 'UNCLASSIFIED'), 'invented rail state')
        require(set(normalized.get('rails', [])) == rails, 'invented rails')
    objects, ports = {}, {}
    for obj in records(p['objects']):
        identity = text(obj['id'])
        require(identity not in objects, 'duplicate placement object')
        objects[identity] = obj
        if obj['object_type'] == 'DEVICE':
            require(identity in components and obj['ref'] == identity, 'unknown DEVICE')
        else:
            require(obj['object_type'] == 'POWER_PORT', 'invalid object type')
            host, pin, net = obj['host_ref'], obj['host_pin'], obj['electrical_net']
            terminal = host + '.' + pin
            require(host in components and pin in components[host]['pins'], 'unknown port host/pin')
            require(net in rails and assignments.get(terminal) == net, 'port without rail incidence')
            require(terminal not in ports, 'duplicate rail projection')
            require(isinstance(obj['provenance'], dict), 'missing port provenance')
            ports[terminal] = identity
    require({i for i, o in objects.items() if o['object_type'] == 'DEVICE'} == set(components), 'DEVICE coverage')
    require(set(ports) == {t for t, n in assignments.items() if n in rails}, 'rail/port bijection')
    expected_signals = {}
    for net, terminals in nets.items():
        if net in rails:
            continue
        members = {t.rsplit('.', 1)[0] for t in terminals}
        for pair in combinations(sorted(members), 2):
            expected_signals.setdefault(pair, {})[net] = sorted(t for t in terminals if t.rsplit('.', 1)[0] in pair)
    signals, soft, attachments, edge_pairs = {}, [], set(), set()
    degrees = dict.fromkeys(objects, 0)
    for edge in records(p['edges']):
        a, b = edge['source'], edge['target']
        require(a in objects and b in objects and a != b, 'invalid edge endpoints')
        pair = tuple(sorted((a, b)))
        require(pair not in edge_pairs, 'duplicate placement edge')
        edge_pairs.add(pair)
        degrees[a] += 1
        degrees[b] += 1
        rels = records(edge['relationships'])
        require(bool(rels), 'empty relationship edge')
        seen = set()
        for r in rels:
            kind = r['kind']
            key = (kind, r.get('ref'), r.get('near'))
            require(key not in seen, 'duplicate placement relationship')
            seen.add(key)
            if kind == 'SIGNAL':
                require(a in components and b in components, 'SIGNAL endpoints must be DEVICE')
                require(r['semantics'] == 'UNORDERED_SHARED_NON_RAIL_NET', 'invalid SIGNAL semantics')
                evidence = {}
                for item in records(r['evidence']):
                    net = item['net']
                    require(net not in evidence, 'duplicate SIGNAL evidence')
                    evidence[net] = sorted(unique_strings(item['terminals']))
                signals[pair] = evidence
            elif kind == 'SOFT_NEARNESS':
                validate_soft_nearness(r)
                require(a in components and b in components and {r['ref'], r['near']} == {a, b}, 'invalid nearness endpoints')
                soft.append(r)
            else:
                require(kind == 'POWER_ATTACHMENT', 'unknown relationship kind')
                terminal = r['host_ref'] + '.' + r['host_pin']
                port = ports.get(terminal)
                require(port is not None and {a, b} == {port, r['host_ref']}, 'invalid power attachment')
                require(r['electrical_net'] == assignments[terminal], 'attachment net mismatch')
                require(terminal not in attachments, 'duplicate attachment')
                require(r['provenance'] == objects[port]['provenance'], 'attachment provenance mismatch')
                require(r['semantics'] == 'LOCAL_PLACEMENT_LEAF_NOT_PHYSICAL_DEVICE_ADJACENCY', 'invalid attachment semantics')
                attachments.add(terminal)
    require(signals == expected_signals, 'lost/invented SIGNAL evidence')
    require(attachments == set(ports) and all(degrees[p] == 1 for p in ports.values()), 'port degree/bijection failure')
    expected_soft = records(source['soft_nearness'])
    for r in expected_soft:
        validate_soft_nearness(r)
    require(sorted(map(digest, soft)) == sorted(map(digest, expected_soft)), 'source nearness changed')
    require(len({(r['ref'], r['near']) for r in soft}) == len(soft), 'duplicate nearness')
    return doc


def validate_semantic(doc):
    try:
        return _semantic(doc)
    except (KeyError, TypeError, AttributeError, IndexError) as exc:
        raise ValueError(f'malformed semantic document: {exc}') from exc


def semantic_views(doc):
    """Deserialize neutral containers; never invoke Stage B or classify rails."""
    validate_semantic(doc)
    import networkx as nx
    e, p = doc['electrical'], doc['placement']
    circuit = SimpleNamespace(name=e['name'], refs=tuple(c['ref'] for c in e['components']),
        pins={c['ref']: tuple(c['pins']) for c in e['components']},
        devices={c['ref']: c['device'] for c in e['components']},
        nets=tuple(SimpleNamespace(name=n['name'], terminals=tuple(n['terminals'])) for n in e['nets']),
        nc=tuple(e['nc']), soft_nearness=deepcopy(doc['source']['soft_nearness']))
    graph = nx.Graph(view='PLACEMENT')
    for obj in p['objects']:
        graph.add_node(obj['id'], **{k: deepcopy(v) for k, v in obj.items() if k != 'id'})
    for edge in p['edges']:
        graph.add_edge(edge['source'], edge['target'], relationships=deepcopy(edge['relationships']))
    placement = SimpleNamespace(topology=nx.freeze(graph), rails=tuple(p['rails']),
        rail_classification=deepcopy(p['rail_classification']), blocks=deepcopy(p['blocks']), groups=deepcopy(p['groups']))
    return circuit, placement


def finite(value, positive=False):
    require(type(value) in (float, int) and math.isfinite(value), 'nonfinite/non-numeric coordinate or dimension')
    if positive:
        require(value > 0, 'nonpositive dimension')
    return value


def _initial(doc, semantic):
    validate_semantic(semantic)
    require(doc['schema'] == INITIAL_SCHEMA, 'wrong initial-placement schema')
    require(doc['input'] == {'schema': SEMANTIC_SCHEMA, 'sha256': digest(semantic)}, 'Stage B identity mismatch')
    refs = {c['ref'] for c in semantic['electrical']['components']}
    method = doc['method']
    require(method['hierarchy'] == 'SOURCE_BLOCKS_ELSE_NONRAIL_BLOCK_CUT' and
            method['geometry'] == 'ONE_DIRECT_DENSITY_AWARE_TIDY_TREE_CONSTRUCTION', 'invalid method')
    require(all(method[k] is False for k in ('search', 'optimization', 'relaxation', 'candidate_family', 'user_hints')),
            'invalid construction method flags')
    require(doc['hierarchy_source'] in ('SOURCE_BLOCKS', 'DERIVED_NONRAIL_BLOCK_CUT'), 'invalid hierarchy source')
    require(doc['hierarchy_source'] == ('SOURCE_BLOCKS' if semantic['placement']['blocks'] else 'DERIVED_NONRAIL_BLOCK_CUT'), 'inconsistent hierarchy source')
    w = finite(doc['device_envelope']['width'], True)
    h = finite(doc['device_envelope']['height'], True)
    require(finite(doc['clearance']) >= 0, 'negative clearance')
    require(finite(doc['territory_gap']) >= 0, 'negative territory gap')
    coordinates = doc['coordinates']
    require(isinstance(coordinates, dict) and set(coordinates) == refs, 'coordinate coverage')
    for p in coordinates.values():
        finite(p['x']); finite(p['y'])
    for a, b in combinations(coordinates.values(), 2):
        require(abs(a['x'] - b['x']) >= w or abs(a['y'] - b['y']) >= h, 'device envelope overlap')
    territories, members = {}, set()
    kinds = {'SOURCE_BLOCK', 'SOURCE_UNASSIGNED', 'POWER_ONLY', 'TREE_COMPONENT', 'CORE', 'BRANCH'}
    for t in records(doc['territories']):
        identity = text(t['id'])
        require(identity not in territories and t['kind'] in kinds, 'invalid territory identity/kind')
        territories[identity] = t
        require(t['source_name'] is None or isinstance(t['source_name'], str), 'invalid source name')
        group = unique_strings(t['devices'])
        require(group and group <= refs and not group.intersection(members), 'invalid/excess territory ownership')
        members.update(group)
        box = t['box']
        x, y = finite(box['x']), finite(box['y'])
        tw, th = finite(box['width'], True), finite(box['height'], True)
        for ref in group:
            p = coordinates[ref]
            require(x <= p['x']-w/2 and p['x']+w/2 <= x+tw and y <= p['y']-h/2 and p['y']+h/2 <= y+th, 'device outside territory')
    require(members == refs, 'territory coverage')
    edges = set()
    for edge in records(doc['territory_edges']):
        require(isinstance(edge, list) and len(edge) == 2, 'invalid territory relationship')
        a, b = edge
        require(a in territories and b in territories and a != b, 'unknown territory relationship')
        pair = tuple(sorted(edge))
        require(pair not in edges, 'duplicate territory relationship')
        edges.add(pair)
    return doc


def validate_initial(doc, semantic):
    try:
        return _initial(doc, semantic)
    except (KeyError, TypeError, AttributeError, IndexError) as exc:
        raise ValueError(f'malformed initial-placement document: {exc}') from exc


def validate_topology(doc, semantic, initial):
    """Stage D validator entry point; admitted B/C contracts are unchanged."""
    from topology_contracts import validate_topology as validate
    return validate(doc, semantic, initial)


def main():
    parser = argparse.ArgumentParser(description='Validate a versioned B, C or D document without running a stage.')
    parser.add_argument('document', type=Path)
    parser.add_argument('--semantic', type=Path, help='required Stage B document for Stage C validation')
    parser.add_argument('--initial', type=Path, help='required Stage C document for Stage D validation')
    args = parser.parse_args()
    try:
        doc = load_document(args.document)
        if args.initial:
            require(args.semantic is not None, 'Stage D validation requires --semantic and --initial')
            validate_topology(doc, load_document(args.semantic), load_document(args.initial))
        elif args.semantic:
            validate_initial(doc, load_document(args.semantic))
        else:
            validate_semantic(doc)
    except (ValueError, OSError) as exc:
        parser.exit(2, f'{exc}\n')
    print('VALID')


if __name__ == '__main__':
    main()
