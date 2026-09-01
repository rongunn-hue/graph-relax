#!/usr/bin/env python3
"""Small deterministic physical-relaxation experiment for .circuit graphs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import shlex
from dataclasses import dataclass
from pathlib import Path


# These constants define the approved first experiment. They are intentionally
# kept here, rather than exposed as a growing configuration framework.
WIDTH = 60.0
HEIGHT = 40.0
BODY_CLEARANCE = 12.0
NET_BODY_CLEARANCE = 10.0
NET_NET_CLEARANCE = 8.0
INITIAL_STEP = 16.0
MIN_STEP = 0.25
MAX_SWEEPS = 500
SEED = 1729
JITTER = 5.0
ROUND_DIGITS = 6
WEIGHTS = {
    "length": 0.02,
    "overlap": 100000.0,
    "clearance": 5000.0,
    "crossing": 2000.0,
    "crowding": 20.0,
}
MOVES = ((1, 0), (-1, 0), (0, 1), (0, -1),
         (1, 1), (1, -1), (-1, 1), (-1, -1))

# Frozen Experiment 2 search parameters. The geometry and objective above are
# shared unchanged with Experiment 1.
ANNEAL_SEED = 314159
ANNEAL_ITERATIONS = 6000
INITIAL_TEMPERATURE = 25000.0
FINAL_TEMPERATURE = 10.0
INITIAL_MOVE_DISTANCE = 140.0
FINAL_MOVE_DISTANCE = 0.5
DEGREE_RING_SPACING = math.hypot(WIDTH, HEIGHT) + BODY_CLEARANCE


@dataclass(frozen=True)
class Net:
    name: str
    terminals: tuple[str, ...]


@dataclass(frozen=True)
class Circuit:
    name: str
    refs: tuple[str, ...]
    pins: dict[str, tuple[str, ...]]
    nets: tuple[Net, ...]


def natural_key(value: str) -> tuple:
    return tuple(int(p) if p.isdigit() else p for p in re.split(r"(\d+)", value))


def parse_circuit(path: Path) -> Circuit:
    title = ""
    refs = []
    nets = []
    pins: dict[str, set[str]] = {}
    used_terminals = set()
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        words = shlex.split(raw, comments=True)
        if not words:
            continue
        if words[0] == "circuit" and len(words) == 2:
            title = words[1]
        elif words[0] == "component" and len(words) == 3:
            ref = words[1]
            if ref in pins:
                raise ValueError(f"line {line_number}: duplicate component {ref}")
            refs.append(ref)
            pins[ref] = set()
        elif words[0] == "net" and len(words) >= 3:
            name = words[1]
            terminals = tuple(words[2:])
            if any(n.name == name for n in nets):
                raise ValueError(f"line {line_number}: duplicate net {name}")
            for terminal in terminals:
                if "." not in terminal:
                    raise ValueError(f"line {line_number}: malformed terminal {terminal}")
                ref, pin = terminal.rsplit(".", 1)
                if ref not in pins:
                    raise ValueError(f"line {line_number}: unknown component {ref}")
                if terminal in used_terminals:
                    raise ValueError(f"line {line_number}: terminal used twice: {terminal}")
                used_terminals.add(terminal)
                pins[ref].add(pin)
            nets.append(Net(name, terminals))
        else:
            raise ValueError(f"line {line_number}: unsupported statement")
    if not title:
        raise ValueError("missing circuit declaration")
    return Circuit(
        title,
        tuple(sorted(refs, key=natural_key)),
        {r: tuple(sorted(pins[r], key=natural_key)) for r in refs},
        tuple(nets),
    )


def initial_positions(circuit: Circuit) -> dict[str, tuple[float, float]]:
    columns = math.ceil(math.sqrt(len(circuit.refs)))
    rng = random.Random(SEED)
    positions = {}
    for index, ref in enumerate(circuit.refs):
        row, column = divmod(index, columns)
        positions[ref] = (
            column * 100.0 + rng.uniform(-JITTER, JITTER),
            row * 80.0 + rng.uniform(-JITTER, JITTER),
        )
    return positions


def component_degrees(circuit: Circuit) -> dict[str, int]:
    connected = {ref: set() for ref in circuit.refs}
    for net in circuit.nets:
        if len(net.terminals) <= 1:
            continue
        for terminal in net.terminals:
            connected[terminal.rsplit(".", 1)[0]].add(net.name)
    return {ref: len(connected[ref]) for ref in circuit.refs}


def degree_seeded_positions(circuit: Circuit) -> dict[str, tuple[float, float]]:
    degrees = component_degrees(circuit)
    bands = {}
    for ref, degree in degrees.items():
        bands.setdefault(degree, []).append(ref)
    ordered_degrees = sorted(bands, reverse=True)
    positions = {}
    previous_radius = 0.0
    first_degree = ordered_degrees[0]
    first_refs = sorted(bands[first_degree])
    if len(first_refs) == 1:
        positions[first_refs[0]] = (0.0, 0.0)
        remaining_degrees = ordered_degrees[1:]
    else:
        remaining_degrees = ordered_degrees
        previous_radius = -DEGREE_RING_SPACING
    for degree in remaining_degrees:
        refs = sorted(bands[degree])
        count = len(refs)
        chord_radius = (DEGREE_RING_SPACING / (2.0 * math.sin(math.pi / count))
                        if count > 1 else 0.0)
        radius = max(previous_radius + DEGREE_RING_SPACING, chord_radius)
        for index, ref in enumerate(refs):
            angle = -math.pi / 2.0 + 2.0 * math.pi * index / count
            positions[ref] = (radius * math.cos(angle), radius * math.sin(angle))
        previous_radius = radius
    return positions


def local_terminal(index: int, count: int) -> tuple[float, float]:
    # Equally spaced perimeter samples, offset half a slot from the top-left.
    perimeter = 2.0 * (WIDTH + HEIGHT)
    distance = (index + 0.5) * perimeter / count
    if distance < WIDTH:
        return (-WIDTH / 2.0 + distance, -HEIGHT / 2.0)
    distance -= WIDTH
    if distance < HEIGHT:
        return (WIDTH / 2.0, -HEIGHT / 2.0 + distance)
    distance -= HEIGHT
    if distance < WIDTH:
        return (WIDTH / 2.0 - distance, HEIGHT / 2.0)
    distance -= WIDTH
    return (-WIDTH / 2.0, HEIGHT / 2.0 - distance)


def terminal_positions(circuit: Circuit, positions: dict[str, tuple[float, float]]):
    result = {}
    for ref in circuit.refs:
        pins = circuit.pins[ref]
        for index, pin in enumerate(pins):
            dx, dy = local_terminal(index, len(pins))
            x, y = positions[ref]
            result[f"{ref}.{pin}"] = (x + dx, y + dy)
    return result


def net_segments(circuit: Circuit, terminals):
    segments = []
    for net in circuit.nets:
        points = [terminals[t] for t in net.terminals]
        if len(points) <= 1:
            continue
        centroid = (sum(p[0] for p in points) / len(points),
                    sum(p[1] for p in points) / len(points))
        for terminal, point in zip(net.terminals, points):
            segments.append((net.name, terminal, point, centroid))
    return segments


def initial_junctions(circuit: Circuit, positions):
    terminals = terminal_positions(circuit, positions)
    junctions = {}
    for net in circuit.nets:
        if len(net.terminals) >= 3:
            points = [terminals[t] for t in net.terminals]
            junctions[net.name] = (
                sum(p[0] for p in points) / len(points),
                sum(p[1] for p in points) / len(points),
            )
    return junctions


def junction_segments(circuit: Circuit, terminals, junctions):
    """Actual Experiment 3 segments, including endpoint component owners."""
    segments = []
    for net in circuit.nets:
        if len(net.terminals) == 2:
            first, second = net.terminals
            owners = (first.rsplit(".", 1)[0], second.rsplit(".", 1)[0])
            segments.append((net.name, first, terminals[first], terminals[second], owners))
        elif len(net.terminals) >= 3:
            junction = junctions[net.name]
            for terminal in net.terminals:
                owner = terminal.rsplit(".", 1)[0]
                segments.append((net.name, terminal, terminals[terminal], junction, (owner,)))
    return segments


def rect_bounds(point):
    x, y = point
    return (x - WIDTH / 2, y - HEIGHT / 2, x + WIDTH / 2, y + HEIGHT / 2)


def rect_gap(a, b):
    ax1, ay1, ax2, ay2 = rect_bounds(a)
    bx1, by1, bx2, by2 = rect_bounds(b)
    dx = max(bx1 - ax2, ax1 - bx2, 0.0)
    dy = max(by1 - ay2, ay1 - by2, 0.0)
    return math.hypot(dx, dy)


def overlap_area(a, b):
    ax1, ay1, ax2, ay2 = rect_bounds(a)
    bx1, by1, bx2, by2 = rect_bounds(b)
    return max(0.0, min(ax2, bx2) - max(ax1, bx1)) * max(0.0, min(ay2, by2) - max(ay1, by1))


def orientation(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def proper_intersection(a, b, c, d):
    eps = 1e-9
    o1, o2 = orientation(a, b, c), orientation(a, b, d)
    o3, o4 = orientation(c, d, a), orientation(c, d, b)
    return ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and \
           ((o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps))


def point_segment_distance(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    length2 = dx * dx + dy * dy
    if length2 == 0:
        return math.hypot(p[0] - a[0], p[1] - a[1])
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2))
    q = (a[0] + t * dx, a[1] + t * dy)
    return math.hypot(p[0] - q[0], p[1] - q[1])


def segment_distance(a, b, c, d):
    if proper_intersection(a, b, c, d):
        return 0.0
    return min(point_segment_distance(a, c, d), point_segment_distance(b, c, d),
               point_segment_distance(c, a, b), point_segment_distance(d, a, b))


def point_rect_distance(p, bounds):
    x1, y1, x2, y2 = bounds
    dx = max(x1 - p[0], p[0] - x2, 0.0)
    dy = max(y1 - p[1], p[1] - y2, 0.0)
    return math.hypot(dx, dy)


def segment_rect_distance(a, b, bounds):
    x1, y1, x2, y2 = bounds
    corners = ((x1, y1), (x2, y1), (x2, y2), (x1, y2))
    edges = tuple(zip(corners, corners[1:] + corners[:1]))
    if x1 <= a[0] <= x2 and y1 <= a[1] <= y2:
        return 0.0
    if x1 <= b[0] <= x2 and y1 <= b[1] <= y2:
        return 0.0
    return min(point_rect_distance(a, bounds), point_rect_distance(b, bounds),
               *(segment_distance(a, b, c, d) for c, d in edges))


def bounds_gap(a, b):
    """Minimum distance between two axis-aligned bounding boxes."""
    dx = max(b[0] - a[2], a[0] - b[2], 0.0)
    dy = max(b[1] - a[3], a[1] - b[3], 0.0)
    return math.hypot(dx, dy)


def segment_bounds(a, b):
    return (min(a[0], b[0]), min(a[1], b[1]),
            max(a[0], b[0]), max(a[1], b[1]))


def evaluate(circuit: Circuit, positions: dict[str, tuple[float, float]]):
    terminals = terminal_positions(circuit, positions)
    segments = net_segments(circuit, terminals)
    length = sum(math.dist(a, b) for _, _, a, b in segments)
    length_raw = sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     for _, _, a, b in segments)
    overlaps = 0
    overlap_raw = 0.0
    clearance_violations = 0
    clearance_raw = 0.0
    for i, ref_a in enumerate(circuit.refs):
        for ref_b in circuit.refs[i + 1:]:
            area = overlap_area(positions[ref_a], positions[ref_b])
            if area > 0:
                overlaps += 1
                overlap_raw += 1.0 + (area / (WIDTH * HEIGHT)) ** 2
            else:
                gap = rect_gap(positions[ref_a], positions[ref_b])
                if gap < BODY_CLEARANCE:
                    clearance_violations += 1
                    clearance_raw += ((BODY_CLEARANCE - gap) / BODY_CLEARANCE) ** 2
    crossings = 0
    crowding_raw = 0.0
    segment_boxes = [segment_bounds(segment[2], segment[3]) for segment in segments]
    component_boxes = {ref: rect_bounds(positions[ref]) for ref in circuit.refs}
    for i, (net_a, term_a, a, b) in enumerate(segments):
        owner = term_a.rsplit(".", 1)[0]
        for ref in circuit.refs:
            if ref != owner:
                if bounds_gap(segment_boxes[i], component_boxes[ref]) >= NET_BODY_CLEARANCE:
                    continue
                distance = segment_rect_distance(a, b, component_boxes[ref])
                if distance < NET_BODY_CLEARANCE:
                    crowding_raw += ((NET_BODY_CLEARANCE - distance) / NET_BODY_CLEARANCE) ** 2
        for offset, (net_b, _, c, d) in enumerate(segments[i + 1:], i + 1):
            if net_a == net_b:
                continue
            if bounds_gap(segment_boxes[i], segment_boxes[offset]) >= NET_NET_CLEARANCE:
                continue
            if proper_intersection(a, b, c, d):
                crossings += 1
            distance = segment_distance(a, b, c, d)
            if distance < NET_NET_CLEARANCE:
                crowding_raw += ((NET_NET_CLEARANCE - distance) / NET_NET_CLEARANCE) ** 2
    raw = {"length": length_raw, "overlap": overlap_raw,
           "clearance": clearance_raw, "crossing": float(crossings),
           "crowding": crowding_raw}
    energy = {name: raw[name] * WEIGHTS[name] for name in WEIGHTS}
    return {
        "total_connection_length": length,
        "component_overlaps": overlaps,
        "clearance_violations": clearance_violations,
        "net_crossings": crossings,
        "energy": energy,
        "objective": sum(energy.values()),
    }


def evaluate_junction_state(circuit: Circuit, positions, junctions):
    terminals = terminal_positions(circuit, positions)
    segments = junction_segments(circuit, terminals, junctions)
    length = sum(math.dist(a, b) for _, _, a, b, _ in segments)
    length_raw = sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     for _, _, a, b, _ in segments)
    overlaps = 0
    overlap_raw = 0.0
    clearance_violations = 0
    clearance_raw = 0.0
    for i, ref_a in enumerate(circuit.refs):
        for ref_b in circuit.refs[i + 1:]:
            area = overlap_area(positions[ref_a], positions[ref_b])
            if area > 0:
                overlaps += 1
                overlap_raw += 1.0 + (area / (WIDTH * HEIGHT)) ** 2
            else:
                gap = rect_gap(positions[ref_a], positions[ref_b])
                if gap < BODY_CLEARANCE:
                    clearance_violations += 1
                    clearance_raw += ((BODY_CLEARANCE - gap) / BODY_CLEARANCE) ** 2
    crossings = 0
    crowding_raw = 0.0
    segment_boxes = [segment_bounds(segment[2], segment[3]) for segment in segments]
    component_boxes = {ref: rect_bounds(positions[ref]) for ref in circuit.refs}
    for i, (net_a, _, a, b, owners) in enumerate(segments):
        for ref in circuit.refs:
            if ref not in owners:
                if bounds_gap(segment_boxes[i], component_boxes[ref]) >= NET_BODY_CLEARANCE:
                    continue
                distance = segment_rect_distance(a, b, component_boxes[ref])
                if distance < NET_BODY_CLEARANCE:
                    crowding_raw += ((NET_BODY_CLEARANCE - distance) / NET_BODY_CLEARANCE) ** 2
        for offset, (net_b, _, c, d, _) in enumerate(segments[i + 1:], i + 1):
            if net_a == net_b:
                continue
            if bounds_gap(segment_boxes[i], segment_boxes[offset]) >= NET_NET_CLEARANCE:
                continue
            if proper_intersection(a, b, c, d):
                crossings += 1
            distance = segment_distance(a, b, c, d)
            if distance < NET_NET_CLEARANCE:
                crowding_raw += ((NET_NET_CLEARANCE - distance) / NET_NET_CLEARANCE) ** 2
    raw = {"length": length_raw, "overlap": overlap_raw,
           "clearance": clearance_raw, "crossing": float(crossings),
           "crowding": crowding_raw}
    energy = {name: raw[name] * WEIGHTS[name] for name in WEIGHTS}
    return {
        "total_connection_length": length,
        "component_overlaps": overlaps,
        "clearance_violations": clearance_violations,
        "net_crossings": crossings,
        "energy": energy,
        "objective": sum(energy.values()),
    }


def relax(circuit: Circuit, start):
    positions = dict(start)
    current = evaluate(circuit, positions)["objective"]
    step = INITIAL_STEP
    sweeps = 0
    while sweeps < MAX_SWEEPS:
        sweeps += 1
        accepted = 0
        for ref in circuit.refs:
            x, y = positions[ref]
            best_point, best_energy = (x, y), current
            for mx, my in MOVES:
                candidate = (x + mx * step, y + my * step)
                positions[ref] = candidate
                candidate_energy = evaluate(circuit, positions)["objective"]
                if candidate_energy < best_energy - 1e-9:
                    best_point, best_energy = candidate, candidate_energy
            positions[ref] = best_point
            if best_energy < current - 1e-9:
                current = best_energy
                accepted += 1
        if accepted == 0:
            if step <= MIN_STEP:
                return positions, sweeps, "no improving move at minimum step"
            step = max(MIN_STEP, step / 2.0)
    return positions, sweeps, "maximum sweeps reached"


def anneal(circuit: Circuit, start, *, seed=ANNEAL_SEED,
           iterations=ANNEAL_ITERATIONS,
           initial_temperature=INITIAL_TEMPERATURE,
           final_temperature=FINAL_TEMPERATURE,
           initial_move_distance=INITIAL_MOVE_DISTANCE,
           final_move_distance=FINAL_MOVE_DISTANCE):
    """Metropolis annealing with geometric cooling and movement schedules."""
    if iterations < 2:
        raise ValueError("annealing requires at least two iterations")
    rng = random.Random(seed)
    state = dict(start)
    energy = evaluate(circuit, state)["objective"]
    best = dict(state)
    best_energy = energy
    downhill = uphill = rejected = 0
    for iteration in range(iterations):
        fraction = iteration / (iterations - 1)
        temperature = initial_temperature * (final_temperature / initial_temperature) ** fraction
        move_distance = initial_move_distance * (final_move_distance / initial_move_distance) ** fraction
        ref = circuit.refs[rng.randrange(len(circuit.refs))]
        old = state[ref]
        state[ref] = (old[0] + rng.uniform(-move_distance, move_distance),
                      old[1] + rng.uniform(-move_distance, move_distance))
        candidate_energy = evaluate(circuit, state)["objective"]
        delta = candidate_energy - energy
        if delta <= 0.0:
            accepted = True
            downhill += 1
        elif rng.random() < math.exp(-delta / temperature):
            accepted = True
            uphill += 1
        else:
            accepted = False
            rejected += 1
        if accepted:
            energy = candidate_energy
            if energy < best_energy:
                best_energy = energy
                best = dict(state)
        else:
            state[ref] = old
    return best, {
        "accepted_downhill_moves": downhill,
        "accepted_uphill_moves": uphill,
        "rejected_moves": rejected,
        "best_energy": best_energy,
    }


def junction_movable_objects(circuit, junctions):
    return tuple(("component", ref) for ref in circuit.refs) + \
           tuple(("junction", name) for name in sorted(junctions, key=natural_key))


def anneal_junction_state(circuit, start_positions, start_junctions, *,
                          seed=ANNEAL_SEED, iterations=ANNEAL_ITERATIONS,
                          initial_temperature=INITIAL_TEMPERATURE,
                          final_temperature=FINAL_TEMPERATURE,
                          initial_move_distance=INITIAL_MOVE_DISTANCE,
                          final_move_distance=FINAL_MOVE_DISTANCE):
    if iterations < 2:
        raise ValueError("annealing requires at least two iterations")
    rng = random.Random(seed)
    positions = dict(start_positions)
    junctions = dict(start_junctions)
    energy = evaluate_junction_state(circuit, positions, junctions)["objective"]
    best_positions, best_junctions, best_energy = dict(positions), dict(junctions), energy
    movable = junction_movable_objects(circuit, junctions)
    downhill = uphill = rejected = 0
    for iteration in range(iterations):
        fraction = iteration / (iterations - 1)
        temperature = initial_temperature * (final_temperature / initial_temperature) ** fraction
        distance = initial_move_distance * (final_move_distance / initial_move_distance) ** fraction
        kind, name = movable[rng.randrange(len(movable))]
        collection = positions if kind == "component" else junctions
        old = collection[name]
        collection[name] = (old[0] + rng.uniform(-distance, distance),
                            old[1] + rng.uniform(-distance, distance))
        candidate = evaluate_junction_state(circuit, positions, junctions)["objective"]
        delta = candidate - energy
        if delta <= 0.0:
            accepted = True
            downhill += 1
        elif rng.random() < math.exp(-delta / temperature):
            accepted = True
            uphill += 1
        else:
            accepted = False
            rejected += 1
        if accepted:
            energy = candidate
            if energy < best_energy:
                best_positions, best_junctions = dict(positions), dict(junctions)
                best_energy = energy
        else:
            collection[name] = old
    return best_positions, best_junctions, {
        "accepted_downhill_moves": downhill,
        "accepted_uphill_moves": uphill,
        "rejected_moves": rejected,
        "best_energy": best_energy,
    }


def relax_junction_state(circuit, start_positions, start_junctions):
    positions, junctions = dict(start_positions), dict(start_junctions)
    current = evaluate_junction_state(circuit, positions, junctions)["objective"]
    movable = junction_movable_objects(circuit, junctions)
    step = INITIAL_STEP
    sweeps = 0
    while sweeps < MAX_SWEEPS:
        sweeps += 1
        accepted = 0
        for kind, name in movable:
            collection = positions if kind == "component" else junctions
            x, y = collection[name]
            best_point, best_energy = (x, y), current
            for mx, my in MOVES:
                candidate = (x + mx * step, y + my * step)
                collection[name] = candidate
                candidate_energy = evaluate_junction_state(circuit, positions, junctions)["objective"]
                if candidate_energy < best_energy - 1e-9:
                    best_point, best_energy = candidate, candidate_energy
            collection[name] = best_point
            if best_energy < current - 1e-9:
                current = best_energy
                accepted += 1
        if accepted == 0:
            if step <= MIN_STEP:
                return positions, junctions, sweeps, "no improving move at minimum step"
            step = max(MIN_STEP, step / 2.0)
    return positions, junctions, sweeps, "maximum sweeps reached"


def rounded(value):
    return round(value, ROUND_DIGITS)


def rounded_metrics(metrics):
    return {
        "total_connection_length": rounded(metrics["total_connection_length"]),
        "component_overlaps": metrics["component_overlaps"],
        "clearance_violations": metrics["clearance_violations"],
        "net_crossings": metrics["net_crossings"],
        "energy": {k: rounded(v) for k, v in metrics["energy"].items()},
        "objective": rounded(metrics["objective"]),
    }


def coordinate_document(circuit, positions):
    terminals = terminal_positions(circuit, positions)
    return {
        "circuit": circuit.name,
        "components": {
            ref: {
                "x": rounded(positions[ref][0]), "y": rounded(positions[ref][1]),
                "width": WIDTH, "height": HEIGHT,
                "terminals": {
                    pin: {"x": rounded(terminals[f"{ref}.{pin}"][0]),
                          "y": rounded(terminals[f"{ref}.{pin}"][1])}
                    for pin in circuit.pins[ref]
                },
            } for ref in circuit.refs
        },
    }


def junction_coordinate_document(circuit, positions, junctions):
    document = coordinate_document(circuit, positions)
    document["junctions"] = {
        name: {"x": rounded(point[0]), "y": rounded(point[1])}
        for name, point in sorted(junctions.items(), key=lambda item: natural_key(item[0]))
    }
    return document


def positions_from_document(circuit, document):
    components = document["components"]
    if set(components) != set(circuit.refs):
        raise ValueError("coordinate document component set does not match circuit")
    return {ref: (float(components[ref]["x"]), float(components[ref]["y"]))
            for ref in circuit.refs}


def write_json(path, document):
    path.write_text(json.dumps(document, sort_keys=True, indent=2,
                               ensure_ascii=False) + "\n", encoding="utf-8")


def escape_xml(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_svg(path, circuit, positions):
    terminals = terminal_positions(circuit, positions)
    segments = net_segments(circuit, terminals)
    xs = [p[0] for p in positions.values()]
    ys = [p[1] for p in positions.values()]
    left, top = min(xs) - WIDTH, min(ys) - HEIGHT
    width, height = max(xs) - min(xs) + 2 * WIDTH, max(ys) - min(ys) + 2 * HEIGHT
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="{left:.3f} {top:.3f} {width:.3f} {height:.3f}">',
        '<rect x="{:.3f}" y="{:.3f}" width="{:.3f}" height="{:.3f}" fill="white"/>'.format(left, top, width, height),
        '<g stroke="#777" stroke-width="1" fill="none">',
    ]
    for _, _, a, b in segments:
        lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    lines.append('</g><g stroke="#111" stroke-width="1.5" fill="#f5f5f5">')
    for ref in circuit.refs:
        x, y = positions[ref]
        lines.append(f'<rect x="{x-WIDTH/2:.3f}" y="{y-HEIGHT/2:.3f}" width="{WIDTH:.3f}" height="{HEIGHT:.3f}"/>')
    lines.append('</g><g fill="#111" font-family="monospace" font-size="10" text-anchor="middle" dominant-baseline="middle">')
    for ref in circuit.refs:
        x, y = positions[ref]
        lines.append(f'<text x="{x:.3f}" y="{y:.3f}">{escape_xml(ref)}</text>')
    lines.append('</g><g fill="#d22" stroke="none">')
    for x, y in terminals.values():
        lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.5"/>')
    lines.append('</g></svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_junction_svg(path, circuit, positions, junctions):
    terminals = terminal_positions(circuit, positions)
    segments = junction_segments(circuit, terminals, junctions)
    all_points = list(positions.values()) + list(junctions.values())
    xs, ys = [p[0] for p in all_points], [p[1] for p in all_points]
    left, top = min(xs) - WIDTH, min(ys) - HEIGHT
    width, height = max(xs) - min(xs) + 2 * WIDTH, max(ys) - min(ys) + 2 * HEIGHT
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="{left:.3f} {top:.3f} {width:.3f} {height:.3f}">',
        '<rect x="{:.3f}" y="{:.3f}" width="{:.3f}" height="{:.3f}" fill="white"/>'.format(left, top, width, height),
        '<g stroke="#777" stroke-width="1" fill="none">',
    ]
    for _, _, a, b, _ in segments:
        lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    lines.append('</g><g stroke="#111" stroke-width="1.5" fill="#f5f5f5">')
    for ref in circuit.refs:
        x, y = positions[ref]
        lines.append(f'<rect x="{x-WIDTH/2:.3f}" y="{y-HEIGHT/2:.3f}" width="{WIDTH:.3f}" height="{HEIGHT:.3f}"/>')
    lines.append('</g><g fill="#111" font-family="monospace" font-size="10" text-anchor="middle" dominant-baseline="middle">')
    for ref in circuit.refs:
        x, y = positions[ref]
        lines.append(f'<text x="{x:.3f}" y="{y:.3f}">{escape_xml(ref)}</text>')
    lines.append('</g><g fill="#d22" stroke="none">')
    for x, y in terminals.values():
        lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.5"/>')
    lines.append('</g><g fill="#1769aa" stroke="#083b66" stroke-width="0.8">')
    for x, y in junctions.values():
        lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3.5"/>')
    lines.append('</g></svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--anneal", action="store_true",
                        help="run frozen Experiment 2 instead of Experiment 1")
    parser.add_argument("--junction-anneal", action="store_true",
                        help="run frozen Experiment 3 movable-junction model")
    parser.add_argument("--degree-seeded", action="store_true",
                        help="run Experiment 4A degree-seeded initialization")
    args = parser.parse_args()
    circuit = parse_circuit(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    initial = initial_positions(circuit)
    initial_metrics = evaluate(circuit, initial)
    if sum((args.anneal, args.junction_anneal, args.degree_seeded)) > 1:
        parser.error("choose only one experiment mode")
    if args.degree_seeded:
        baseline_one = json.loads((args.output_dir / "metrics.json").read_text(encoding="utf-8"))
        baseline_two = json.loads((args.output_dir / "annealed_metrics.json").read_text(encoding="utf-8"))
        baseline_three = json.loads(
            (args.output_dir / "junction_annealed_metrics.json").read_text(encoding="utf-8"))
        seeded = degree_seeded_positions(circuit)
        seeded_metrics = evaluate(circuit, seeded)
        best, stats = anneal(circuit, seeded)
        final, quench_sweeps, reason = relax(circuit, best)
        final_metrics = evaluate(circuit, final)
        write_svg(args.output_dir / "degree_seeded_initial.svg", circuit, seeded)
        write_svg(args.output_dir / "degree_seeded_annealed.svg", circuit, final)
        coordinates_path = args.output_dir / "degree_seeded_coordinates.json"
        write_json(coordinates_path, coordinate_document(circuit, final))
        degrees = component_degrees(circuit)
        degree_table = [{"ref": ref, "degree": degrees[ref]}
                        for ref in sorted(circuit.refs,
                                          key=lambda ref: (-degrees[ref], ref))]
        maximum_degree = max(degrees.values())
        report = {
            "input": {"components": len(circuit.refs), "nets": len(circuit.nets),
                      "singleton_nets": sum(len(n.terminals) == 1 for n in circuit.nets)},
            "degree_table": degree_table,
            "maximum_degree": maximum_degree,
            "maximum_degree_components": sorted(
                ref for ref, degree in degrees.items() if degree == maximum_degree),
            "ring_spacing": DEGREE_RING_SPACING,
            "degree_seeded_initial": rounded_metrics(seeded_metrics),
            "states": {
                "original_initial": baseline_one["initial"],
                "experiment_1_final": baseline_one["final"],
                "experiment_2_annealed_and_quenched":
                    baseline_two["states"]["experiment_2_annealed_and_quenched"],
                "experiment_3_junction_annealed_and_quenched":
                    baseline_three["states"]["experiment_3_junction_annealed_and_quenched"],
                "experiment_4a_degree_seeded_annealed_and_quenched":
                    rounded_metrics(final_metrics),
            },
            "annealing": {
                "seed": ANNEAL_SEED, "initial_temperature": INITIAL_TEMPERATURE,
                "final_temperature": FINAL_TEMPERATURE, "cooling_schedule": "geometric",
                "initial_move_distance": INITIAL_MOVE_DISTANCE,
                "final_move_distance": FINAL_MOVE_DISTANCE,
                "move_distance_schedule": "geometric", "iterations": ANNEAL_ITERATIONS,
                "accepted_downhill_moves": stats["accepted_downhill_moves"],
                "accepted_uphill_moves": stats["accepted_uphill_moves"],
                "rejected_moves": stats["rejected_moves"],
                "best_energy_encountered": rounded(stats["best_energy"]),
            },
            "quench": {"final_energy": rounded(final_metrics["objective"]),
                       "sweeps": quench_sweeps, "stopping_reason": reason},
            "degree_seeded_coordinate_sha256":
                hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
        }
        write_json(args.output_dir / "degree_seeded_metrics.json", report)
        print(json.dumps(report, sort_keys=True, indent=2))
        return
    if args.junction_anneal:
        baseline_one = json.loads((args.output_dir / "metrics.json").read_text(encoding="utf-8"))
        baseline_two = json.loads((args.output_dir / "annealed_metrics.json").read_text(encoding="utf-8"))
        junctions = initial_junctions(circuit, initial)
        junction_initial_metrics = evaluate_junction_state(circuit, initial, junctions)
        best_positions, best_junctions, stats = anneal_junction_state(circuit, initial, junctions)
        final, final_junctions, quench_sweeps, reason = relax_junction_state(
            circuit, best_positions, best_junctions)
        final_metrics = evaluate_junction_state(circuit, final, final_junctions)
        coordinates_path = args.output_dir / "junction_annealed_coordinates.json"
        write_json(coordinates_path, junction_coordinate_document(circuit, final, final_junctions))
        write_junction_svg(args.output_dir / "junction_annealed.svg", circuit, final, final_junctions)
        report = {
            "input": {"components": len(circuit.refs), "nets": len(circuit.nets),
                      "singleton_nets": sum(len(n.terminals) == 1 for n in circuit.nets),
                      "movable_junctions": len(junctions)},
            "states": {
                "initial_junction_model": rounded_metrics(junction_initial_metrics),
                "experiment_1_final": baseline_one["final"],
                "experiment_2_annealed_and_quenched":
                    baseline_two["states"]["experiment_2_annealed_and_quenched"],
                "experiment_3_junction_annealed_and_quenched": rounded_metrics(final_metrics),
            },
            "annealing": {
                "seed": ANNEAL_SEED, "initial_temperature": INITIAL_TEMPERATURE,
                "final_temperature": FINAL_TEMPERATURE, "cooling_schedule": "geometric",
                "initial_move_distance": INITIAL_MOVE_DISTANCE,
                "final_move_distance": FINAL_MOVE_DISTANCE,
                "move_distance_schedule": "geometric", "iterations": ANNEAL_ITERATIONS,
                "accepted_downhill_moves": stats["accepted_downhill_moves"],
                "accepted_uphill_moves": stats["accepted_uphill_moves"],
                "rejected_moves": stats["rejected_moves"],
                "best_energy_encountered": rounded(stats["best_energy"]),
            },
            "quench": {"final_energy": rounded(final_metrics["objective"]),
                       "sweeps": quench_sweeps, "stopping_reason": reason},
            "junction_annealed_coordinate_sha256":
                hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
        }
        write_json(args.output_dir / "junction_annealed_metrics.json", report)
        print(json.dumps(report, sort_keys=True, indent=2))
        return
    if args.anneal:
        baseline_path = args.output_dir / "final_coordinates.json"
        baseline_metrics_path = args.output_dir / "metrics.json"
        if not baseline_path.exists():
            raise ValueError(f"Experiment 1 baseline is missing: {baseline_path}")
        if not baseline_metrics_path.exists():
            raise ValueError(f"Experiment 1 metrics are missing: {baseline_metrics_path}")
        baseline_document = json.loads(baseline_path.read_text(encoding="utf-8"))
        baseline = positions_from_document(circuit, baseline_document)
        # Validate that the coordinate artifact still matches the circuit, but
        # report the frozen baseline metrics rather than recomputing from its
        # deliberately rounded serialized coordinates.
        evaluate(circuit, baseline)
        baseline_report = json.loads(baseline_metrics_path.read_text(encoding="utf-8"))
        best, anneal_stats = anneal(circuit, initial)
        final, quench_sweeps, quench_reason = relax(circuit, best)
        final_metrics = evaluate(circuit, final)
        coordinates_path = args.output_dir / "annealed_coordinates.json"
        write_json(coordinates_path, coordinate_document(circuit, final))
        write_svg(args.output_dir / "annealed.svg", circuit, final)
        report = {
            "input": {"components": len(circuit.refs), "nets": len(circuit.nets),
                      "singleton_nets": sum(len(n.terminals) == 1 for n in circuit.nets)},
            "states": {
                "initial": rounded_metrics(initial_metrics),
                "experiment_1_final": baseline_report["final"],
                "experiment_2_annealed_and_quenched": rounded_metrics(final_metrics),
            },
            "annealing": {
                "seed": ANNEAL_SEED,
                "initial_temperature": INITIAL_TEMPERATURE,
                "final_temperature": FINAL_TEMPERATURE,
                "cooling_schedule": "geometric",
                "initial_move_distance": INITIAL_MOVE_DISTANCE,
                "final_move_distance": FINAL_MOVE_DISTANCE,
                "move_distance_schedule": "geometric",
                "iterations": ANNEAL_ITERATIONS,
                "accepted_downhill_moves": anneal_stats["accepted_downhill_moves"],
                "accepted_uphill_moves": anneal_stats["accepted_uphill_moves"],
                "rejected_moves": anneal_stats["rejected_moves"],
                "best_energy_encountered": rounded(anneal_stats["best_energy"]),
            },
            "quench": {"final_energy": rounded(final_metrics["objective"]),
                       "sweeps": quench_sweeps, "stopping_reason": quench_reason},
            "annealed_coordinate_sha256": hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
        }
        write_json(args.output_dir / "annealed_metrics.json", report)
        print(json.dumps(report, sort_keys=True, indent=2))
        return
    final, sweeps, reason = relax(circuit, initial)
    final_metrics = evaluate(circuit, final)
    write_svg(args.output_dir / "initial.svg", circuit, initial)
    write_svg(args.output_dir / "final.svg", circuit, final)
    coordinates_path = args.output_dir / "final_coordinates.json"
    write_json(coordinates_path, coordinate_document(circuit, final))
    report = {
        "input": {"components": len(circuit.refs), "nets": len(circuit.nets),
                  "singleton_nets": sum(len(n.terminals) == 1 for n in circuit.nets)},
        "parameters": {"seed": SEED, "weights": WEIGHTS, "max_sweeps": MAX_SWEEPS,
                       "initial_step": INITIAL_STEP, "minimum_step": MIN_STEP},
        "initial": rounded_metrics(initial_metrics),
        "final": rounded_metrics(final_metrics),
        "connection_length_change_percent": rounded(
            100.0 * (final_metrics["total_connection_length"] - initial_metrics["total_connection_length"])
            / initial_metrics["total_connection_length"]),
        "stopping_reason": reason,
        "sweeps": sweeps,
        "final_coordinate_sha256": hashlib.sha256(coordinates_path.read_bytes()).hexdigest(),
    }
    write_json(args.output_dir / "metrics.json", report)
    print(json.dumps(report, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
