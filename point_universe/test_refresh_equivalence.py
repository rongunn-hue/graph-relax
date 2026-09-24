"""The fast empty-space refresh must give exactly the same matrix as the authoritative pairwise definition
(PointUniverse._empty_point_usable). Deterministic: fixed seeds, no randomness beyond them."""

import random
import time

from .universe import Line, LineEndpoint, PointState, PointUniverse, Square


def reference_states(universe):
    """Recompute every point's WHITE/UNUSABLE flag directly from _empty_point_usable."""
    out = {}
    for point in universe.all_points():
        state = universe.point_state(point)
        occupied = bool(state & (PointState.SQUARE | PointState.LINE | PointState.ATTACHMENT))
        state &= ~(PointState.WHITE | PointState.UNUSABLE)
        if not occupied:
            state |= PointState.WHITE if universe._empty_point_usable(point) else PointState.UNUSABLE
        out[point] = state
    return out


def random_universe(seed):
    rng = random.Random(seed)
    width, height = rng.randint(8, 26), rng.randint(8, 26)
    universe = PointUniverse(width, height, line_separation=rng.randint(0, 3), square_clearance=rng.randint(0, 3))
    for i in range(rng.randint(0, 6)):
        w, h = rng.randint(2, 5), rng.randint(2, 5)
        square = Square(f"S{i}", rng.randint(0, max(0, width - w - 1)), rng.randint(0, max(0, height - h - 1)), w, h,
                        rotation=rng.choice((0, 90, 180, 270)))
        try:
            universe.add_square(square)
        except ValueError:
            pass
    for i in range(rng.randint(0, 6)):
        x, y = rng.randint(0, width - 1), rng.randint(0, height - 1)
        x2, y2 = rng.randint(0, width - 1), rng.randint(0, height - 1)
        points = tuple(dict.fromkeys(((x, y), (x2, y), (x2, y2))))
        if len(points) < 2:
            continue
        try:
            universe.add_line(Line(f"L{i}", f"N{i}", points))
        except ValueError:
            pass
    return universe


def test_fast_refresh_matches_pairwise_definition_on_many_universes():
    checked = 0
    for seed in range(300):
        universe = random_universe(seed)
        expected = reference_states(universe)
        for point in universe.all_points():
            assert universe.point_state(point) == expected[point], (seed, point)
        checked += 1
    assert checked == 300


def test_fast_refresh_matches_after_every_commit():
    universe = PointUniverse(20, 14, line_separation=1, square_clearance=2)
    steps = [
        ("square", Square("A", 3, 3, 4, 3)),
        ("square", Square("B", 12, 6, 3, 4)),
        ("line", Line("L1", "N1", ((7, 4), (10, 4), (10, 8), (12, 8)))),
    ]
    for kind, obj in steps:
        (universe.add_square if kind == "square" else universe.add_line)(obj)
        expected = reference_states(universe)
        for point in universe.all_points():
            assert universe.point_state(point) == expected[point], (kind, point)


def test_refresh_is_fast_enough_for_a_schematic_sized_grid():
    # HRNG-scale: ~130 x 90 points. The pairwise definition would need ~10^8+ operations per commit.
    universe = PointUniverse(130, 90, line_separation=1, square_clearance=1)
    start = time.time()
    for i in range(20):
        universe.add_square(Square(f"S{i}", 4 + (i % 5) * 24, 4 + (i // 5) * 20, 8, 6))
    elapsed = time.time() - start
    assert elapsed < 20.0, elapsed
