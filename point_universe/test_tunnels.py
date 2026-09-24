"""Tunnels: short splice wires that pass only over wires, at 90 degrees, and can never be crossed."""

import pytest

from .universe import Line, PointUniverse, Square


def base():
    u = PointUniverse(14, 14, line_separation=0, square_clearance=0)
    u.add_line(Line("W", "A", ((6, 2), (6, 10))))
    return u


def test_tunnel_over_one_wire_is_2q_long_and_legal():
    u = base()
    u.add_line(Line("T", "B", ((5, 5), (7, 5)), tunnel=True))
    assert u.point_lines[(6, 5)] == {"W", "T"}          # two lines at the crossing point, no more


def test_tunnel_over_two_adjacent_wires_is_3q_long():
    u = base()
    u.add_line(Line("W2", "C", ((7, 2), (7, 10))))
    u.add_line(Line("T", "B", ((5, 5), (8, 5)), tunnel=True))
    assert u.point_lines[(7, 5)] == {"W2", "T"}


def test_tunnel_may_not_cross_a_device():
    u = base()
    u.add_square(Square("D", 9, 4, 2, 2))
    with pytest.raises(ValueError):
        u.add_line(Line("T", "B", ((5, 5), (11, 5)), tunnel=True))


def test_tunnel_may_not_cross_at_a_bend_or_wire_end():
    u = PointUniverse(14, 14, line_separation=0, square_clearance=0)
    u.add_line(Line("W", "A", ((6, 2), (6, 5), (9, 5))))
    with pytest.raises(ValueError):
        u.add_line(Line("T", "B", ((5, 5), (7, 5)), tunnel=True))   # bend at (6,5)


def test_a_tunnel_can_never_be_crossed_by_a_wire():
    u = base()
    u.add_line(Line("T", "B", ((5, 5), (7, 5)), tunnel=True))
    with pytest.raises(ValueError):
        u.add_line(Line("X", "C", ((5, 3), (5, 8))))               # runs through the tunnel's end point


def test_a_tunnel_can_never_be_crossed_by_a_tunnel():
    u = base()
    u.add_line(Line("T", "B", ((3, 5), (9, 5)), tunnel=True)) if False else None
    u.add_line(Line("T", "B", ((5, 5), (7, 5)), tunnel=True))
    with pytest.raises(ValueError):
        u.add_line(Line("T2", "C", ((6, 4), (6, 6)), tunnel=True))


def test_a_tunnel_may_not_cross_its_own_net():
    u = base()
    with pytest.raises(ValueError):
        u.add_line(Line("T", "A", ((5, 5), (7, 5)), tunnel=True))


def test_tunnel_is_one_straight_run():
    with pytest.raises(ValueError):
        Line("T", "B", ((5, 5), (7, 5), (7, 7)), tunnel=True)
    with pytest.raises(ValueError):
        Line("T", "B", ((5, 5), (6, 5)), tunnel=True)
