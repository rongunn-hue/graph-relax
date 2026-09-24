"""Discrete point-universe prototype.

The authoritative geometry is a finite matrix of integer-coordinate points.
Square and Line objects describe proposed geometry; PointUniverse validates and
commits that geometry into the matrix.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from enum import IntFlag, auto
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Set, Tuple


Point = Tuple[int, int]


class PointState(IntFlag):
    SQUARE = auto()
    LINE = auto()
    ATTACHMENT = auto()
    WHITE = auto()
    UNUSABLE = auto()


@dataclass(frozen=True)
class Attachment:
    attachment_id: str
    square_id: str
    side: str
    offset: int


@dataclass
class Square:
    square_id: str
    x: int
    y: int
    width: int
    height: int
    rotation: int = 0
    attachments: List[Attachment] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.width < 2 or self.height < 2:
            raise ValueError("minimum square size is 2q x 2q")
        if self.rotation not in (0, 90, 180, 270):
            raise ValueError("rotation must be one of 0, 90, 180, 270")

    @property
    def effective_width(self) -> int:
        return self.height if self.rotation in (90, 270) else self.width

    @property
    def effective_height(self) -> int:
        return self.width if self.rotation in (90, 270) else self.height

    def occupied_points(self) -> Set[Point]:
        return {
            (px, py)
            for px in range(self.x, self.x + self.effective_width + 1)
            for py in range(self.y, self.y + self.effective_height + 1)
        }

    def legal_attachment_offsets(self, side: str) -> range:
        side_length = self.width if side in ("top", "bottom") else self.height
        if side not in ("top", "right", "bottom", "left"):
            raise ValueError(f"unknown side: {side}")
        return range(1, side_length)

    def add_attachment(self, attachment_id: str, side: str, offset: int) -> Attachment:
        if offset not in self.legal_attachment_offsets(side):
            raise ValueError("attachment must be at least 1q from either side vertex")
        attachment = Attachment(attachment_id, self.square_id, side, offset)
        self.attachments.append(attachment)
        return attachment

    def attachment_point(self, attachment: Attachment) -> Point:
        if attachment.square_id != self.square_id:
            raise ValueError("attachment belongs to another square")
        local = self._unrotated_local_attachment(attachment.side, attachment.offset)
        rotated = self._rotate_local_point(local)
        return (self.x + rotated[0], self.y + rotated[1])

    def _unrotated_local_attachment(self, side: str, offset: int) -> Point:
        if offset not in self.legal_attachment_offsets(side):
            raise ValueError("illegal attachment offset")
        if side == "top":
            return (offset, 0)
        if side == "right":
            return (self.width, offset)
        if side == "bottom":
            return (offset, self.height)
        if side == "left":
            return (0, offset)
        raise ValueError(f"unknown side: {side}")

    def _rotate_local_point(self, point: Point) -> Point:
        x, y = point
        if self.rotation == 0:
            return (x, y)
        if self.rotation == 90:
            return (self.height - y, x)
        if self.rotation == 180:
            return (self.width - x, self.height - y)
        return (y, self.width - x)


@dataclass(frozen=True)
class LineEndpoint:
    attachment_id: str = ""
    line_id: str = ""


@dataclass(frozen=True)
class Line:
    line_id: str
    net_id: str
    points: Tuple[Point, ...]
    endpoints: Tuple[LineEndpoint, LineEndpoint] = (LineEndpoint(), LineEndpoint())
    tunnel: bool = False

    def __post_init__(self) -> None:
        if len(self.points) < 2:
            raise ValueError("a line needs at least two points")
        if len(set(self.points)) != len(self.points):
            raise ValueError("duplicate consecutive or repeated line points are not allowed")
        for point in self.points:
            if not all(isinstance(value, int) for value in point):
                raise ValueError("line points must be integer grid points")
        for a, b in zip(self.points, self.points[1:]):
            if a[0] != b[0] and a[1] != b[1]:
                raise ValueError("line segments must be orthogonal")
        if self.tunnel:
            # A tunnel is a splice: one straight, very short wire, at least 2q long (it must pass over a wire).
            if len(self.points) != 2:
                raise ValueError("a tunnel is one straight run: exactly two points")
            a, b = self.points
            if abs(a[0] - b[0]) + abs(a[1] - b[1]) < 2:
                raise ValueError("a tunnel is at least 2q long")

    def occupied_points(self) -> Set[Point]:
        out: Set[Point] = set()
        for a, b in zip(self.points, self.points[1:]):
            if a[0] == b[0]:
                x = a[0]
                y0, y1 = sorted((a[1], b[1]))
                out.update((x, y) for y in range(y0, y1 + 1))
            else:
                y = a[1]
                x0, x1 = sorted((a[0], b[0]))
                out.update((x, y) for x in range(x0, x1 + 1))
        return out

    def endpoint_points(self) -> Tuple[Point, Point]:
        return (self.points[0], self.points[-1])


@dataclass(frozen=True)
class Connectivity:
    line_to_attachment: Dict[Tuple[str, int], str] = field(default_factory=dict)
    line_to_line: Dict[Tuple[str, int], str] = field(default_factory=dict)
    point_to_lines: Dict[Point, Set[str]] = field(default_factory=dict)
    net_to_lines: Dict[str, Set[str]] = field(default_factory=dict)


class PointUniverse:
    def __init__(self, width: int, height: int, line_separation: int = 1, square_clearance: int = 1):
        if width <= 0 or height <= 0:
            raise ValueError("universe dimensions must be positive")
        self.width = width
        self.height = height
        self.line_separation = line_separation
        self.square_clearance = square_clearance
        self.matrix: List[List[PointState]] = [
            [PointState.WHITE for _x in range(width)] for _y in range(height)
        ]
        self.squares: Dict[str, Square] = {}
        self.lines: Dict[str, Line] = {}
        self.attachments: Dict[str, Attachment] = {}
        self.attachment_points: Dict[str, Point] = {}
        self.line_points: Dict[str, Set[Point]] = {}
        self.point_lines: Dict[Point, Set[str]] = {}
        self.line_to_attachment: Dict[Tuple[str, int], str] = {}
        self.line_to_line: Dict[Tuple[str, int], str] = {}
        self.dependencies: Dict[str, Set[str]] = {}
        self.net_lines: Dict[str, Set[str]] = {}
        self._refresh_empty_space()

    def in_bounds(self, point: Point) -> bool:
        x, y = point
        return 0 <= x < self.width and 0 <= y < self.height

    def all_points(self) -> Iterable[Point]:
        for y in range(self.height):
            for x in range(self.width):
                yield (x, y)

    def point_state(self, point: Point) -> PointState:
        self._require_in_bounds(point)
        x, y = point
        return self.matrix[y][x]

    def layer(self, state: PointState) -> Set[Point]:
        return {
            (x, y)
            for y, row in enumerate(self.matrix)
            for x, value in enumerate(row)
            if value & state
        }

    @property
    def connectivity(self) -> Connectivity:
        return Connectivity(
            line_to_attachment=dict(self.line_to_attachment),
            line_to_line=dict(self.line_to_line),
            point_to_lines={point: set(lines) for point, lines in self.point_lines.items()},
            net_to_lines={net: set(lines) for net, lines in self.net_lines.items()},
        )

    def add_square(self, square: Square) -> None:
        if square.square_id in self.squares:
            raise ValueError(f"duplicate square: {square.square_id}")
        occupied = square.occupied_points()
        for point in occupied:
            self._require_in_bounds(point)
            if self.point_state(point) & PointState.SQUARE:
                raise ValueError(f"square overlap at {point}")
        attachment_points = {}
        for attachment in square.attachments:
            if attachment.square_id != square.square_id:
                raise ValueError("attachment belongs to another square")
            point = square.attachment_point(attachment)
            if point not in occupied:
                raise ValueError("attachment is not on the rotated square boundary")
            if point not in self._boundary_points(occupied):
                raise ValueError("attachment is not on square boundary")
            attachment_points[attachment.attachment_id] = point
        self.squares[square.square_id] = square
        for point in occupied:
            self._add_state(point, PointState.SQUARE)
        for attachment in square.attachments:
            self.attachments[attachment.attachment_id] = attachment
            self.attachment_points[attachment.attachment_id] = attachment_points[attachment.attachment_id]
            self._add_state(attachment_points[attachment.attachment_id], PointState.ATTACHMENT)
        self._maybe_refresh()

    def add_line(self, line: Line) -> None:
        if line.line_id in self.lines:
            raise ValueError(f"duplicate line: {line.line_id}")
        segment_points = expanded_segment_points(line)
        self._validate_line_self_intersections(line, segment_points)
        points = line.occupied_points()
        for point in points:
            self._require_in_bounds(point)
        endpoint_lookup = self._validate_line_endpoints(line)
        for point in points:
            if line.tunnel and point not in line.endpoint_points():
                self._validate_tunnel_interior_point(line, point)
            else:
                self._validate_line_point(line, point, endpoint_lookup)
        self.lines[line.line_id] = line
        self.line_points[line.line_id] = set(points)
        self.net_lines.setdefault(line.net_id, set()).add(line.line_id)
        for point in points:
            self._add_state(point, PointState.LINE)
            self.point_lines.setdefault(point, set()).add(line.line_id)
        for index, target in endpoint_lookup.items():
            endpoint = line.endpoints[index]
            if endpoint.attachment_id:
                self.line_to_attachment[(line.line_id, index)] = target
            elif endpoint.line_id:
                self.line_to_line[(line.line_id, index)] = target
            self.dependencies.setdefault(target, set()).add(line.line_id)
        self._maybe_refresh()

    def _validate_line_endpoints(self, line: Line) -> Dict[int, str]:
        endpoint_lookup: Dict[int, str] = {}
        endpoint_points = line.endpoint_points()
        for index, endpoint in enumerate(line.endpoints):
            point = endpoint_points[index]
            target = endpoint.attachment_id or endpoint.line_id
            if not target:
                continue
            if endpoint.attachment_id:
                expected = self.attachment_points.get(endpoint.attachment_id)
                if expected is None:
                    raise ValueError(f"unknown attachment endpoint: {endpoint.attachment_id}")
                if point != expected:
                    raise ValueError("line endpoint does not coincide with declared attachment")
            if endpoint.line_id:
                if endpoint.line_id not in self.lines:
                    raise ValueError(f"unknown line endpoint: {endpoint.line_id}")
                if point not in self.line_points[endpoint.line_id]:
                    raise ValueError("line endpoint does not coincide with declared line")
                other = self.lines[endpoint.line_id]
                if other.net_id != line.net_id:
                    raise ValueError("line-to-line endpoint must be same net")
            endpoint_lookup[index] = target
        return endpoint_lookup

    def _validate_line_point(self, line: Line, point: Point, endpoint_lookup: Dict[int, str]) -> None:
        state = self.point_state(point)
        endpoint_points = line.endpoint_points()
        is_endpoint = point in endpoint_points
        endpoint_indices = {index for index, endpoint_point in enumerate(endpoint_points) if endpoint_point == point}
        if state & PointState.SQUARE:
            if not (is_endpoint and state & PointState.ATTACHMENT):
                raise ValueError(f"line passes through square interior or boundary at {point}")
            if not any(index in endpoint_lookup for index in endpoint_indices):
                raise ValueError("line may enter square only at declared attachment")
        existing_lines = self.point_lines.get(point, set())
        if existing_lines:
            for other_id in existing_lines:
                other = self.lines[other_id]
                if other.net_id != line.net_id:
                    raise ValueError(f"unrelated line overlap at {point}")
            if len(existing_lines | {line.line_id}) > 2:
                raise ValueError(f"more than two lines at point {point}")
            if not is_endpoint:
                raise ValueError("line may join another line only at a grid endpoint")
        if line.net_id:
            for neighbor in points_within_manhattan(point, self.line_separation):
                if not self.in_bounds(neighbor):
                    continue
                for other_id in self.point_lines.get(neighbor, set()):
                    other = self.lines[other_id]
                    if other.net_id != line.net_id:
                        raise ValueError(f"line-spacing violation near {point}")

    def _validate_tunnel_interior_point(self, line: Line, point: Point) -> None:
        """A tunnel may pass only over WIRES (or empty space between them): never a square, an attachment or another
        tunnel; only at 90 degrees on a straight stretch of a wire of another net; and at most two lines per point."""
        state = self.point_state(point)
        if state & (PointState.SQUARE | PointState.ATTACHMENT):
            raise ValueError(f"a tunnel may cross only wires (device at {point})")
        existing = self.point_lines.get(point, set())
        if not existing:
            return
        if len(existing) != 1:
            raise ValueError(f"more than two lines at point {point}")
        (other_id,) = existing
        other = self.lines[other_id]
        if other.tunnel:
            raise ValueError(f"a tunnel can never be crossed (at {point})")
        if other.net_id == line.net_id:
            raise ValueError(f"a tunnel crosses wires of other nets only (at {point})")
        a, b = line.points
        horizontal = a[1] == b[1]
        along = (1, 0) if horizontal else (0, 1)
        across = (0, 1) if horizontal else (1, 0)
        pts = self.line_points[other_id]
        if not ((point[0] + across[0], point[1] + across[1]) in pts and (point[0] - across[0], point[1] - across[1]) in pts):
            raise ValueError(f"a tunnel crosses a wire only at 90 degrees on a straight stretch (at {point})")
        if (point[0] + along[0], point[1] + along[1]) in pts or (point[0] - along[0], point[1] - along[1]) in pts:
            raise ValueError(f"a tunnel may not lie along a wire (at {point})")

    def _validate_line_self_intersections(self, line: Line, segment_points: List[Set[Point]]) -> None:
        point_segments: Dict[Point, Set[int]] = {}
        for index, points in enumerate(segment_points):
            for point in points:
                point_segments.setdefault(point, set()).add(index)
        for point, indices in point_segments.items():
            if len(indices) <= 1:
                continue
            ordered = sorted(indices)
            if len(ordered) == 2 and ordered[1] == ordered[0] + 1 and point == line.points[ordered[1]]:
                continue
            raise ValueError(f"line self-intersection or self-retrace at {point}")

    defer_refresh = False

    def _maybe_refresh(self) -> None:
        if not self.defer_refresh:
            self._refresh_empty_space()

    def _refresh_empty_space(self) -> None:
        # Same rule as _empty_point_usable (which remains the authoritative definition and is what the
        # equivalence test compares against): an empty point is UNUSABLE iff it lies within
        # square_clearance (Manhattan) of a SQUARE point or within line_separation of a LINE point.
        # Computed for the whole matrix at once instead of scanning every point against every occupied point.
        near_square = self._within_manhattan(self.layer(PointState.SQUARE), self.square_clearance)
        near_line = self._within_manhattan(self.layer(PointState.LINE), self.line_separation)
        for y in range(self.height):
            row = self.matrix[y]
            for x in range(self.width):
                state = row[x]
                occupied = bool(state & (PointState.SQUARE | PointState.LINE | PointState.ATTACHMENT))
                state &= ~(PointState.WHITE | PointState.UNUSABLE)
                if not occupied:
                    point = (x, y)
                    state |= PointState.UNUSABLE if (point in near_square or point in near_line) else PointState.WHITE
                row[x] = state

    def _within_manhattan(self, sources: Set[Point], radius: int) -> Set[Point]:
        """Every in-bounds point whose Manhattan distance to some source point is <= radius.

        Multi-source breadth-first expansion. In a full rectangle of integer points, four-neighbor
        step count equals Manhattan distance, so this is exactly the set the pairwise definition gives."""
        if radius < 0 or not sources:
            return set()
        reached = set(sources)
        frontier = list(sources)
        for _ in range(radius):
            next_frontier = []
            for point in frontier:
                for neighbor in four_neighbors(point):
                    if neighbor not in reached and self.in_bounds(neighbor):
                        reached.add(neighbor)
                        next_frontier.append(neighbor)
            frontier = next_frontier
        return reached

    def _empty_point_usable(self, point: Point) -> bool:
        for square_point in self.layer(PointState.SQUARE):
            if manhattan(point, square_point) <= self.square_clearance:
                return False
        for line_point in self.layer(PointState.LINE):
            if manhattan(point, line_point) <= self.line_separation:
                return False
        return True

    def _boundary_points(self, occupied: Set[Point]) -> Set[Point]:
        return {
            point
            for point in occupied
            if any(neighbor not in occupied for neighbor in four_neighbors(point))
        }

    def _require_in_bounds(self, point: Point) -> None:
        if not self.in_bounds(point):
            raise ValueError(f"point outside universe: {point}")

    def _add_state(self, point: Point, state: PointState) -> None:
        self._require_in_bounds(point)
        x, y = point
        self.matrix[y][x] |= state

    def _set_state(self, point: Point, state: PointState) -> None:
        self._require_in_bounds(point)
        x, y = point
        self.matrix[y][x] = state


@dataclass(frozen=True)
class HypotheticalWhiteConsequence:
    R_before: int
    R_after: int
    delta_R: int
    removed_white_points: int
    resulting_region_sizes: Tuple[int, ...]


@dataclass(frozen=True)
class HypotheticalTerminalCapacityConsequence:
    capacity_before: int
    capacity_after: Optional[int]
    delta_capacity: Optional[int]
    B_routable_after: bool
    finite_cut_after: bool
    minimum_cuts_before: Tuple[Tuple[Point, ...], ...]
    minimum_cuts_after: Tuple[Tuple[Point, ...], ...]


@dataclass(frozen=True)
class TerminalPair:
    pair_id: str
    start: Point
    goal: Point


@dataclass(frozen=True)
class MultiTerminalCapacityConsequence:
    candidate_points: Tuple[Point, ...]
    pair_consequences: Tuple[
        Tuple[str, HypotheticalTerminalCapacityConsequence],
        ...
    ]


@dataclass(frozen=True)
class AggregateTerminalCapacityConsequence:
    pair_count: int
    routable_after_count: int
    disconnected_after_count: int
    unchanged_capacity_pair_ids: Tuple[str, ...]
    reduced_positive_capacity_pair_ids: Tuple[str, ...]
    disconnected_pair_ids: Tuple[str, ...]
    total_capacity_before: int
    total_capacity_after: int
    total_delta_capacity: int


@dataclass
class _ResidualEdge:
    target: Tuple[Point, str]
    reverse_index: int
    capacity: int


class WhiteSpaceViewer:
    def __init__(self, universe: PointUniverse):
        self.universe = universe
        self._region_at: Dict[Point, int] = {}
        self._sizes: Dict[int, int] = {}
        self._build_regions()

    def _build_regions(self) -> None:
        white = set(self.universe.layer(PointState.WHITE))
        next_id = 1
        while white:
            start = min(white, key=lambda p: (p[1], p[0]))
            queue = deque([start])
            white.remove(start)
            size = 0
            while queue:
                point = queue.popleft()
                self._region_at[point] = next_id
                size += 1
                for neighbor in four_neighbors(point):
                    if neighbor in white:
                        white.remove(neighbor)
                        queue.append(neighbor)
            self._sizes[next_id] = size
            next_id += 1

    def white_region_count(self) -> int:
        return len(self._sizes)

    def white_region_size(self, region_id: int) -> int:
        return self._sizes[region_id]

    def white_region_at(self, x: int, y: int) -> Optional[int]:
        return self._region_at.get((x, y))

    def same_white_region(self, point_a: Point, point_b: Point) -> bool:
        region_a = self._region_at.get(point_a)
        return region_a is not None and region_a == self._region_at.get(point_b)

    def white_region_boundary(self, region_id: int) -> Set[Point]:
        return {
            point
            for point, point_region in self._region_at.items()
            if point_region == region_id and self._is_boundary_point(point)
        }

    def boundary_components(self, region_id: int) -> Dict[int, Set[Point]]:
        boundary = set(self.white_region_boundary(region_id))
        components: Dict[int, Set[Point]] = {}
        next_id = 1
        while boundary:
            start = min(boundary, key=lambda p: (p[1], p[0]))
            queue = deque([start])
            boundary.remove(start)
            component: Set[Point] = set()
            while queue:
                point = queue.popleft()
                component.add(point)
                for neighbor in four_neighbors(point):
                    if neighbor in boundary:
                        boundary.remove(neighbor)
                        queue.append(neighbor)
            components[next_id] = component
            next_id += 1
        return components

    def same_boundary_component(self, region_id: int, point_a: Point, point_b: Point) -> bool:
        for component in self.boundary_components(region_id).values():
            if point_a in component and point_b in component:
                return True
        return False

    def boundary_path(self, region_id: int, point_a: Point, point_b: Point) -> Optional[List[Point]]:
        boundary = self.white_region_boundary(region_id)
        if point_a not in boundary or point_b not in boundary:
            return None
        queue = deque([point_a])
        previous: Dict[Point, Optional[Point]] = {point_a: None}
        while queue:
            point = queue.popleft()
            if point == point_b:
                break
            for neighbor in sorted(four_neighbors(point), key=lambda p: (p[1], p[0])):
                if neighbor in boundary and neighbor not in previous:
                    previous[neighbor] = point
                    queue.append(neighbor)
        if point_b not in previous:
            return None
        path: List[Point] = []
        cursor: Optional[Point] = point_b
        while cursor is not None:
            path.append(cursor)
            cursor = previous[cursor]
        path.reverse()
        return path

    def hypothetical_white_consequence(self, candidate_points: Iterable[Point]) -> HypotheticalWhiteConsequence:
        points = tuple(candidate_points)
        self._validate_hypothetical_points(points)
        white_after = set(self._region_at) - set(points)
        sizes = self._white_component_sizes(white_after)
        r_before = self.white_region_count()
        r_after = len(sizes)
        return HypotheticalWhiteConsequence(
            R_before=r_before,
            R_after=r_after,
            delta_R=r_after - r_before,
            removed_white_points=len(points),
            resulting_region_sizes=tuple(sorted(sizes)),
        )

    def terminal_separating_cut_points(self, start: Point, goal: Point) -> Tuple[Point, ...]:
        self._validate_white_terminal(start, "start")
        self._validate_white_terminal(goal, "goal")
        if start == goal:
            return ()
        if not self.same_white_region(start, goal):
            raise ValueError("start and goal must be in the same WHITE region")

        white_points = set(self._region_at)
        discovery: Dict[Point, int] = {}
        low: Dict[Point, int] = {}
        cut_points: Set[Point] = set()
        time = 0

        def visit(point: Point, parent: Optional[Point]) -> bool:
            nonlocal time
            discovery[point] = time
            low[point] = time
            time += 1
            subtree_contains_goal = point == goal
            for neighbor in sorted(four_neighbors(point), key=lambda p: (p[1], p[0])):
                if neighbor not in white_points:
                    continue
                if neighbor == parent:
                    continue
                if neighbor not in discovery:
                    child_contains_goal = visit(neighbor, point)
                    subtree_contains_goal = subtree_contains_goal or child_contains_goal
                    low[point] = min(low[point], low[neighbor])
                    if (
                        point != start
                        and point != goal
                        and child_contains_goal
                        and low[neighbor] >= discovery[point]
                    ):
                        cut_points.add(point)
                else:
                    low[point] = min(low[point], discovery[neighbor])
            return subtree_contains_goal

        visit(start, None)
        return tuple(sorted(cut_points, key=lambda p: (p[1], p[0])))

    def terminal_minimum_vertex_cuts(self, start: Point, goal: Point) -> Tuple[Tuple[Point, ...], ...]:
        self._validate_white_terminal(start, "start")
        self._validate_white_terminal(goal, "goal")
        if start == goal:
            return ()
        if not self.same_white_region(start, goal):
            raise ValueError("start and goal must be in the same WHITE region")
        if manhattan(start, goal) == 1:
            return ()

        white_points = set(self._region_at)
        return self._terminal_minimum_vertex_cuts_for_white_points(start, goal, white_points)

    def hypothetical_terminal_capacity_consequence(
        self,
        candidate_points: Iterable[Point],
        *,
        B_start: Point,
        B_goal: Point,
    ) -> HypotheticalTerminalCapacityConsequence:
        candidate = tuple(candidate_points)
        self._validate_hypothetical_points(candidate)
        candidate_set = set(candidate)
        minimum_cuts_before = self._validate_terminal_capacity_pair(
            B_start,
            B_goal,
            candidate_set,
            start_label="B_start",
            goal_label="B_goal",
            distinct_message="B_start and B_goal must be distinct",
            same_region_message="B_start and B_goal must be in the same WHITE region",
            candidate_start_message="hypothetical candidate may not include B_start",
            candidate_goal_message="hypothetical candidate may not include B_goal",
        )
        return self._hypothetical_terminal_capacity_consequence_for_validated_pair(
            candidate_set,
            B_start,
            B_goal,
            minimum_cuts_before,
        )

    def hypothetical_multi_terminal_capacity_consequence(
        self,
        candidate_points: Iterable[Point],
        terminal_pairs: Iterable[TerminalPair],
    ) -> MultiTerminalCapacityConsequence:
        candidate = tuple(candidate_points)
        self._validate_hypothetical_points(candidate)
        candidate_set = set(candidate)

        pairs = tuple(terminal_pairs)
        if not pairs:
            raise ValueError("at least one TerminalPair is required")
        seen_pair_ids: Set[str] = set()
        validated: List[Tuple[TerminalPair, Tuple[Tuple[Point, ...], ...]]] = []
        for pair in pairs:
            if not isinstance(pair.pair_id, str):
                raise ValueError("pair_id must be a string")
            if not pair.pair_id:
                raise ValueError("pair_id must be non-empty")
            if pair.pair_id in seen_pair_ids:
                raise ValueError(f"duplicate pair_id: {pair.pair_id}")
            seen_pair_ids.add(pair.pair_id)
            minimum_cuts_before = self._validate_terminal_capacity_pair(
                pair.start,
                pair.goal,
                candidate_set,
                start_label="start",
                goal_label="goal",
                distinct_message="terminal pair endpoints must be distinct",
                same_region_message="terminal pair endpoints must be in the same WHITE region",
                candidate_start_message="hypothetical candidate may not include terminal pair start",
                candidate_goal_message="hypothetical candidate may not include terminal pair goal",
            )
            validated.append((pair, minimum_cuts_before))

        return MultiTerminalCapacityConsequence(
            candidate_points=candidate,
            pair_consequences=tuple(
                (
                    pair.pair_id,
                    self._hypothetical_terminal_capacity_consequence_for_validated_pair(
                        candidate_set,
                        pair.start,
                        pair.goal,
                        minimum_cuts_before,
                    ),
                )
                for pair, minimum_cuts_before in validated
            ),
        )

    def aggregate_terminal_capacity_consequence(
        self,
        consequence: MultiTerminalCapacityConsequence,
    ) -> AggregateTerminalCapacityConsequence:
        if not isinstance(consequence, MultiTerminalCapacityConsequence):
            raise ValueError("consequence must be a MultiTerminalCapacityConsequence")
        if not consequence.pair_consequences:
            raise ValueError("pair_consequences must be non-empty")

        seen_pair_ids: Set[str] = set()
        unchanged: List[str] = []
        reduced_positive: List[str] = []
        disconnected: List[str] = []
        total_before = 0
        total_after = 0
        total_delta = 0
        routable_count = 0
        disconnected_count = 0

        for pair_id, pair_consequence in consequence.pair_consequences:
            if not isinstance(pair_id, str):
                raise ValueError("pair_id must be a string")
            if not pair_id:
                raise ValueError("pair_id must be non-empty")
            if pair_id in seen_pair_ids:
                raise ValueError(f"duplicate pair_id: {pair_id}")
            seen_pair_ids.add(pair_id)
            if not isinstance(pair_consequence, HypotheticalTerminalCapacityConsequence):
                raise ValueError("pair consequence must be a HypotheticalTerminalCapacityConsequence")
            if not isinstance(pair_consequence.capacity_before, int) or pair_consequence.capacity_before <= 0:
                raise ValueError("capacity_before must be a positive integer")
            if pair_consequence.capacity_after is None:
                raise ValueError("capacity_after must be an integer")
            if not isinstance(pair_consequence.capacity_after, int) or pair_consequence.capacity_after < 0:
                raise ValueError("capacity_after must be a non-negative integer")
            if pair_consequence.delta_capacity is None:
                raise ValueError("delta_capacity must be an integer")
            if not isinstance(pair_consequence.delta_capacity, int):
                raise ValueError("delta_capacity must be an integer")
            if pair_consequence.delta_capacity != pair_consequence.capacity_after - pair_consequence.capacity_before:
                raise ValueError("delta_capacity must equal capacity_after - capacity_before")
            if not pair_consequence.minimum_cuts_before:
                raise ValueError("minimum_cuts_before must be non-empty")

            if not pair_consequence.B_routable_after:
                if pair_consequence.capacity_after != 0:
                    raise ValueError("disconnected consequence must have capacity_after == 0")
                if pair_consequence.finite_cut_after:
                    raise ValueError("disconnected consequence must have finite_cut_after == False")
                if pair_consequence.minimum_cuts_after:
                    raise ValueError("disconnected consequence must have empty minimum_cuts_after")
                disconnected.append(pair_id)
                disconnected_count += 1
            else:
                if pair_consequence.capacity_after == 0:
                    raise ValueError("routable consequence must have positive capacity_after")
                if not pair_consequence.finite_cut_after:
                    raise ValueError("routable consequence must have finite_cut_after == True")
                if not pair_consequence.minimum_cuts_after:
                    raise ValueError("routable consequence must have non-empty minimum_cuts_after")
                if pair_consequence.delta_capacity > 0:
                    raise ValueError("routable capacity increase is not an accepted aggregate class")
                if pair_consequence.delta_capacity == 0:
                    unchanged.append(pair_id)
                elif pair_consequence.delta_capacity < 0:
                    reduced_positive.append(pair_id)
                routable_count += 1

            total_before += pair_consequence.capacity_before
            total_after += pair_consequence.capacity_after
            total_delta += pair_consequence.delta_capacity

        pair_count = len(consequence.pair_consequences)
        if pair_count != routable_count + disconnected_count:
            raise ValueError("routable and disconnected counts must cover every pair")
        if pair_count != len(unchanged) + len(reduced_positive) + len(disconnected):
            raise ValueError("aggregate classes must cover every pair exactly once")
        if total_delta != total_after - total_before:
            raise ValueError("total_delta_capacity must equal total_capacity_after - total_capacity_before")

        return AggregateTerminalCapacityConsequence(
            pair_count=pair_count,
            routable_after_count=routable_count,
            disconnected_after_count=disconnected_count,
            unchanged_capacity_pair_ids=tuple(unchanged),
            reduced_positive_capacity_pair_ids=tuple(reduced_positive),
            disconnected_pair_ids=tuple(disconnected),
            total_capacity_before=total_before,
            total_capacity_after=total_after,
            total_delta_capacity=total_delta,
        )

    def _validate_terminal_capacity_pair(
        self,
        start: Point,
        goal: Point,
        candidate_set: Set[Point],
        *,
        start_label: str,
        goal_label: str,
        distinct_message: str,
        same_region_message: str,
        candidate_start_message: str,
        candidate_goal_message: str,
    ) -> Tuple[Tuple[Point, ...], ...]:
        self._validate_white_terminal(start, start_label)
        self._validate_white_terminal(goal, goal_label)
        if start == goal:
            raise ValueError(distinct_message)
        if not self.same_white_region(start, goal):
            raise ValueError(same_region_message)
        if start in candidate_set:
            raise ValueError(candidate_start_message)
        if goal in candidate_set:
            raise ValueError(candidate_goal_message)

        white_before = set(self._region_at)
        minimum_cuts_before = self._terminal_minimum_vertex_cuts_for_white_points(start, goal, white_before)
        if not minimum_cuts_before:
            raise ValueError("a finite non-terminal terminal cut is required")
        return minimum_cuts_before

    def _hypothetical_terminal_capacity_consequence_for_validated_pair(
        self,
        candidate_set: Set[Point],
        B_start: Point,
        B_goal: Point,
        minimum_cuts_before: Tuple[Tuple[Point, ...], ...],
    ) -> HypotheticalTerminalCapacityConsequence:
        capacity_before = len(minimum_cuts_before[0])

        white_before = set(self._region_at)
        white_after = white_before - candidate_set
        if not self._white_terminals_connected_after_removal(B_start, B_goal, white_before, candidate_set):
            return HypotheticalTerminalCapacityConsequence(
                capacity_before=capacity_before,
                capacity_after=0,
                delta_capacity=-capacity_before,
                B_routable_after=False,
                finite_cut_after=False,
                minimum_cuts_before=minimum_cuts_before,
                minimum_cuts_after=(),
            )

        minimum_cuts_after = self._terminal_minimum_vertex_cuts_for_white_points(B_start, B_goal, white_after)
        if not minimum_cuts_after:
            return HypotheticalTerminalCapacityConsequence(
                capacity_before=capacity_before,
                capacity_after=None,
                delta_capacity=None,
                B_routable_after=True,
                finite_cut_after=False,
                minimum_cuts_before=minimum_cuts_before,
                minimum_cuts_after=(),
            )
        capacity_after = len(minimum_cuts_after[0])
        return HypotheticalTerminalCapacityConsequence(
            capacity_before=capacity_before,
            capacity_after=capacity_after,
            delta_capacity=capacity_after - capacity_before,
            B_routable_after=True,
            finite_cut_after=True,
            minimum_cuts_before=minimum_cuts_before,
            minimum_cuts_after=minimum_cuts_after,
        )

    def _terminal_minimum_vertex_cuts_for_white_points(
        self,
        start: Point,
        goal: Point,
        white_points: Set[Point],
    ) -> Tuple[Tuple[Point, ...], ...]:
        nonterminal_points = tuple(sorted(white_points - {start, goal}, key=lambda p: (p[1], p[0])))
        infinite_capacity = len(nonterminal_points) + 1
        minimum_cardinality = self._terminal_vertex_cut_cardinality(
            start,
            goal,
            white_points,
            infinite_capacity,
        )
        if minimum_cardinality >= infinite_capacity:
            return ()

        cuts: List[Tuple[Point, ...]] = []
        for candidate in combinations(nonterminal_points, minimum_cardinality):
            removed = set(candidate)
            if not self._white_terminals_connected_after_removal(start, goal, white_points, removed):
                cuts.append(tuple(sorted(candidate, key=lambda p: (p[1], p[0]))))
        return tuple(sorted(cuts, key=lambda cut: tuple((point[1], point[0]) for point in cut)))

    def _terminal_vertex_cut_cardinality(
        self,
        start: Point,
        goal: Point,
        white_points: Set[Point],
        infinite_capacity: int,
    ) -> int:
        graph: Dict[Tuple[Point, str], List[_ResidualEdge]] = {}

        def add_edge(source: Tuple[Point, str], target: Tuple[Point, str], capacity: int) -> None:
            graph.setdefault(source, [])
            graph.setdefault(target, [])
            forward = _ResidualEdge(
                target=target,
                reverse_index=len(graph[target]),
                capacity=capacity,
            )
            reverse = _ResidualEdge(
                target=source,
                reverse_index=len(graph[source]),
                capacity=0,
            )
            graph[source].append(forward)
            graph[target].append(reverse)

        for point in white_points:
            capacity = infinite_capacity if point in (start, goal) else 1
            add_edge((point, "in"), (point, "out"), capacity)
        for point in white_points:
            for neighbor in four_neighbors(point):
                if neighbor in white_points:
                    add_edge((point, "out"), (neighbor, "in"), infinite_capacity)

        source = (start, "out")
        sink = (goal, "in")
        flow = 0
        while True:
            parent: Dict[Tuple[Point, str], Tuple[Tuple[Point, str], int]] = {}
            seen = {source}
            queue = deque([source])
            while queue and sink not in seen:
                node = queue.popleft()
                indexed_edges = sorted(
                    enumerate(graph[node]),
                    key=lambda item: (item[1].target[0][1], item[1].target[0][0], item[1].target[1], item[0]),
                )
                for edge_index, edge in indexed_edges:
                    if edge.capacity > 0 and edge.target not in seen:
                        seen.add(edge.target)
                        parent[edge.target] = (node, edge_index)
                        queue.append(edge.target)
            if sink not in seen:
                return flow

            path_capacity = infinite_capacity
            cursor = sink
            while cursor != source:
                previous, edge_index = parent[cursor]
                path_capacity = min(path_capacity, graph[previous][edge_index].capacity)
                cursor = previous

            cursor = sink
            while cursor != source:
                previous, edge_index = parent[cursor]
                reverse_index = graph[previous][edge_index].reverse_index
                graph[previous][edge_index].capacity -= path_capacity
                graph[cursor][reverse_index].capacity += path_capacity
                cursor = previous
            flow += path_capacity

    def _white_terminals_connected_after_removal(
        self,
        start: Point,
        goal: Point,
        white_points: Set[Point],
        removed: Set[Point],
    ) -> bool:
        if start in removed or goal in removed:
            return False
        remaining = set(white_points) - removed
        if start not in remaining or goal not in remaining:
            return False
        queue = deque([start])
        seen = {start}
        while queue:
            point = queue.popleft()
            if point == goal:
                return True
            for neighbor in four_neighbors(point):
                if neighbor in remaining and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        return False

    def _validate_white_terminal(self, point: Point, label: str) -> None:
        if not self.universe.in_bounds(point):
            raise ValueError(f"{label} point outside universe: {point}")
        if not (self.universe.point_state(point) & PointState.WHITE):
            raise ValueError(f"{label} point is not WHITE: {point}")

    def _validate_hypothetical_points(self, points: Tuple[Point, ...]) -> None:
        if not points:
            raise ValueError("hypothetical candidate is empty")
        seen: Set[Point] = set()
        for point in points:
            if point in seen:
                raise ValueError(f"hypothetical candidate repeats point {point}")
            seen.add(point)
            if not self.universe.in_bounds(point):
                raise ValueError(f"hypothetical candidate point outside universe: {point}")
            if not (self.universe.point_state(point) & PointState.WHITE):
                raise ValueError(f"hypothetical candidate point is not WHITE: {point}")
        for a, b in zip(points, points[1:]):
            if manhattan(a, b) != 1:
                raise ValueError("hypothetical candidate points must be four-neighbor adjacent")

    def _white_component_sizes(self, white_points: Set[Point]) -> List[int]:
        remaining = set(white_points)
        sizes: List[int] = []
        while remaining:
            start = min(remaining, key=lambda p: (p[1], p[0]))
            queue = deque([start])
            remaining.remove(start)
            size = 0
            while queue:
                point = queue.popleft()
                size += 1
                for neighbor in four_neighbors(point):
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        queue.append(neighbor)
            sizes.append(size)
        return sizes

    def _is_boundary_point(self, point: Point) -> bool:
        if point not in self._region_at:
            return False
        for neighbor in four_neighbors(point):
            if not self.universe.in_bounds(neighbor):
                return True
            if not (self.universe.point_state(neighbor) & PointState.WHITE):
                return True
        return False


def four_neighbors(point: Point) -> Tuple[Point, Point, Point, Point]:
    x, y = point
    return ((x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y))


def expanded_segment_points(line: Line) -> List[Set[Point]]:
    segments: List[Set[Point]] = []
    for a, b in zip(line.points, line.points[1:]):
        if a[0] == b[0]:
            x = a[0]
            y0, y1 = sorted((a[1], b[1]))
            segments.append({(x, y) for y in range(y0, y1 + 1)})
        else:
            y = a[1]
            x0, x1 = sorted((a[0], b[0]))
            segments.append({(x, y) for x in range(x0, x1 + 1)})
    return segments


def points_within_manhattan(point: Point, distance: int) -> Iterable[Point]:
    x, y = point
    for dx in range(-distance, distance + 1):
        for dy in range(-distance, distance + 1):
            if abs(dx) + abs(dy) <= distance:
                yield (x + dx, y + dy)


def manhattan(a: Point, b: Point) -> int:
    return abs(a[0] - b[0]) + abs(a[1] - b[1])
