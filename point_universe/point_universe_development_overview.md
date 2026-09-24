# Point Universe Development Overview

## Scope

This isolated experiment lives under:

`/home/rgunn/Chat-Projects/graph-relax/experiments/point_universe/`

It is separate from the existing schematic router, placement code, HRNG circuit work, and production modules. The purpose was to build a tiny deterministic mathematical routing universe step by step, with tests proving each behavior before moving to the next layer.

No HRNG routing, production-router modification, placement modification, optimization, backtracking, or global routing policy was added in this experiment.

## Core model established

The experiment defines a finite integer-grid point universe. The authoritative geometry is the point matrix, not rendered graphics or inferred continuous geometry.

The universe currently represents:

- `SQUARE` occupied objects
- `LINE` occupied paths
- `ATTACHMENT` points
- usable `WHITE` space
- empty but unusable `UNUSABLE` space

Squares, lines, line occupancy, connectivity, same-net/unrelated-line behavior, white/unusable classification, and validation are handled by `PointUniverse` in `universe.py`.

## Viewer capabilities established

`WhiteSpaceViewer` now perceives usable WHITE space directly from the authoritative matrix.

It supports:

- four-neighbor WHITE connected components
- deterministic region IDs and region sizes
- region membership queries
- white-space boundary extraction
- boundary components
- deterministic boundary paths
- hypothetical pre-commit WHITE-space consequence analysis

The key perception method added was:

`WhiteSpaceViewer.hypothetical_white_consequence(candidate_points)`

It accepts an already-expanded ordered sequence of WHITE grid points and reports what would happen if those points ceased to be WHITE, without mutating the universe.

It returns:

- `R_before`
- `R_after`
- `delta_R`
- `removed_white_points`
- sorted `resulting_region_sizes`

This is perception only. It does not route, rank, choose, or commit anything.

## Boundary and candidate-consequence result

A same-connection experiment demonstrated that two legal candidate paths between the same endpoints can have different WHITE-space consequences.

For the established candidate-consequence geometry:

- boundary-following candidate: `R_before=1`, `R_after=1`, `delta_R=0`
- interior-cut candidate: `R_before=1`, `R_after=2`, `delta_R=+1`
- after the interior candidate, point `(2,4)` remains WHITE but becomes an isolated WHITE region of size `1`

This proved the viewer can perceive a topological difference before commitment.

## First PointRouter

`router.py` introduced `PointRouter`.

The router:

- enumerates legal simple paths through current WHITE points
- uses four-neighbor moves only
- uses fixed neighbor order: UP, LEFT, RIGHT, DOWN
- is deterministic
- is bounded by `max_paths` complete returned paths
- attaches `WhiteSpaceViewer.hypothetical_white_consequence(...)` to each returned candidate

The router does not:

- commit lines
- select a route
- rank candidates
- optimize length or bends
- use Dijkstra, A*, NetworkX, randomness, or the production schematic router

The returned object is `RouteCandidate(points, consequence)`.

## First RouteSelector

`route_selector.py` introduced `RouteSelector`.

The selector receives already-generated `RouteCandidate` objects and chooses one by exactly this lexicographic key:

```python
(
    candidate.consequence.delta_R,
    len(candidate.points),
    candidate.points,
)
```

Lowest tuple wins.

This means:

1. smaller `delta_R` wins
2. if tied, fewer occupied path points wins
3. if tied, lexicographically smaller complete point tuple wins

The selector does not generate routes, call the viewer, call the router, commit lines, mutate candidates, or create replacement candidates.

The verified established case selected a longer 31-point candidate with `delta_R=0` over a shorter 29-point candidate with `delta_R=+1`, proving that preserving WHITE-space connectivity precedes path length.

## First RouteCommitter

`route_committer.py` introduced `RouteCommitter`.

It commits one already-selected candidate by constructing exactly:

```python
Line(line_id, net_id, candidate.points)
```

and passing it to:

```python
PointUniverse.add_line(line)
```

The committer does not alter, simplify, reverse, reorder, or repair candidate points. It does not catch `add_line` failures. Existing `PointUniverse` legality remains authoritative.

The complete one-connection cycle was verified:

```text
Universe
  -> PointRouter generates candidates
  -> WhiteSpaceViewer attaches hypothetical consequences
  -> RouteSelector selects one candidate
  -> RouteCommitter commits the exact selected path
  -> WhiteSpaceViewer observes the updated universe
```

The tests confirmed that the pre-commit hypothetical consequence matched the post-commit actual WHITE-region consequence.

## Final experiment: two-connection sequential mechanism

The last experiment added:

`test_two_connection_sequence.py`

This was a mechanism test, not a new router. It composed the existing pieces:

- `PointUniverse`
- `WhiteSpaceViewer`
- `PointRouter`
- `RouteSelector`
- `RouteCommitter`

The primary geometry was:

```text
PointUniverse(7, 5, line_separation=0, square_clearance=0)
Square("CENTER", 2, 1, 2, 2)

A_start = (0,2)
A_goal  = (6,2)

B_start = (0,4)
B_goal  = (6,4)
```

The test established the first two-connection sequential state transition:

```text
Route connection A against the initial universe
  -> select A
  -> commit A
  -> the authoritative universe changes
  -> route connection B against that same mutated universe
  -> select B
  -> commit B
  -> observe final universe
```

The critical proposition was that B must see the consequences of committed A. The test directly checked:

- A exists in `universe.lines` before B candidate generation
- B's router uses the same universe object
- B candidate `R_before` equals the post-A WHITE-region count
- B candidates do not occupy A's line points
- A survives B routing and commit unchanged

The test also included a blocked second-connection case:

```text
PointUniverse(3, 3, line_separation=0, square_clearance=0)
A path = ((1,0), (1,1), (1,2))
B_start = (0,1)
B_goal  = (2,1)
```

After A is committed, B returns no candidates. The test verifies that this is reported truthfully: no B route is invented, A is not removed, and A is not rerouted.

## Verified final test state

The final approved execution report for this stage was written to:

`/mnt/ZYXEL/DLINK/dump/two_connection_sequence_execution.txt`

It reported:

- tests collected: `40`
- tests passed: `40`
- tests failed: `0`
- errors: `0`
- final result: `PASS`

The full test suite at that point included all previous point-universe tests plus the final two-connection sequence test.

## What this experiment proves

The point-universe prototype now proves a minimal sequential pipeline:

1. represent exact discrete geometry in an authoritative point matrix
2. perceive connected WHITE space and boundary structure
3. evaluate hypothetical candidate consequences without mutation
4. enumerate deterministic legal candidate paths
5. select a candidate by a fixed consequence-first law
6. commit the exact selected path through existing universe validation
7. route the next connection against the mutated authoritative universe
8. truthfully report when a later connection has no legal candidate

## What this experiment does not claim

This does not yet prove:

- global optimality
- that greedy `delta_R` selection is always sufficient
- that boundary-following always wins
- that all future routes remain possible
- that this solves HRNG
- that this replaces the production schematic router

No multi-connection optimization, lookahead, backtracking, route ordering policy, HRNG integration, or production routing change was implemented.

## Main source files

- `universe.py` — authoritative point universe and viewer
- `router.py` — deterministic candidate-path enumerator
- `route_selector.py` — deterministic candidate selector
- `route_committer.py` — exact selected-candidate commit operation
- `test_two_connection_sequence.py` — final sequential two-connection mechanism test

## Current status

The isolated point-universe development chain is internally consistent and test-verified through the first two-connection sequential mechanism.
