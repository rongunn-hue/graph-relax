# Point Universe -> 2D router: what to try next (written 2026-09-20)

## Guiding idea (the user's, the reason this project exists)
The 2D router only knows where a wire starts and ends. It knows nothing about the
world in between. The Point Universe is a way to give an algorithm the global view
("from outside the plane"): what regions exist, where space is scarce, and what a
commit would do to the whole picture before it is made. Judge every idea by one
question: does it give the router knowledge of the global space it lacks?

## Already shown (real router, HRNG, fixed placement, subclass only)
- Router is blind to crossings (only parallel spacing). 27 crossings, order-insensitive (23-27).
- Hard "no line on line" + truthful fallback: 27 -> 8 (default order), best 5 with adaptive
  order; verify() and pin-level equivalence PASS. Cost: length +39%, 2-4x slower.
  Files: ../no_cross_router_experiment.py, ../no_cross_adaptive.py, ../net_order_experiment.py
- Forced-connection list = free crossing certificate (which nets have no crossing-free route).

## To try, in this order
1. PREDICT forced crossings topologically (flood fill / small max-flow on the router's
   free grid) instead of routing. Goal: a millisecond placement oracle vs 22-90s per route.
   Test: correlate predicted disconnected-pending-net count with the router's actual
   forced count over the ~40 placements logged by local_crossing_repair.py.
   Kill it if the correlation is weak.
2. MIN-CUT RESERVATION COST: cells in a pending pair's minimum vertex cut cost extra for
   other nets (protect the corridor a future net needs); re-run Dijkstra with choke cells
   penalised to get alternative paths. Replaces the noisy forced-first reordering
   (8,20,5,11,9). Limit: per-pair cuts cannot see two pending nets competing for one gap.
3. ROTATE/MIRROR as a placement move picked by the crossing certificate (the remaining 5
   crossings are op-amp feedback nets on U1/U2). Only if orientation is legal for that
   symbol (DIP notch orientation may forbid it; generic boxes / 2-pin parts likely fine).
4. Lower value: horizontal-run cell graph for speed (spacing/congestion costs depend on
   placed wires, so compression is awkward).

## Also open
- Finite crossing penalty instead of hard block (may recover length).
- Run IAMP and mono741, not just HRNG.
- Two point_universe tests fail in original and copy (future_routability fixture,
  capacity_aware_comparison finds no fixture).

## RESULT of idea 1 (2026-09-20): KILLED as specified
Predictor = NoCrossRouter with rt.cheap=True (BFS paths, same hard occupancy) vs exact (Dijkstra),
24 placements (baseline + stride-3 sample of the 74 legal territory swaps). Files:
forced_predictor_test.py, static_proxy_check.py, forced_predictor_results.json.
- cheap runs 1.9s vs exact 51.8s (27x faster) but Spearman cheap_forced vs exact_forced = 0.19;
  ranks pairs correctly 40%, wrongly 21%, ties 40%. Cheap BFS over-forces (10 vs 4): its staircase
  paths block space differently, and it misses the swaps that actually matter.
- Free static proxies also useless: HPWL -0.08, net-bbox overlap 0.24. The damaging swaps
  (C13<->C2 forced 9, C13<->R4 8, C14<->* 6) have equal/LOWER HPWL: they break the decoupling caps'
  pin-aligned placement next to U1 pins 8/4, which centre-based metrics cannot see.
- KEY FACT: the exact hard-mode landscape over single swaps is FLAT. 18 of 24 swaps -> identical
  forced=4 / crossings 8; the rest are worse (up to 9 / 16); none improved on baseline forced=4.
  So swap-based local repair has nothing to gain here; the 4 forced connections (U1A_OUT, U2_MINUS,
  U2_PLUS) are structural (op-amp feedback vs fixed DIP pin sides), not placement-order noise.
  One exact run timed out (R3<->U1): pathological router search again.
=> Evidence now favours idea 3 (orientation/pin-side freedom) over idea 2. A better predictor, if
   ever wanted, must be PIN-aware (real alloc() pin coordinates), not centre-based.

## 2026-09-20 (later) — the "child's view" is quantifiable; first real signal
User: the viewer does not see holistically; a child sees how to move things to untangle; nobody knows how to quantify it.
What worked: the RUBBER-BAND PICTURE (connected things as straight bands pin->junction, real pin offsets, decoupling
rail wires included) and its crossing count, computed from geometry only (rubberband.py, rubberband_check.py).
- Predictors of the exact hard-mode router's forced crossings over 24 placements (Spearman): cheap-BFS sim 0.19,
  HPWL -0.08, net-bbox overlap 0.24, pin-direction mismatch 0.21, RUBBER-BAND CROSSINGS 0.73 (0.60 vs crossings).
- Child's move made computable (rubberband_untangle.py): for each device in band crossings, scan the whole table
  (grid 30, radius 400) for the spot where ITS bands cross least; greedy. 21 -> 0 band crossings in 13 moves, 7 seconds.
  Real router once on the result: blind (existing) router 27 -> 17 crossings, length 9184 -> 8984 (better);
  hard-no-cross router 8 -> 14, forced 4 -> 7 (worse). verify() and pin-level equivalence PASS.
- Why not zero: bands ignore BODIES and PIN SIDES, so bands that are straight on paper still wrap around boxes. Also no
  nearness term: U1 was carried 400 units and C14 (its -9V decoupling cap) was stranded far from U1 pin 4.
- Topology (sound): with real template pin orders (U1, U2, Q1) and free junction rotations a crossing-free embedding
  EXISTS for all mirror choices (rotation_genus.py; signal graph has only 5 independent cycles). So HRNG's crossings are
  not topologically forced; my earlier "op-amp feedback is structural" guess was wrong.
  (actual_rotation_check.py is INCONCLUSIVE: could not read 6-7 of 13 junctions from drawings; ignore its verdict.)
NEXT: add bodies-as-obstacles + pin-side (band must leave the pin outward) + nearness/SOFT_NEARNESS + displacement
cost to the picture; keep block structure; re-test on the 24 placements (does correlation hold/improve?) then route once.
Also try IAMP/mono741 for generality.
