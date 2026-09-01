#!/usr/bin/env python3
"""Build the requested single-file archive of graph-relax experiments."""

from pathlib import Path
from datetime import datetime, timezone


ROOT = Path(__file__).parent
OUTPUT = ROOT / "output"
DESTINATION = ROOT / "graph-relax_experiments_1_through_4R.txt"
MANIFEST = ROOT / "graph-relax_checkpoint_manifest.sha256"


COMMON = """COMMON PROJECT CONSTRAINTS

Project: graph-relax
Authoritative input: test_data/iamp.circuit (IAMP LM1458)
Input interpretation: component REF DEVICE_ID declares an anonymous object REF;
DEVICE_ID has no semantics.  A net line connects named REF.PIN terminals.

Unless an experiment explicitly changed one variable, the following stayed frozen:
- no KiCad and no deterministic_schematic imports;
- no device, symbol, electrical-function, power-net, signal-flow, or block knowledge;
- components are identical movable 60 x 40 rectangles with deterministic generic
  perimeter terminals;
- no schematic-specific aesthetics or Manhattan routing;
- deterministic ordering and fixed PRNG seeds;
- prior experiment artifacts and behavior are immutable;
- weights: squared length 0.02, overlap 100000, clearance 5000,
  crossing 2000, crowding 20;
- ordinary Experiment 2 annealing: seed 314159, temperature 25000 to 10
  geometrically, movement 140 to 0.5 geometrically, 6000 iterations,
  Metropolis acceptance, best-state retention, followed by the frozen downhill
  coordinate-descent quench;
- first predetermined configuration only: no tuning after viewing results.

Verified IAMP input: 19 components, 15 named nets, 44 terminal incidences,
and three singleton nets (R1_W, R4_W, RA_W).
"""


EXPERIMENTS = [
    ("Experiment 1", """PURPOSE AND INSTRUCTIONS
Test whether an electrical connectivity graph can organize itself into a useful
2-D layout by deterministic physical relaxation.  Multi-terminal nets use a
derived arithmetic-centroid star; singleton nets exert no force.  Optimize only
connection length, body overlap, minimum body separation, net crossings, and
unrelated net/body or net/net crowding.  Use deterministic generic initial
terminal placement and downhill coordinate descent, with fixed-seed perturbation
only if used.  Stop when no improving move exists at the minimum step or the
predetermined sweep limit is reached.  Produce initial/final coordinate JSON,
initial.svg, final.svg with identical rendering rules, metrics and complete
energy breakdown.  Run unchanged input ten times and require byte-identical
final coordinate JSON.  Before implementation the input, centroid-star rule,
energy terms, movement, noise, and stopping criterion were stated and approved.
The corrected input gate required exactly 15 nets and three named singleton nets.
""", ["metrics.json"]),

    ("Experiment 2", """PURPOSE AND INSTRUCTIONS
Freeze Experiment 1 and change only the search mechanism.  Replace deterministic
downhill-only relaxation with deterministic simulated annealing using seed
314159, conventional temperature-dependent Metropolis uphill acceptance,
temperature 25000 to 10, movement 140 to 0.5, geometric schedules, and 6000
iterations.  Components remain the only movable objects and multi-terminal nets
remain derived arithmetic centroids.  Retain the best annealing state and quench
it using the exact Experiment 1 coordinate descent.  Preserve Experiment 1
artifacts and create annealed SVG/coordinates/metrics.  Compare initial,
Experiment 1, and Experiment 2; report move counts, best energy, final energy,
quench sweeps, and ten-run hash determinism.  No parameter or weight tuning.
""", ["annealed_metrics.json"]),

    ("Experiment 3", """PURPOSE AND INSTRUCTIONS
Test a literal physical model while freezing weights and the Experiment 2
annealing schedule.  Two-terminal nets become direct terminal-to-terminal
segments.  Every 3+-terminal net receives exactly one movable, massless junction
initialized at its terminal arithmetic centroid; singleton nets remain metadata.
Annealing chooses deterministically from components and junctions; component
moves translate terminals and junction moves alter only the junction coordinate.
Crossing/crowding calculations use the actual segments.  Retain the best state
and quench both components and junctions.  Render visible junction dots, create
separate junction artifacts, compare Experiments 1-3, and require ten identical
coordinate hashes.  Do not add hard constraints, routing, or heuristics.
""", ["junction_annealed_metrics.json"]),

    ("Experiment 4A", """PURPOSE AND INSTRUCTIONS
Test degree-seeded initialization only, using the Experiment 2 centroid-star
physical model and unchanged optimizer.  Degree is the number of distinct
non-singleton nets incident to a component.  Sort by descending degree then
lexical REF.  Put a unique maximum-degree component at the exact center; tied
maxima go symmetrically on the smallest practical ring.  Place lower-degree bands
on progressively larger deterministic rings, using spacing derived once from
component dimensions and clearance, with no initial overlaps.  Once annealing
starts, release every component completely: no centrality force, pinning, or
restraint.  Save degree-seeded initial/final SVGs, coordinates, metrics, the full
degree table, annealing/quench statistics, comparisons, and ten-run hashes.
""", ["degree_seeded_metrics.json"]),

    ("Experiment 4B", """PURPOSE AND INSTRUCTIONS
Test degree-seeded 3-D relaxation followed by a fixed, deterministic search for
the best raw 2-D orthographic projection.  Begin from the exact 4A XY rings with
z=0.  Components translate in x/y/z without rotation; terminals share component
z and multi-terminal centroids are 3-D arithmetic means.  Use 3-D Euclidean
length and deterministic 3-D body/clearance/crowding generalizations.  Do not use
projected 2-D crossings during 3-D relaxation; penalize only actual 3-D segment
proximity using exact segment minimum distance.  Use a mechanically derived
symmetric ZMAX, no flattening penalty, the unchanged seed/schedules/6000 moves,
and a 3-D coordinate-descent quench.  Freeze the resulting 3-D configuration.
Search a predetermined spherical set of 128 viewing directions and 8 roll angles
(1024 projections), with no adaptive refinement.  Rigidly project component
centers, recreate ordinary axis-aligned 2-D rectangles and terminals, recompute
centroids, and score with the exact 4A objective.  Select by minimum objective and
enumeration-order tie break.  Do not perform a final 2-D quench.  Save 3-D state,
canonical XY/XZ/YZ SVGs, search record, best projection, comparisons, and ten-run
hashes.
""", ["degree_seeded_3d_metrics.json", "3d_projection_search.json",
             "best_3d_projection_metrics.json"]),

    ("Experiment 4C", """PURPOSE AND INSTRUCTIONS
Change only initialization to a topology-derived closeness organization.  Build
the unweighted component graph: every non-singleton electrical net makes all
participating components pairwise adjacent; duplicate adjacencies collapse and
singletons are excluded.  Compute all-pairs shortest paths, degree, distance sum,
and conventional closeness.  Select ROOT by highest closeness, lowest distance
sum, highest degree, lexical REF; place ROOT at (0,0).  Assign shortest-path
layers.  Give non-root components deterministic provisional ring positions, then
perform synchronous topology-only barycentric initialization toward graph
neighbors while keeping components in deterministic layer radial bands.  Resolve
overlap mechanically without consulting crossings, length, crowding, or final
energy.  Stop at convergence or the predetermined initialization sweep limit.
Save and score the frozen initial state, then release ROOT/layers and run the
unchanged Experiment 2 annealer/quench.  Save graph/centrality/layer data,
initial/final SVGs, coordinates, metrics, comparisons to 4A, and ten-run hashes.
""", ["closeness_graph.json", "closeness_metrics.json"]),

    ("Experiment 4D", """PURPOSE AND INSTRUCTIONS
Test neighborhood-affinity initialization only.  Reuse the 4C component graph
and closeness root (calculated, expected U1).  For non-root A,B define direct=1
when adjacent; shared=number of common neighbors excluding ROOT; affinity is
direct+shared.  Start non-root components in a deterministic outer arrangement
and synchronously move each toward the affinity-weighted barycenter of other
non-root components.  ROOT is excluded from the barycenter.  Preserve one weak,
mechanically derived radial band and resolve overlaps by deterministic minimum
displacement.  Initialization may use only graph, affinity, geometry, clearance,
and tie-breaking—not the final objective.  Use a fixed 100-sweep cap and
convergence test.  Then release root, radial, and affinity constraints and run the
unchanged 4C annealer/quench.  Save affinity graph, initial/final SVGs,
coordinates, metrics, comparisons, and ten-run hashes.
""", ["affinity_graph.json", "affinity_metrics.json"]),

    ("Experiment 4E", """PURPOSE AND INSTRUCTIONS
Start from the exact stored 4C closeness_initial coordinates; do not regenerate
them.  Verify 23 crossings and zero overlap.  Preprocess horizontally only: all Y
coordinates fixed, X slides and swaps only, overlap forbidden.  The sole search
objective is lexicographic (crossings, ordinary connection length); ordinary
energy is measurement only.  Initial X step is width+clearance=72 and minimum is
72/256=0.28125.  At each state enumerate every legal +/- slide and every legal X
swap whose vertical spans overlap, all from the same state; apply exactly the
best improvement, with slide/swap and lexical/direction tie breaks.  Halve step
when stuck and stop at the minimum step.  Log every move.  Release fixed Y and all
preprocessing knowledge, then run the exact 4C annealer/quench.  Save horizontal
and final SVG/coordinates/metrics/moves and require ten identical hashes.
""", ["closeness_horizontal_metrics.json", "closeness_horizontal_moves.json",
             "closeness_horizontal_annealed_metrics.json"]),

    ("Experiment 4F", """PURPOSE AND INSTRUCTIONS
Start again from exact 4C closeness_initial, not 4E output.  Alternate complete
horizontal and vertical one-dimensional best-improvement passes.  H passes are
exactly 4E.  V passes are the transpose: fixed X, +/-Y slides, Y swaps only for
components with overlapping horizontal spans, initial step height+clearance=52,
minimum 52/256=0.203125.  At every instant only one axis may move; no diagonal or
general XY candidate.  A pass continues through all step halvings.  Alternate
H,V until the latest pass on both axes accepts zero moves, with a fixed 20-pass
safety bound.  Record every move and every complete pass state, metrics, and
crossing sequence.  Release all axis restrictions and run the unchanged 4C/4E
annealer/quench.  Save preprocessing/final artifacts and require ten identical
preprocessed and final hashes.
""", ["closeness_hv_metrics.json", "closeness_hv_passes.json",
             "closeness_hv_moves.json", "closeness_hv_annealed_metrics.json"]),

    ("Experiment 4G", """PURPOSE AND INSTRUCTIONS
Perform diagnosis only on the exact frozen 4F final coordinates.  Reconstruct the
same components, terminals, centroid segments, and strict proper-crossing
predicate; verify seven crossings and zero overlap.  Enumerate crossings X001...
by intersection Y, X, nets, and endpoint identities.  Report exact intersections,
segments, anchors, complete net memberships, relevant component unions, and
locality boxes.  For each crossing independently, test hypothetical single-axis
slides and eligible swaps among relevant components using the frozen 4E/4F step
schedules, with no overlaps and no committed movement.  Rank target-eliminating
candidates by resulting crossings, newly introduced crossings, length,
displacement, and deterministic ties.  Classify X, Y, XY, or STRUCTURAL under the
specified rules.  Build the 7x7 E/P/R/N-A interaction matrix and component/net
involvement summaries.  Save JSON/text and an unchanged-geometry annotated SVG.
Run ten identical diagnoses.  Never apply a candidate, anneal, quench, or route.
""", ["residual_crossings.json", "residual_crossing_interactions.json",
             "residual_crossing_report.txt"]),

    ("Experiment 4H", """PURPOSE AND INSTRUCTIONS
Freeze every 4F component coordinate and allow only one massless junction for
each 3+-terminal net to move.  Singleton nets have no optimization; two-terminal
nets use direct segments.  Initialize seven movable junctions at exact arithmetic
centroids.  Optimize lexicographically (crossings, total connection length) only;
weighted energy, clearance, and crowding are measurements.  Axis moves are -X,
+X,-Y,+Y in lexical net order.  Initial step is max(width,height)+clearance=72,
minimum 72/256=0.28125.  Per-net bounds are its terminal box expanded by 72.
Junctions may not lie inside unrelated component bodies.  Use deterministic
best-improvement and step halving.  Save unchanged component coordinates,
initial/final junctions, bounds, move log, metrics, and diagnostic SVG.  Compare
X001-X007 identities and new net-pair crossings.  Stop before any annealing or
joint quench and require ten identical coordinate/move-log hashes.
""", ["junction_untangled_coordinates.json", "junction_untangled_moves.json",
             "junction_untangled_metrics.json"]),

    ("Experiment 4P", """PURPOSE AND INSTRUCTIONS
Certify exactly where planarity is lost, using source connectivity only and no
geometry optimization.  A is the bipartite component/net incidence graph, with
singleton incidences retained.  B explicitly inserts terminal vertices and
component cores, using component-terminal-net for every incidence.  Test A and B
with an established exact planarity implementation, save cyclic orders, and
independently verify graph equality, embedding structure, faces, and Euler
characteristic; save Kuratowski certificates if nonplanar.  C constrains the
physical cyclic order of terminals around every component disk; a conventional
unconstrained test is insufficient.  The exact reduction adds immutable rim
cycles to component terminal stars.  D, conditional on C nonplanarity, tests the
topologically unique current/180 and mirror-reversed orders.  If A were
nonplanar, search incidence-edge planarization through k=3; it was planar, so that
search was not applicable.  Record NetworkX version/function, exact graphs,
embeddings/certificates, distinctions among planarity variants, tests, and ten
identical hashes.
""", ["planarity_report.txt", "planarity_report.json",
             "planarity_A_embedding.json", "planarity_B_embedding.json",
             "planarity_C_nonplanarity_certificate.json", "planarity_D_result.json"]),

    ("Experiment 4Q", """PURPOSE AND INSTRUCTIONS
Starting from the exact mechanically reconstructed 4P wheel-augmented C graph,
find the exact minimum number of real electrical terminal-to-net incidences that
must be exempted from the primary constrained planar embedding.  Exemptions
remain electrically connected metadata for special graphical treatment.  Only
terminal-net electrical-incidence edges are removable; component-terminal
spokes, wheel/rim constraints, vertices, nets, and terminals are immutable.
Report all 44 candidates.  Exhaustively test combinations in increasing k,
bounded at k<=4; at the first successful k finish all combinations and do not
search k+1.  Exact NetworkX planarity tests are authoritative; optional
Kuratowski pruning must be sound (the completed search used no pruning).  Save
all minimum sets, verified embeddings, frequency/mandatory incidence-component-
net analysis, original and subsequent obstruction certificates, and the exact
search proof.  Attempt exact constrained reinsertion crossing cost only if an
established local implementation exists; otherwise report not computed.  Do not
equate exemption cardinality with crossing number and do not alter electrical
connectivity.  Run all prior tests and ten byte-identical complete diagnostics.
""", ["constrained_planarization_candidates.json",
             "constrained_planarization_minimum_sets.json",
             "constrained_planarization_frequency.json",
             "constrained_planarization_obstructions.json",
             "constrained_planarization_report.txt",
             "constrained_planarization_report.json",
             "constrained_planarization_embedding_SET001.json",
             "constrained_planarization_embedding_SET002.json",
             "constrained_planarization_embedding_SET003.json",
             "constrained_planarization_embedding_SET004.json"]),

    ("Experiment 4R", """PURPOSE AND INSTRUCTIONS
Attempt to turn each exact saved 4Q constrained-planar embedding into a compact
2-D physical realization, with topology authoritative before geometry.  Load and
verify SET001-SET004 exactly; do not regenerate their vertex sets, edge sets, or
rotation systems.  Replace mathematical wheel gadgets by ordinary anonymous
rectangles and fixed-order terminals.  Layer 0 must preserve the exact embedding,
have zero primary crossings, no body intersections, no false junctions, and no
component overlap.  Straight edges should be attempted first and deterministic
polylines used only where required.  Compact lexicographically while treating
zero crossings as a hard invariant, without annealing.  Only after Layer 0 is
frozen, restore the two exempt incidences as visibly distinct Layer 1 geometry;
Layer 1 may cross Layer 0 but may not enter component bodies.  Stage A must
realize the exact saved embedding.  Stage B may examine other embeddings only
after Stage A succeeds, and must say NOT EXHAUSTIVELY SEARCHED if complete exact
enumeration is unavailable.  Independently reconstruct and validate all geometry
and electrical incidence.  Compare all four sets using the predetermined ranking.
Failure is valid and must not be repaired by terminal reordering, component
mirroring, altered exemptions, a Layer-0 crossing, or a changed embedding.

Stage A input validation found that the saved wheel embeddings independently
reverse selected component terminal cycles: SET001-SET003 reverse J_PWR and U1;
SET004 reverses J_PWR, J_R4, and U1, while other multi-terminal components retain
current order.  A global reflection cannot repair this mixed chirality and
individual mirroring is forbidden.  Consequently all four realizations stopped
before generating Layer-0 or Layer-1 connections.  The failure artifacts use a
deterministic diagnostic grid and make no claim of geometric realization.
""", ["planar_realization_comparison.json", "planar_realization_report.txt",
             "planar_realization_SET001_coordinates.json",
             "planar_realization_SET001_geometry.json",
             "planar_realization_SET001_metrics.json",
             "planar_realization_SET001_validation.json",
             "planar_realization_SET002_coordinates.json",
             "planar_realization_SET002_geometry.json",
             "planar_realization_SET002_metrics.json",
             "planar_realization_SET002_validation.json",
             "planar_realization_SET003_coordinates.json",
             "planar_realization_SET003_geometry.json",
             "planar_realization_SET003_metrics.json",
             "planar_realization_SET003_validation.json",
             "planar_realization_SET004_coordinates.json",
             "planar_realization_SET004_geometry.json",
             "planar_realization_SET004_metrics.json",
             "planar_realization_SET004_validation.json"]),
]


def main():
    parts = ["GRAPH-RELAX EXPERIMENT ARCHIVE", "=" * 80,
             "Generated: " + datetime.now(timezone.utc).isoformat(),
             "Workspace: " + str(ROOT), "",
             "This is a consolidated, faithful instruction record followed by the",
             "complete saved textual/machine-readable results for every experiment.",
             "The instruction sections consolidate the governing user directives;",
             "the result artifacts are embedded verbatim.", "", COMMON]
    for title, instructions, files in EXPERIMENTS:
        parts.extend(["", "=" * 80, title.upper(), "=" * 80, "", instructions.strip(),
                      "", "COMPLETED RESULT ARTIFACTS (VERBATIM)"])
        for name in files:
            path = OUTPUT / name
            parts.extend(["", "-" * 80, f"FILE: output/{name}", "-" * 80])
            if path.exists():
                parts.append(path.read_text(encoding="utf-8").rstrip())
            else:
                parts.append("[MISSING ARTIFACT]")
    parts.extend(["", "=" * 80, "END OF ARCHIVE", "=" * 80, ""])
    DESTINATION.write_text("\n".join(parts), encoding="utf-8")
    import hashlib
    checkpoint_files = ([ROOT / "test_data" / "iamp.circuit", ROOT / "graph_relax.py"] +
                        sorted(ROOT.glob("experiment_*.py")) + sorted(ROOT.glob("test_*.py")) +
                        sorted(path for path in OUTPUT.iterdir() if path.is_file()) + [DESTINATION])
    unique_files = sorted(set(checkpoint_files), key=lambda path: str(path.relative_to(ROOT)))
    manifest_lines = []
    for path in unique_files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_lines.append(f"{digest}  {path.relative_to(ROOT)}")
    MANIFEST.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print(DESTINATION)
    print(MANIFEST)


if __name__ == "__main__":
    main()
