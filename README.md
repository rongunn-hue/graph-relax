# Graph-Relax

Graph-Relax is experimental research into automatically generating electrically
correct, human-understandable schematic layouts from circuit connectivity.

This repository preserves exploratory placement research, including graph
construction from circuit connectivity, graph planarity experiments, planar
embeddings, connectivity- and nearness-based placement, crossing reduction,
compaction, radial and tent-pole experiments, circular-envelope measurement,
and circular-envelope bulge correction.

The current architectural direction is:

```text
circuit/netlist
    -> abstract graph
    -> determine planar/nonplanar
    -> analyze connectivity/nearness
    -> intelligent initial placement
    -> geometric relaxation/compaction
    -> eventual schematic representation
```

This is research and prototype code, not a finished schematic generator.

## Graph placement invariant

Every new graph coordinate generator must pass its raw coordinates through
`placement_geometry.enforce_minimum_node_spacing` before generic geometry
validation and serialization. The shared `MIN_NODE_DISTANCE` is enforced by a
single uniform scale about the placement anchor, preserving angles, relative
distance ratios, aspect ratio, cyclic orientation, and crossing topology.
Coincident distinct vertices are a placement failure because uniform scaling
cannot separate them. Vertex-on-unrelated-edge cases are reported separately;
minimum node spacing does not prevent collinearity.
