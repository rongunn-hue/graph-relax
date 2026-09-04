# Local de-crowding checkpoint

**THIS CHECKPOINT IS BEFORE GLOBAL DE-CROWDING EXPERIMENTS.**

Pure nearness MDS produced a useful initial node placement. Direct analytic endpoint rotation then converted the IAMP test drawing into a zero-crossing straight-line drawing while preserving connectivity.

Crowding detection was generalized around a common local spacing scale, `D`, for node-node, node-edge, unrelated edge-edge, and shared-endpoint edge-overlay relationships. Deterministic direct local corrections materially improved the drawing and reduced measured crowding. However, sequential local correction is fundamentally insufficient: correcting one relation can displace congestion elsewhere. Protecting already-satisfied relations prevents some regressions, but does not solve that underlying problem.

The next research direction is to **keep the local crowding measurements, but combine their demands into a coordinated global node displacement** while preserving the planar zero-crossing embedding.

The local work must be retained because it supplies the detector, local `d` and `D` measurements, closest-point geometry, separation directions, validity checks, and deterministic candidate mathematics. If the global approach fails, this checkpoint is the known recovery point.
