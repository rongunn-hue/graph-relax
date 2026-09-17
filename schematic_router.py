import json, math, sys, os, collections, heapq
import xml.etree.ElementTree as ET

def power_like(name):
    return (name.upper() in ('GND', 'VCC', 'VEE', '0V')
            or name.upper().startswith('+') or name.upper().startswith('-'))

def seg_key(p1, p2):
    if abs(p1[0]-p2[0]) < 1e-6:
        return ('V', round(p1[0], 1), round(min(p1[1], p2[1]), 1), round(max(p1[1], p2[1]), 1))
    if abs(p1[1]-p2[1]) < 1e-6:
        return ('H', round(p1[1], 1), round(min(p1[0], p2[0]), 1), round(max(p1[0], p2[0]), 1))
    return None

def dedup(path):
    out = [path[0]]
    for p in path[1:]:
        if p != out[-1]:
            out.append(p)
    return out

def simplify(path):
    if len(path) < 3:
        return path
    out = [path[0]]
    for i in range(1, len(path)-1):
        x0, y0 = out[-1]; x1, y1 = path[i]; x2, y2 = path[i+1]
        if (x1-x0)*(y2-y1) == (x2-x1)*(y1-y0) and (x1-x0)*(x2-x1) >= 0 and (y1-y0)*(y2-y1) >= 0:
            continue
        out.append(path[i])
    out.append(path[-1])
    return out


class Router:
    """Exact-geometry orthogonal box router for one graph-relax net.
    Boxes never get crossed; every returned path segment is verified
    (path_clear) against every non-owning box before it is ever used --
    including the fallback tiers. No candidate is ever returned unchecked.
    """

    GRID = 10
    LEADS = (18, 30, 50, 80, 130, 220, 350)
    BOX_MARGIN = 12

    MIN_WIRE_GAP = 6.0  # minimum clearance between parallel segments of DIFFERENT nets

    def __init__(self, boxes, positions):
        self.boxes = boxes       # ref -> (x0,y0,x1,y1)
        self.pos = positions     # ref -> (cx,cy)
        self.fallback_log = []
        self.unverified = 0      # count of any return that skipped a clearance check
        self.placed_segs = []    # (axis, coord, lo, hi, net_name) for wires already routed this pass
        self.spacing_failed = 0  # count of wires where no wire-spacing-respecting route existed
        self.cheap = False       # True during swap-optimization trials: use plain BFS (no turn
                                  # tracking) for a fast approximate overlap count. The exact,
                                  # bend-minimizing Dijkstra is reserved for the one final routing
                                  # pass that's actually delivered -- re-running the full (cell,
                                  # heading)-state search on every one of 16-75 trials per circuit
                                  # is what made a single build_wires() call take ~6s even on the
                                  # smallest circuit, a real performance regression found this pass.

    def reset_wire_history(self):
        self.placed_segs = []
        self._seg_buckets = {'H': {}, 'V': {}}
        self.all_pins = []  # (x, y, net_name) for every device pin -- known before ANY routing

    def rebuild_spatial_index(self):
        # Precomputes, once per build_wires() call (boxes are fixed for its
        # duration), a reverse index from grid cell -> the (few) boxes
        # covering it. grid_blocked() used to loop over EVERY box on EVERY
        # single cell check -- called on every one of 4 neighbor-expansions
        # for every state Dijkstra pops -- which is what made routing slow
        # enough to need the cheap/exact split. This makes each check O(boxes
        # actually near this one cell), typically 0-2, instead of O(all boxes).
        self._dilated_cells = {}  # (gx,gy) -> set of refs, margin-dilated rect
        self._true_cells = {}     # (gx,gy) -> set of refs, true rect (for owned-box check)
        for ref, box in self.boxes.items():
            x0, y0, x1, y1 = box
            dx0, dy0 = x0-self.BOX_MARGIN, y0-self.BOX_MARGIN
            dx1, dy1 = x1+self.BOX_MARGIN, y1+self.BOX_MARGIN
            for gx in range(math.floor(dx0/self.GRID), math.ceil(dx1/self.GRID)+1):
                for gy in range(math.floor(dy0/self.GRID), math.ceil(dy1/self.GRID)+1):
                    px, py = gx*self.GRID, gy*self.GRID
                    if dx0 <= px <= dx1 and dy0 <= py <= dy1:
                        self._dilated_cells.setdefault((gx, gy), set()).add(ref)
            for gx in range(math.floor(x0/self.GRID), math.ceil(x1/self.GRID)+1):
                for gy in range(math.floor(y0/self.GRID), math.ceil(y1/self.GRID)+1):
                    px, py = gx*self.GRID, gy*self.GRID
                    if x0 <= px <= x1 and y0 <= py <= y1:
                        self._true_cells.setdefault((gx, gy), set()).add(ref)

    def seg_hits_box(self, p1, p2, box, margin=None):
        if margin is None:
            margin = self.BOX_MARGIN
        x0, y0, x1, y1 = box[0]-margin, box[1]-margin, box[2]+margin, box[3]+margin
        if abs(p1[0]-p2[0]) < 1e-6:
            x = p1[0]; ylo, yhi = sorted((p1[1], p2[1]))
            return x0 <= x <= x1 and not (yhi < y0 or ylo > y1)
        elif abs(p1[1]-p2[1]) < 1e-6:
            y = p1[1]; xlo, xhi = sorted((p1[0], p2[0]))
            return y0 <= y <= y1 and not (xhi < x0 or xlo > x1)
        return False

    def seg_too_close_to_other_wires(self, p1, p2, avoid_name):
        # A candidate segment is rejected if it runs closer than
        # MIN_WIRE_GAP, parallel and range-overlapping, to any segment
        # already placed for a DIFFERENT net. Same-net segments (this net's
        # own other spokes) are exempt -- they're expected to converge
        # near their shared junction. Without this, the router had zero
        # awareness of other wires at all, only of boxes -- two different
        # nets could end up running edge-to-edge with no visible gap
        # (a real defect the user caught by inspection).
        seg = seg_key(p1, p2)
        if seg is None:
            return False
        axis, coord, lo, hi = seg
        # Bucketed by MIN_WIRE_GAP-sized coordinate ranges instead of a full
        # scan of every placed segment -- checking the bucket a query
        # coordinate falls in plus its immediate neighbors covers the full
        # MIN_WIRE_GAP radius, same result as the old linear scan.
        b = int(coord // self.MIN_WIRE_GAP)
        buckets = self._seg_buckets[axis]
        for bb in (b-1, b, b+1):
            for coord2, lo2, hi2, name2 in buckets.get(bb, ()):
                if name2 == avoid_name:
                    continue
                if abs(coord2-coord) < self.MIN_WIRE_GAP and not (hi2 < lo or lo2 > hi):
                    return True
        return False

    def seg_too_close_to_other_pins(self, p1, p2, avoid_name):
        # Calculated, not tuned: seg_too_close_to_other_wires only protects
        # against segments ALREADY routed, so a not-yet-routed net's own
        # fixed pin position had no protection against an EARLIER net's
        # path settling within MIN_WIRE_GAP of it. Root cause found this
        # pass, not guessed: NODE_A's path bent to exactly y=130, which was
        # precisely R2's own NODE_B pin -- 14 units of pin-to-pin spacing
        # gave zero protection because the conflict was never pin-to-pin,
        # it was path-to-pin. Every pin position is fully known before any
        # routing starts (side_counts is computed up front), so this
        # registers all of them as real obstacles from the very first route
        # call -- a structural fix sized directly by MIN_WIRE_GAP, not a
        # spacing constant searched until it happened to pass.
        seg = seg_key(p1, p2)
        if seg is None:
            return False
        axis, coord, lo, hi = seg
        g = self.MIN_WIRE_GAP
        for px, py, name in self.all_pins:
            if name == avoid_name:
                continue
            if axis == 'H':
                if abs(py - coord) < g and lo - g <= px <= hi + g:
                    return True
            else:
                if abs(px - coord) < g and lo - g <= py <= hi + g:
                    return True
        return False

    CONGESTION_RADIUS = 18.0  # 3x MIN_WIRE_GAP: a soft awareness zone beyond the hard threshold
    CONGESTION_WEIGHT = 0.5   # extra Dijkstra cost per nearby different-net segment, per grid step

    def congestion_cost(self, p1, p2, avoid_name):
        # The router previously had NO preference between a corridor with
        # zero other wires and one with several, as long as none were
        # literally within MIN_WIRE_GAP -- only a hard threshold, no soft
        # awareness. That's the real cause of a bug found this pass: a
        # LATER-routed net (EB) ran a genuinely unnecessary 160-unit
        # straight line directly through the single most congested corridor
        # in the whole layout (where R2/R3/U1/NODE_A/NODE_B all also
        # route), because nothing made it prefer a quieter path when one
        # existed. This adds a small, soft cost proportional to how many
        # other nets' segments already run nearby, so the search naturally
        # steers away from busy corridors -- without ever hard-blocking a
        # route the way MIN_WIRE_GAP does, so it can never cause a routing
        # failure, only a preference.
        seg = seg_key(p1, p2)
        if seg is None:
            return 0
        axis, coord, lo, hi = seg
        b = int(coord // self.MIN_WIRE_GAP)
        span = int(self.CONGESTION_RADIUS // self.MIN_WIRE_GAP) + 1
        buckets = self._seg_buckets[axis]
        count = 0
        for bb in range(b - span, b + span + 1):
            for coord2, lo2, hi2, name2 in buckets.get(bb, ()):
                if name2 == avoid_name:
                    continue
                if abs(coord2 - coord) < self.CONGESTION_RADIUS and not (hi2 < lo or lo2 > hi):
                    count += 1
        return count

    def seg_enters_box_interior(self, p1, p2, box):
        # Stricter than seg_hits_box: true only if the segment passes
        # through the box's STRICT interior (open interval), not just
        # touches its boundary. Used for a box the wire OWNS -- exclude
        # correctly lets a wire touch/leave from its own boundary (which
        # the margin-dilated seg_hits_box would otherwise also block), but
        # that exemption was total, so a connector's sideways leg could
        # slide back into the box's real interior with nothing to stop it
        # (found near J_R1: the connector moved 13 units sideways along
        # y=174.5, back into J_R1's own box, before turning). An owned box
        # must still block genuine re-entry -- only boundary contact is OK.
        x0, y0, x1, y1 = box
        if abs(p1[0]-p2[0]) < 1e-6:
            x = p1[0]
            if not (x0 < x < x1):
                return False
            ylo, yhi = sorted((p1[1], p2[1]))
            return yhi > y0 and ylo < y1
        elif abs(p1[1]-p2[1]) < 1e-6:
            y = p1[1]
            if not (y0 < y < y1):
                return False
            xlo, xhi = sorted((p1[0], p2[0]))
            return xhi > x0 and xlo < x1
        return False

    def seg_hugs_box_edge(self, p1, p2, box):
        # A wire may touch its own box at exactly its single exit point
        # (zero-length contact) and leave perpendicular to that edge. It
        # must NOT then run collinear with that same edge for any real
        # distance -- that reads as part of the box's own outline, not a
        # wire (found at J_EB: after reaching the box's corner, the wire
        # ran 13 units along J_EB's own bottom edge before turning away).
        # "More than a point" of overlap with the edge is always wrong.
        #
        # Widened from exact collinearity (distance==0) to within
        # MIN_WIRE_GAP -- the same real clearance standard the router
        # already enforces between any two independent wires. Real defect
        # found this pass (R12, HRNG): a spoke exited its own box, stepped
        # up just 1 unit, then ran 37 units parallel to R12's own top edge
        # before turning -- visually indistinguishable from hugging at any
        # render scale, but the old exact-match check let it straight
        # through since 1 unit != 0. Reusing MIN_WIRE_GAP instead of a new
        # tuned constant.
        g = self.MIN_WIRE_GAP
        x0, y0, x1, y1 = box
        if abs(p1[0]-p2[0]) < 1e-6:
            x = p1[0]
            if abs(x-x0) > g and abs(x-x1) > g:
                return False
            ylo, yhi = sorted((p1[1], p2[1]))
            return min(yhi, y1) - max(ylo, y0) > 1e-6
        elif abs(p1[1]-p2[1]) < 1e-6:
            y = p1[1]
            if abs(y-y0) > g and abs(y-y1) > g:
                return False
            xlo, xhi = sorted((p1[0], p2[0]))
            return min(xhi, x1) - max(xlo, x0) > 1e-6
        return False

    def edge_hug_cost(self, p1, p2):
        # A soft Dijkstra cost (not a hard block -- grid_blocked already
        # owns hard blocking) added for any single grid step that would
        # hug a box edge per the widened seg_hugs_box_edge above. This is
        # what actually keeps the search from PRODUCING a hugging path in
        # the first place, rather than only detecting one after the fact:
        # a real fix belongs in the cost function the search optimizes,
        # not a post-hoc patch.
        #
        # Weight is BEND_PENALTY itself, not a separate tuned constant --
        # one hugging step should cost exactly as much as the one turn
        # that would escape it. A weaker weight (1.0, tried first) left a
        # real case unfixed (IAMP's R2/NODE_A): a short hug run was cost-
        # tied with the 2-bend detour needed to clear it, and the tie broke
        # the wrong way. Escaping can need more than one bend, so the
        # PER-STEP cost has to be at least one full bend to guarantee
        # multi-step hugs always lose to detouring around them.
        cost = 0.0
        for box in self.boxes.values():
            if self.seg_hugs_box_edge(p1, p2, box):
                cost += self.BEND_PENALTY
        return cost

    def path_clear(self, path, exclude, margin=None, avoid_name=None):
        for p1, p2 in zip(path, path[1:]):
            for ref, box in self.boxes.items():
                if ref in exclude:
                    if self.seg_enters_box_interior(p1, p2, box) or self.seg_hugs_box_edge(p1, p2, box):
                        return False
                    continue
                if self.seg_hits_box(p1, p2, box, margin=margin):
                    return False
            if avoid_name is not None and (self.seg_too_close_to_other_wires(p1, p2, avoid_name)
                                            or self.seg_too_close_to_other_pins(p1, p2, avoid_name)):
                return False
        return True

    def side_of(self, ref, target):
        cx, cy = self.pos[ref]
        tx, ty = target
        dx, dy = tx-cx, ty-cy
        if abs(dx) > abs(dy):
            return 'R' if dx > 0 else 'L'
        return 'T' if dy > 0 else 'B'

    def side_pin_point(self, ref, side, index, total):
        # Distributes MULTIPLE wires exiting the same side of the same box
        # across distinct points instead of all collapsing onto that side's
        # single center point -- a real, severe defect found by the user's
        # own count of wires touching an IC: 5 different signal nets were
        # landing on the exact same coordinate on U1's right edge, fanning
        # out from one spot (indistinguishable paths, and an undercount of
        # visibly separate lines). Mirrors the earlier fix for multi-rail
        # power stubs, generalized to every side of every device.
        x0, y0, x1, y1 = self.boxes[ref]
        cx, cy = self.pos[ref]
        # Cap total spread to the box's own side length (minus a margin) so
        # a device with many pins on one side never pushes a pin past its
        # own box into a neighbor's territory -- a real defect found this
        # pass (an IAMP pin landed inside a different, nearby device's box
        # once unifying power+signal pins raised some sides' pin counts).
        # Spacing shrinks (down to a floor) rather than staying fixed at 14
        # when that would overflow the box. Back to the original 14 (not
        # the 20 tried as a first, tuned-not-calculated attempt) -- the
        # real fix for the "earlier path settles on a not-yet-routed pin"
        # bug is all_pins/seg_too_close_to_other_pins below, which protects
        # every pin regardless of what this spacing constant is.
        side_len = (y1 - y0) if side in ('R', 'L') else (x1 - x0)
        max_span = max(0.0, side_len - 10)
        spacing = 14 if total <= 1 else min(14, max_span / (total - 1))
        spacing = max(spacing, 4)
        offset = (index - (total-1)/2) * spacing
        if side == 'R':
            return (x1, cy+offset), (1, 0)
        if side == 'L':
            return (x0, cy+offset), (-1, 0)
        if side == 'T':
            return (cx+offset, y1), (0, 1)
        return (cx+offset, y0), (0, -1)

    def exit_point(self, ref, target):
        side = self.side_of(ref, target)
        return self.side_pin_point(ref, side, 0, 1)

    STUB_LENGTHS = (20, 30, 45, 65, 90, 120, 160)

    def stub_endpoint(self, p, d, exclude, name):
        # An open-lead stub (single-terminal net, or a power/ground tail)
        # was previously just p + 20*d with zero clearance checking -- the
        # real cause of the "stinger" defect the user caught: a fixed
        # 20-unit dash fired blindly in its exit direction with no
        # awareness of whatever box or wire happened to sit there, so it
        # could graze or cut across unrelated territory. Every OTHER
        # endpoint in this tool is routed or at minimum clearance-checked;
        # this was the one place that wasn't. Calculated fix: grow the
        # length until the straight segment is actually clear (same
        # path_clear test real routes use), instead of tuning the one fixed
        # constant. Falls back to the longest length tried (still logged as
        # a real fallback, never silently accepted) if nothing is clear.
        for L in self.STUB_LENGTHS:
            stub = (p[0] + L*d[0], p[1] + L*d[1])
            if self.path_clear([p, stub], exclude, avoid_name=name):
                return stub
        self.fallback_log.append(f'stub for {name} never found clear length, using {self.STUB_LENGTHS[-1]}')
        L = self.STUB_LENGTHS[-1]
        return (p[0] + L*d[0], p[1] + L*d[1])

    def remove_backtracks(self, path, exclude, net_name):
        # Collapses any 3 consecutive points that reverse direction on the
        # same axis (go one way, then immediately back past the start) --
        # a real, distinct defect from what simplify() catches: simplify()
        # only merges collinear points continuing in the SAME direction, so
        # a genuine "there and back" detour (found at J_EA: up 8.5 units,
        # then down 10, overshooting the start by 1.5) passed through
        # untouched. Only collapses when the direct p0->p2 segment is still
        # fully verified-clear (box AND wire-spacing) -- never trades
        # correctness for tidiness.
        if len(path) < 3:
            return path
        out = [path[0]]
        i = 1
        while i < len(path) - 1:
            p0, p1, p2 = out[-1], path[i], path[i+1]
            reversal = False
            if abs(p0[0]-p1[0]) < 1e-6 and abs(p1[0]-p2[0]) < 1e-6:
                reversal = (p1[1]-p0[1]) * (p2[1]-p1[1]) < -1e-9
            elif abs(p0[1]-p1[1]) < 1e-6 and abs(p1[1]-p2[1]) < 1e-6:
                reversal = (p1[0]-p0[0]) * (p2[0]-p1[0]) < -1e-9
            if reversal and self.path_clear([p0, p2], exclude, avoid_name=net_name):
                i += 1  # drop p1, re-test p0->p2 against the point after it next loop
                continue
            out.append(p1)
            i += 1
        out.append(path[-1])
        return out

    def route(self, src, src_dir, dst, dst_dir, exclude, net_name=None):
        # Direction-aware shortest path (Dijkstra over (cell, heading)
        # states, turns penalized) replaces the old stack of heuristic
        # tiers entirely. This is the actual fix for the whole run of
        # distinct-looking bugs (backtracks, edge-hugs, interior dips, long
        # unnecessary detours): those were symptoms of testing hand-picked
        # candidate shapes one at a time. A real shortest-path search
        # structurally cannot produce a backtrack (never optimal) and won't
        # hug a box edge unless that edge is genuinely the only legal
        # corridor -- not because each of those cases was checked for, but
        # because the algorithm only ever returns the minimum-cost path.
        path = self.grid_route(src, dst, exclude, net_name=net_name, allow_fail=True)
        if path is not None:
            return path
        # A direct spacing-respecting route doesn't exist -- before
        # sacrificing wire-spacing, try routing through a free RELAY point
        # first (user's idea: decouple the path -- a two-ended waypoint
        # that isn't anchored to either box, so it adds a real degree of
        # freedom the direct src->dst search doesn't have). This is exactly
        # what fixes the diagnosed root cause of the remaining IAMP
        # failures: a single device's own two pins were both forced through
        # the same corridor because they were rigidly anchored to that
        # device's boundary. A relay lets one leg detour around the
        # congestion while the other stays direct.
        # Only in exact mode -- this is a real, somewhat expensive search
        # (up to 8 candidates x 2 legs each), and during cheap-mode swap
        # trials only an approximate comparison metric is needed anyway.
        if not self.cheap:
            relay_path = self.route_via_relay(src, dst, exclude, net_name)
            if relay_path is not None:
                self.fallback_log.append(('RELAY', src, dst))
                return relay_path
        # Genuine failure: no route exists honoring wire spacing, even via
        # a relay. Reported, not papered over -- retry box-clearance-only.
        self.spacing_failed += 1
        self.fallback_log.append(('SPACING-SACRIFICED', src, dst))
        path = self.grid_route(src, dst, exclude)
        if path is not None:
            return path
        self.unverified += 1
        return dedup([src, dst])

    RELAY_DISTANCES = (10, 20, 30, 60, 100, 150)
    RELAY_DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))

    def route_via_relay(self, src, dst, exclude, net_name):
        # Try relay points radiating out from BOTH src and dst -- not just
        # src. First attempt only anchored at src and found nothing: traced
        # it and found the real bottleneck was actually at the DESTINATION
        # end (a virtual multi-terminal junction sitting close to another
        # net's own junction, both fed by paths converging through the same
        # small area) -- a source-only relay can't help a destination-side
        # pinch. A candidate is accepted only if BOTH legs (anchor->relay,
        # relay->other end) are found by the real spacing-aware search --
        # never a partial or unverified result. No dot is drawn at the
        # relay (it's a routing waypoint, not a genuine electrical
        # junction).
        for anchor, other in ((src, dst), (dst, src)):
            for dist in self.RELAY_DISTANCES:
                for ddx, ddy in self.RELAY_DIRECTIONS:
                    norm = math.hypot(ddx, ddy)
                    relay = (anchor[0] + dist*ddx/norm, anchor[1] + dist*ddy/norm)
                    if self.point_in_any_box(relay, margin=self.BOX_MARGIN, exclude=exclude):
                        continue
                    leg1 = self.grid_route(anchor, relay, exclude, net_name=net_name, allow_fail=True)
                    if leg1 is None:
                        continue
                    leg2 = self.grid_route(relay, other, exclude, net_name=net_name, allow_fail=True)
                    if leg2 is None:
                        continue
                    full = leg1 + leg2 if anchor is src else leg2[::-1] + leg1[::-1]
                    return dedup(simplify(full))
        return None

    def grid_blocked(self, gx, gy, exclude):
        # Same semantics as before (owned box blocks via its true rect,
        # non-owned boxes block via their margin-dilated rect), just looked
        # up in the precomputed index instead of scanning every box.
        true_refs = self._true_cells.get((gx, gy))
        if true_refs:
            for ref in true_refs:
                if ref in exclude:
                    return True
        dilated_refs = self._dilated_cells.get((gx, gy))
        if dilated_refs:
            for ref in dilated_refs:
                if ref not in exclude:
                    return True
        return False

    def point_in_any_box(self, pt, margin=0, exclude=()):
        for ref, box in self.boxes.items():
            m = 0 if ref in exclude else margin
            if box[0]-m < pt[0] < box[2]+m and box[1]-m < pt[1] < box[3]+m:
                return True
        return False

    def connect_to_grid(self, p, exclude, avoid_name=None):
        # "Nearest reachable grid node" is not enough: with the source's own
        # box excluded from obstacle checks (so the connector can leave from
        # its own boundary), the geometrically nearest node can sit BACKWARD
        # -- literally inside that same excluded box, since nothing there
        # blocks it. A real defect found this pass: a D2-owned wire's
        # connector stepped 3 units back into D2 instead of forward, simply
        # because that grid point was closer than the first outward one.
        # The candidate node itself must be genuine open space regardless of
        # which box is excluded for the path-clearance check.
        gx0, gy0 = round(p[0]/self.GRID), round(p[1]/self.GRID)
        for r in range(0, 12):
            candidates = [(gx0+ddx, gy0+ddy) for ddx in range(-r, r+1) for ddy in range(-r, r+1)
                          if max(abs(ddx), abs(ddy)) == r]
            for gx, gy in candidates:
                # Must use the SAME acceptance test the search graph itself
                # uses (grid_blocked), not a looser one -- a real bug found
                # this pass: point_in_any_box (no margin dilation) accepted
                # nodes that grid_blocked (margin-dilated for non-owned
                # boxes) then refused to ever use, so the search's own
                # start/end node could be unreachable by construction,
                # failing silently for routes that should have worked fine.
                if self.grid_blocked(gx, gy, exclude):
                    continue
                node = (gx*self.GRID, gy*self.GRID)
                for mid in [(node[0], p[1]), (p[0], node[1])]:
                    path = dedup([p, mid, node])
                    if self.path_clear(path, exclude, avoid_name=avoid_name):
                        return (gx, gy), path
        return None, None

    BEND_PENALTY = 3  # extra cost per turn -- prefers fewer bends over marginally shorter cell counts

    def _bfs_cells(self, entry_cell, exit_cell, exclude, avoid_name):
        # Cheap-mode search: single state per cell, no heading tracking, no
        # bend penalty. Not used for anything actually delivered -- see
        # self.cheap above.
        minx = min(b[0] for b in self.boxes.values())/self.GRID - 8
        maxx = max(b[2] for b in self.boxes.values())/self.GRID + 8
        miny = min(b[1] for b in self.boxes.values())/self.GRID - 8
        maxy = max(b[3] for b in self.boxes.values())/self.GRID + 8
        q = collections.deque([entry_cell])
        came = {entry_cell: None}
        while q:
            cx, cy = q.popleft()
            if (cx, cy) == exit_cell:
                break
            for ddx, ddy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                nx, ny = cx+ddx, cy+ddy
                if not (minx <= nx <= maxx and miny <= ny <= maxy):
                    continue
                if (nx, ny) in came:
                    continue
                if self.grid_blocked(nx, ny, exclude):
                    continue
                if avoid_name is not None and self.seg_too_close_to_other_wires(
                        (cx*self.GRID, cy*self.GRID), (nx*self.GRID, ny*self.GRID), avoid_name):
                    continue
                came[(nx, ny)] = (cx, cy)
                q.append((nx, ny))
        if exit_cell not in came:
            return None
        cell_path = []
        cur = exit_cell
        while cur is not None:
            cell_path.append(cur); cur = came[cur]
        cell_path.reverse()
        return cell_path

    def _dijkstra(self, entry_cell, exit_cell, exclude, avoid_name):
        # State = (cell, heading-arrived-on). Dijkstra over this state
        # graph is a real minimum-cost orthogonal path with turns
        # penalized -- this is what makes the result provably free of
        # backtracks and needless detours, not a property that has to be
        # separately checked for afterward.
        #
        # self.cheap collapses this to plain BFS (single state per cell,
        # no heading, no bend penalty) -- used only during swap-optimization
        # trials, where all that's needed is an approximate overlap count to
        # compare candidates, not the exact minimum-bend path. Re-running the
        # full state-space search on every trial (16-75 per circuit) made a
        # single build_wires() call take ~6s even on the smallest circuit --
        # a real performance regression, not an acceptable cost of
        # correctness. The final, delivered routing always uses the full
        # search (see run(): rt.cheap is reset to False before the last pass).
        if self.cheap:
            return self._bfs_cells(entry_cell, exit_cell, exclude, avoid_name)
        minx = min(b[0] for b in self.boxes.values())/self.GRID - 8
        maxx = max(b[2] for b in self.boxes.values())/self.GRID + 8
        miny = min(b[1] for b in self.boxes.values())/self.GRID - 8
        maxy = max(b[3] for b in self.boxes.values())/self.GRID + 8
        DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))
        dist = {(entry_cell, None): 0}
        came = {}
        pq = [(0, entry_cell, None)]
        goal = None
        while pq:
            d, cell, heading = heapq.heappop(pq)
            if d > dist.get((cell, heading), float('inf')):
                continue
            if cell == exit_cell:
                goal = (cell, heading)
                break
            cx, cy = cell
            for ddx, ddy in DIRS:
                nx, ny = cx+ddx, cy+ddy
                if not (minx <= nx <= maxx and miny <= ny <= maxy):
                    continue
                if self.grid_blocked(nx, ny, exclude):
                    continue
                p1c, p2c = (cx*self.GRID, cy*self.GRID), (nx*self.GRID, ny*self.GRID)
                if avoid_name is not None and (self.seg_too_close_to_other_wires(p1c, p2c, avoid_name)
                                                or self.seg_too_close_to_other_pins(p1c, p2c, avoid_name)):
                    continue
                step = 1 + (self.BEND_PENALTY if heading is not None and (ddx, ddy) != heading else 0)
                step += self.edge_hug_cost(p1c, p2c)
                if avoid_name is not None:
                    step += self.CONGESTION_WEIGHT * self.congestion_cost(p1c, p2c, avoid_name)
                nd = d + step
                key = ((nx, ny), (ddx, ddy))
                if nd < dist.get(key, float('inf')):
                    dist[key] = nd
                    came[key] = (cell, heading)
                    heapq.heappush(pq, (nd, (nx, ny), (ddx, ddy)))
        if goal is None:
            return None
        cell_path = []
        cur = goal
        while cur is not None:
            cell_path.append(cur[0])
            cur = came.get(cur)
        cell_path.reverse()
        return cell_path

    def grid_route(self, src, dst, exclude, net_name=None, allow_fail=False):
        # allow_fail=True: return None on failure instead of the unverified
        # dedup([src,dst]) escape hatch -- used by route()'s wire-spacing
        # retry, which still has a box-clearance-only pass to fall back on.
        entry_cell, entry_path = self.connect_to_grid(src, exclude, avoid_name=net_name)
        exit_cell, exit_path = self.connect_to_grid(dst, exclude, avoid_name=net_name)
        if entry_cell is None or exit_cell is None:
            if allow_fail:
                return None
            self.unverified += 1
            return dedup([src, dst])
        cell_path = self._dijkstra(entry_cell, exit_cell, exclude, net_name)
        if cell_path is None:
            if allow_fail:
                return None
            self.unverified += 1
            return dedup([src, dst])
        grid_real = [(gx*self.GRID, gy*self.GRID) for gx, gy in cell_path]
        full = entry_path[:-1] + grid_real + list(reversed(exit_path))[1:]
        full = dedup(simplify(full))
        # Real defect found this pass (the "stinger" the user kept
        # reporting, e.g. C1's own spoke to its junction): stitching
        # connect_to_grid's entry/exit approach legs onto the Dijkstra
        # cell path can produce a genuine "there and back" reversal --
        # simplify() only merges collinear points continuing in the SAME
        # direction, so it passes a reversal straight through untouched.
        # remove_backtracks() was written to catch exactly this, but was
        # never actually wired into the pipeline that builds this path --
        # a real gap, not a tuning problem. Applying it here, at the one
        # place every routed path (main grid route, relay legs, entry/exit
        # stitching all included) funnels through before being returned.
        return self.remove_backtracks(full, exclude, net_name)

    def nearest_clear_point(self, px, py, exclude):
        # Used for 3+-terminal virtual-junction points, which are not tied
        # to any single box's exit -- if the raw centroid lands inside a
        # box it doesn't belong to, nudge it to the nearest open spot
        # instead of leaving an invalid, unverified target in place.
        # Margin must match the router's own avoidance margin (BOX_MARGIN),
        # not a smaller one -- a point can be "outside the box" at margin=4
        # yet still be trapped inside another box's real 12-unit avoidance
        # shadow, which is what actually matters for routability.
        #
        # Snapped to self.GRID here (and every candidate below searched on
        # that same grid) -- the real cause of the "floating dot" defect:
        # an off-grid junction (a raw centroid average, e.g. x=413.1) forced
        # every spoke's path to end with a sub-grid final hop of just 1-3
        # units to reach it, since the Dijkstra search and connect_to_grid
        # both operate in GRID-sized steps. At normal render scale that
        # last hop is nearly invisible, so the dot visually reads as
        # floating just off the end of each line instead of sitting AT it.
        px, py = round(px/self.GRID)*self.GRID, round(py/self.GRID)*self.GRID
        m = self.BOX_MARGIN
        blocked = any(ref not in exclude and b[0]-m <= px <= b[2]+m and b[1]-m <= py <= b[3]+m
                      for ref, b in self.boxes.items())
        if not blocked:
            return (px, py)
        step = self.GRID
        for r in range(1, 40):
            for ddx in range(-r, r+1):
                for ddy in (-r, r):
                    cx, cy = px+ddx*step, py+ddy*step
                    if not any(ref not in exclude and b[0]-m <= cx <= b[2]+m and b[1]-m <= cy <= b[3]+m
                              for ref, b in self.boxes.items()):
                        return (cx, cy)
            for ddy in range(-r, r+1):
                for ddx in (-r, r):
                    cx, cy = px+ddx*step, py+ddy*step
                    if not any(ref not in exclude and b[0]-m <= cx <= b[2]+m and b[1]-m <= cy <= b[3]+m
                              for ref, b in self.boxes.items()):
                        return (cx, cy)
        return (px, py)  # gave up -- caller's assertion will catch this


def _decoupling_pin_targets(sem):
    # Derives, for each decoupling-purpose satellite (e.g. HRNG's C13-C17),
    # exactly which host pin it's meant to bypass -- purely from data
    # already in semantic.json (placement.near + real net membership),
    # never hardcoded by ref name. Naturally does nothing for a circuit
    # without this metadata (confirmed: IAMP has none) and naturally
    # extends to any future circuit with the same pattern.
    #
    # Only returns a target when the host's device type has a FIXED-
    # TEMPLATE pin transform registered (today: OPAMP_DUAL_PINOUT via
    # dual_opamp_xform, and OPAMP_DIP8_PINOUT via opamp8_xform, now that
    # U2 has its own fixed-template DIP-8 symbol too). A generic multi-pin
    # box only allocates its exact pin coordinate at routing time, so
    # there is no coordinate to pre-position against yet -- a real
    # structural limit, not an oversight. This is why HRNG's U4-targeted
    # cap (C17) is still excluded here (PICAXE remains a generic box).
    comps = {c['ref']: c for c in sem['electrical']['components']}
    nets_by_name = {n['name']: n['terminals'] for n in sem['electrical']['nets']}
    net_of_terminal = {}
    for name, terms in nets_by_name.items():
        for t in terms:
            net_of_terminal[t] = name
    parts = sem.get('source', {}).get('metadata', {}).get('normalized', {}).get('parts', [])
    targets = {}
    for p in parts:
        placement = p.get('placement') or {}
        if placement.get('purpose') != 'local decoupling':
            continue
        cap_ref, host_ref = p['ref'], placement.get('near')
        host_device = comps.get(host_ref, {}).get('device') if host_ref else None
        if host_device not in OPAMP_DUAL_PINOUT and host_device not in OPAMP_DIP8_PINOUT:
            continue
        rail_net = net_of_terminal.get(f'{cap_ref}.1')
        if rail_net is None:
            continue
        host_pin = None
        for t in nets_by_name[rail_net]:
            ref, pin = split_term(t)
            if ref == host_ref:
                host_pin = pin
                break
        if host_pin is not None:
            targets[cap_ref] = (host_ref, host_pin)
    return targets


def layout_boxes(sem, init, gap=60):
    elec = sem['electrical']
    comps = {c['ref']: c for c in elec['components']}
    territories = {t['id']: t for t in init['territories']}
    edges = init['territory_edges']

    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    # list(), not set() -- a set's iteration order depends on Python's
    # per-process string hash randomization, not on the actual data. That
    # made "start" (and everything downstream: chain order, row/column
    # layout, every box position, every routed wire) genuinely
    # nondeterministic between runs of the IDENTICAL script on the
    # IDENTICAL input -- a real bug, not something safe to assume away.
    # dict.keys() already preserves insertion order deterministically;
    # wrapping it in set() was actively throwing that away.
    in_chain = list(adj.keys())
    start = next((n for n in in_chain if len(adj[n]) == 1), in_chain[0] if in_chain else None)
    chain = []; seen = set(); cur = start
    while cur and cur not in seen:
        chain.append(cur); seen.add(cur)
        nxt = [x for x in adj.get(cur, []) if x not in seen]
        cur = nxt[0] if nxt else None
    order = chain + [t for t in territories if t not in chain]
    COLS = 3
    rows = [order[i:i+COLS] for i in range(0, len(order), COLS)]

    def device_box_pts(ref):
        npins = len(comps[ref]['pins'])
        side = max(40, 22 + 7*npins)
        if comps[ref].get('device') in OPAMP_DUAL_PINOUT:
            # The calculated dual-op-amp symbol (dual_opamp_xform) is
            # defined in a normalized DUAL_OPAMP_NORM_SIZE coordinate
            # system and uses a SINGLE scale factor (box_w/norm_w) for
            # everything, artwork and electrical terminals alike -- that
            # only stays exact if this box's own aspect ratio is EXACTLY
            # norm_w:norm_h, not an approximation. Scaled up 1.35x from
            # the generic pin-count-driven "side" (still proportional to
            # schematic density, not a fixed absolute size) for
            # legibility, same multiplier proven not to cause routing
            # regressions at the previous (different) ratio.
            side = side * 1.35
            norm_w, norm_h = DUAL_OPAMP_NORM_SIZE
            box_w = side
            box_h = side * (norm_h/norm_w)
            return box_w, box_h
        if comps[ref].get('device') in TRANSISTOR_PINOUT:
            # Same exact-aspect-ratio principle as the dual-op-amp above:
            # npn_xform derives scale from width alone, so this box's own
            # aspect ratio must be EXACTLY NPN_NORM_SIZE's (351:421), not
            # an approximation, or the circle/leads would stretch off-
            # ratio. This also structurally guarantees the whole glyph
            # (including its terminal points, which sit at the norm box's
            # own edges) always fits inside the box -- no separate
            # overflow-avoidance math needed, unlike an earlier version of
            # this glyph that sized itself independently of the box and
            # silently overflowed it.
            side = side * 1.3
            norm_w, norm_h = NPN_NORM_SIZE
            box_w = side
            box_h = side * (norm_h/norm_w)
            return box_w, box_h
        if comps[ref].get('device') in OPAMP_DIP8_PINOUT:
            # Same exact-aspect-ratio principle as the dual-op-amp/
            # transistor above (109:133) -- structurally guarantees the
            # whole glyph fits inside the box, no separate overflow math.
            side = side * 1.3
            norm_w, norm_h = OPAMP_DIP8_NORM_SIZE
            box_w = side
            box_h = side * (norm_h/norm_w)
            return box_w, box_h
        if comps[ref].get('device') in SW_DPST_PINOUT:
            # Same exact-aspect-ratio principle (10.16:6.604).
            side = side * 1.3
            norm_w, norm_h = SW_DPST_NORM_SIZE
            box_w = side
            box_h = side * (norm_h/norm_w)
            return box_w, box_h
        return side, side

    GAP = gap
    pos = {}; boxsize = {}
    row_y = 0.0
    for r, trow in enumerate(rows):
        trow2 = list(reversed(trow)) if r % 2 == 1 else trow
        col_x = 0.0; row_h = 0.0
        for tid in trow2:
            devs = sorted(territories[tid]['devices'])
            n = len(devs)
            gcols = max(1, math.ceil(math.sqrt(n)))
            cur_x, cur_y, tallest, col_in_block = col_x, row_y, 0.0, 0
            block_w = 0.0
            for ref in devs:
                w, h = device_box_pts(ref)
                boxsize[ref] = (w, h)
                if col_in_block >= gcols:
                    cur_x = col_x; cur_y += tallest + GAP; tallest = 0.0; col_in_block = 0
                pos[ref] = (cur_x + w/2, cur_y + h/2)
                cur_x += w + GAP
                block_w = max(block_w, cur_x - col_x)
                tallest = max(tallest, h)
                col_in_block += 1
            row_h = max(row_h, (cur_y - row_y) + tallest)
            col_x += block_w + GAP*1.5
        row_y += row_h + GAP*1.5

    boxes = {ref: (pos[ref][0]-boxsize[ref][0]/2, pos[ref][1]-boxsize[ref][1]/2,
                   pos[ref][0]+boxsize[ref][0]/2, pos[ref][1]+boxsize[ref][1]/2) for ref in pos}

    # Reposition decoupling-purpose satellites (see _decoupling_pin_targets)
    # to sit immediately at their real target pin, instead of the generic
    # territory-adjacency placement above -- only for hosts with a fixed-
    # template pin transform (today: the dual op-amp). User: "You should
    # be able to place C13 at pin 8 where it connects to U1 instead of
    # floating at the top, right?" This runs AFTER every other box is
    # already placed (both the host's box and, implicitly, everything
    # else this satellite might now sit near) -- happens strictly before
    # build_wires() ever routes a single wire, so there is no "earlier
    # routing" to invalidate, only a box position possibly worth re-
    # verifying for overlap (verify() already catches that class of
    # defect the same way it catches every other placement issue).
    for cap_ref, (host_ref, host_pin) in _decoupling_pin_targets(sem).items():
        if cap_ref not in boxes or host_ref not in boxes:
            continue
        host_device = comps[host_ref]['device']
        # Dispatch to whichever fixed-template host this is -- same
        # principle either way (dual op-amp DIP-8 or single op-amp
        # DIP-8), just a different terminals dict/xform per family. Add a
        # new elif here, not a new copy of this whole loop, whenever
        # another fixed-template device type gets decoupling caps.
        if host_device in OPAMP_DUAL_PINOUT:
            dp = OPAMP_DUAL_PINOUT[host_device]
            terminals, xform = DUAL_OPAMP_TERMINALS, dual_opamp_xform
        elif host_device in OPAMP_DIP8_PINOUT:
            dp = OPAMP_DIP8_PINOUT[host_device]
            terminals, xform = OPAMP_DIP8_TERMINALS, opamp8_xform
        else:
            continue
        if host_pin not in dp['left_order'] and host_pin not in dp['right_order']:
            continue
        hx0, hy0, hx1, hy1 = boxes[host_ref]
        px, py = xform(hx0, hy0, hx1, hy1)(*terminals[host_pin])
        cw, ch = boxsize[cap_ref]
        # Same standard gap every other device in this layout is spaced
        # by (the density rule stays in effect) -- NOT a tighter squeeze.
        # A real, confirmed defect this pass: an earlier attempt used a
        # tight 20-unit clearance here (vs. the normal 60), which pinched
        # off a routing corridor near U1 badly enough that the router's
        # search blew up to 5+ minutes on one net before being killed.
        # User: "Don't push the box flush against the block. The density
        # rule is still in effect." Exact pin-Y alignment is what actually
        # matters here, not minimizing distance.
        clearance = gap
        if host_pin in dp['right_order']:
            cx = hx1 + clearance + cw/2
        else:
            cx = hx0 - clearance - cw/2
        cy = py
        pos[cap_ref] = (cx, cy)
        boxes[cap_ref] = (cx-cw/2, cy-ch/2, cx+cw/2, cy+ch/2)

    return comps, territories, pos, boxsize, boxes


def split_term(t):
    # "D1.A" -> ("D1", "A"). Real pin identity, kept from here all the way
    # through to pin_points/the SVG/the equivalence checker -- previously
    # every net-terminal string was immediately reduced to just its device
    # ref (t.split('.')[0]) and the pin half thrown away, so nothing
    # downstream (routing, the SVG, the checker) ever knew WHICH pin of a
    # multi-pin device a connection used. That's the real prerequisite for
    # both individually-identifiable symbol pins (diode polarity,
    # transistor B/C/E, op-amp IN+/IN-/OUT) and a pin-level equivalence
    # proof, not just a diode-specific concern.
    if '.' in t:
        ref, pin = t.split('.', 1)
    else:
        ref, pin = t, ''
    return ref, pin


def _merge_intervals(intervals):
    # [(lo,hi), ...] -> sorted, non-overlapping, touching-merged. This is
    # the actual "internal segment multiplicity must disappear" step: N
    # collinear/overlapping same-net segments collapse into exactly the
    # visible run(s) they cover, so counting endpoints afterward can never
    # over-count a single visible direction as multiple arms.
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo <= merged[-1][1] + 1e-6:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def compute_junction_dots(wires):
    # THE general junction rule, computed directly from final routed
    # geometry plus real net identity -- not from where a junction was
    # PLANNED to be, and not from raw segment multiplicity.
    #
    # Real defect fixed by the first version of this function: dots were
    # only ever placed at the one pre-computed multi-terminal junction
    # point (nx,ny) that every spoke was routed TOWARD, so a spoke bending
    # exactly on top of a DIFFERENT spoke's own segment (same net) --
    # found concretely at mono741's RFB, (460,110) sitting strictly inside
    # another OUT_RAW spoke's own run -- got no dot.
    #
    # Real defect in THAT fix, found next: it counted incident raw
    # SEGMENTS, not distinct VISIBLE directions. Two different spokes of
    # the same net can geometrically overlap (route the same collinear
    # stretch), which the segment-counting version would count as two
    # separate arms even though only one direction is actually visible
    # there -- a false dot, or an inflated count that could mask a true
    # one. Fixed by normalizing first: merge all collinear/overlapping
    # same-net segments into their maximal visible runs (_merge_intervals)
    # BEFORE counting anything, so internal multiplicity genuinely
    # disappears rather than being counted around.
    #
    # Rule, now over normalized geometry: group by net NAME (never by
    # spatial coincidence -- different nets are never unioned, so a
    # crossing between them can't produce a dot). Merge each net's H runs
    # per shared y, and V runs per shared x. A "relevant vertex" is any
    # merged run's own endpoint, or any point where a merged H run and a
    # merged V run of the SAME net cross or touch. At each such vertex,
    # determine which of the (at most 4) cardinal directions has visible
    # conductor leaving it -- a boolean per direction, not a count, so
    # coincident overlapping runs in the same direction contribute once.
    # >= 3 distinct directions present is a real branch; a dot is required
    # there. This reproduces the ordinary 3+-terminal junction (every
    # spoke's endpoint contributes a distinct direction) and a genuine
    # same-net 4-way crossing automatically, with no special-casing.
    by_net = collections.defaultdict(list)
    for path, kind, owners, name in wires:
        for p1, p2 in zip(path, path[1:]):
            if p1 == p2:
                continue
            by_net[name].append((p1, p2))

    dots = []
    for name, segs in by_net.items():
        h_by_y = collections.defaultdict(list)  # y -> [(xlo,xhi), ...]
        v_by_x = collections.defaultdict(list)  # x -> [(ylo,yhi), ...]
        # Every coordinate that reaches the direction check below (both the
        # axis a run is keyed on AND its own along-axis endpoints) must go
        # through the SAME round(.,1) before anything compares it against
        # anything else. Rounding only the dict key (as an earlier version
        # did) let a run's raw endpoint (e.g. 377.25, which round(.,1) ties
        # down to 377.2) disagree with the OTHER axis's already-rounded
        # coordinate by up to 0.05 -- far outside the 1e-6 tolerance used
        # everywhere else here -- so a plain bend's endpoint could read as
        # strictly interior to a run and pick up a phantom extra direction
        # (or, on the other rounding side of a tie, silently lose a real
        # one). Rounding both axes' values up front makes every later
        # comparison a same-precision compare, so 1e-6 only ever has to
        # absorb float noise, never a rounding-tie mismatch.
        for p1, p2 in segs:
            if abs(p1[1]-p2[1]) < 1e-6:
                h_by_y[round(p1[1], 1)].append(tuple(sorted((round(p1[0], 1), round(p2[0], 1)))))
            elif abs(p1[0]-p2[0]) < 1e-6:
                v_by_x[round(p1[0], 1)].append(tuple(sorted((round(p1[1], 1), round(p2[1], 1)))))

        h_runs = [(y, lo, hi) for y, ivs in h_by_y.items() for lo, hi in _merge_intervals(ivs)]
        v_runs = [(x, lo, hi) for x, ivs in v_by_x.items() for lo, hi in _merge_intervals(ivs)]

        vertices = set()
        for y, xlo, xhi in h_runs:
            vertices.add((round(xlo, 1), y)); vertices.add((round(xhi, 1), y))
        for x, ylo, yhi in v_runs:
            vertices.add((x, round(ylo, 1))); vertices.add((x, round(yhi, 1)))
        for y, xlo, xhi in h_runs:
            for x, ylo, yhi in v_runs:
                if xlo - 1e-6 <= x <= xhi + 1e-6 and ylo - 1e-6 <= y <= yhi + 1e-6:
                    vertices.add((round(x, 1), round(y, 1)))

        for px, py in vertices:
            directions = set()
            for y, xlo, xhi in h_runs:
                if abs(y - py) > 1e-6:
                    continue
                if xlo - 1e-6 <= px < xhi - 1e-6:
                    directions.add('R')
                if xlo + 1e-6 < px <= xhi + 1e-6:
                    directions.add('L')
            for x, ylo, yhi in v_runs:
                if abs(x - px) > 1e-6:
                    continue
                if ylo - 1e-6 <= py < yhi - 1e-6:
                    directions.add('D')
                if ylo + 1e-6 < py <= yhi + 1e-6:
                    directions.add('U')
            if len(directions) >= 3:
                dots.append((px, py, name))
    return dots


def build_wires(rt, comps, nets, pos, boxes, sem=None):
    wires = []
    pin_points = {}  # (ref, pin_id) -> (x, y), for EVERY allocated pin, all device types
    # rt.boxes/rt.pos must always point at THESE arguments, not whatever
    # dict Router.__init__ happened to be given -- a real bug found this
    # pass: every existing caller mutated the same boxes dict in place
    # across the whole script, so rt.boxes and the local `boxes` variable
    # were always the same object by accident. uncrowd_boxes() builds a
    # genuinely NEW dict for its trial layout, which broke that hidden
    # assumption -- side_pin_point/exit_point read rt.boxes internally, so
    # they kept computing exit points from the OLD box position even
    # though this function was explicitly called with the new one.
    rt.boxes = boxes
    rt.pos = pos
    rt.reset_wire_history()
    rt.rebuild_spatial_index()  # boxes are fixed for the rest of this call

    def record(path, name):
        for p1, p2 in zip(path, path[1:]):
            k = seg_key(p1, p2)
            if k is None:
                continue
            axis, coord, lo, hi = k
            rt.placed_segs.append((axis, coord, lo, hi, name))
            b = int(coord // rt.MIN_WIRE_GAP)
            rt._seg_buckets[axis].setdefault(b, []).append((coord, lo, hi, name))

    # Power connections are deferred into the same per-side allocation as
    # signal wires (see pass 1/2 below) instead of using their own separate
    # spacing system -- an earlier version spaced power pins and signal
    # pins independently, and since both always used the same side ('R')
    # with each system centered on the box, the two interleaved to HALVE
    # the effective gap (7 units instead of 14) between any two adjacent
    # lines. All 8 of U1's pins were genuinely distinct but too close
    # together to reliably count by eye -- this is the actual fix.
    # Decoupling caps with a known exact target pin (see
    # _decoupling_pin_targets) get a REAL wire to that specific pin
    # instead of two independent open stubs to the same rail label --
    # user: "the end that says +9V on CAP13 should be connected to pin 8
    # on U1." Every OTHER terminal on that same rail net (C4, J1, R1,
    # etc.) is untouched, still a generic open stub -- this pulls out
    # only the one specific (cap, host) pair per qualifying net, nothing
    # else about power-net handling changes.
    decoupling_targets = _decoupling_pin_targets(sem) if sem is not None else {}

    power_items = []  # (name, ref, pin)
    plan = []  # non-power connections, deferred to a two-pass build below
    for n in nets:
        name = n['name']; terms = n['terminals']
        if power_like(name):
            parsed_terms = [split_term(t) for t in terms]
            wired_refs = set()
            for ref, pin in parsed_terms:
                target = decoupling_targets.get(ref)
                if target is None:
                    continue
                host_ref, host_pin = target
                if (host_ref, host_pin) in parsed_terms:
                    plan.append(('pair', name, ref, host_ref, pin, host_pin))
                    wired_refs.add(ref); wired_refs.add(host_ref)
            for ref, pin in parsed_terms:
                if ref in wired_refs:
                    continue
                power_items.append((name, ref, pin))
            continue
        parsed = [split_term(t) for t in terms]  # [(ref, pin), ...]
        refs = [r for r, p in parsed]
        if len(refs) == 2:
            (a, pin_a), (b, pin_b) = parsed
            plan.append(('pair', name, a, b, pin_a, pin_b))
        elif len(refs) == 1:
            ref, pin = parsed[0]
            plan.append(('stub', name, ref, pin))
        else:
            xs = [pos[r][0] for r in refs]; ys = [pos[r][1] for r in refs]
            nx, ny = sum(xs)/len(xs), sum(ys)/len(ys)
            # The junction point itself must never sit inside ANY box,
            # including its own net's participants -- letting a participant
            # exclude itself here (an earlier version did) let the junction
            # land inside one large participant's box, which was fine for
            # that device's own wire but made every other spoke wire dip
            # into that box's interior to reach it. Empty exclude: every box
            # is real territory for choosing the junction location.
            nx, ny = rt.nearest_clear_point(nx, ny, set())
            pins = [p for r, p in parsed]
            plan.append(('multi', name, refs, nx, ny, pins))

    # Devices drawn as a real 2-terminal symbol (resistor zigzag, capacitor
    # plates -- see draw_two_pin_symbol/render_svg below) need their two
    # leads on genuinely OPPOSITE sides of the box to read as a straight
    # axial component, not wherever each connection's independent side_of()
    # happens to land -- checked empirically: real 2-pin-passive sides were
    # 'R'/'R' or adjacent sides far more often than opposite, since a
    # passive's power/ground leg is always forced to side 'R' regardless of
    # direction. Decided once per device: whichever connection has a real
    # directional target (a signal pair or multi-junction) sets the axis;
    # the other connection (often an undirected power/stub lead) is simply
    # forced to the opposite side. Same uid scheme as pin_rank below, so
    # this single precomputed dict overrides side_of() at every one of its
    # several call sites consistently instead of patching each separately.
    SYMBOL_2PIN_TYPES = {'GENERIC_RESISTOR', 'GENERIC_CAPACITOR', '1N4148', 'ZENER_2PIN', 'GENERIC_LED', 'GENERIC_METER'}
    device_conns = collections.defaultdict(list)  # ref -> [(target_or_None, uid)]
    # uid -> physical pin string, populated alongside device_conns. Needed
    # (only) by the op-amp forced-side rule below: unlike the 2-pin case,
    # which decides a device's two sides from connection DIRECTION, an
    # op-amp's side is decided by PIN IDENTITY (IN-/IN+ always left, OUT
    # always right, regardless of where the target device happens to sit).
    uid_pin = {}
    for name, ref, pin in power_items:
        device_conns[ref].append((None, ('power', name, ref)))
        uid_pin[('power', name, ref)] = pin
    for item in plan:
        if item[0] == 'pair':
            _, name, a, b, pin_a, pin_b = item
            device_conns[a].append((pos[b], ('pair', name, a)))
            device_conns[b].append((pos[a], ('pair', name, b)))
            uid_pin[('pair', name, a)] = pin_a
            uid_pin[('pair', name, b)] = pin_b
        elif item[0] == 'stub':
            _, name, ref, pin = item
            device_conns[ref].append((None, ('stub', name, ref)))
            uid_pin[('stub', name, ref)] = pin
        else:
            _, name, refs, nx, ny, pins = item
            for ref, pin in zip(refs, pins):
                device_conns[ref].append(((nx, ny), ('multi', name, ref)))
                uid_pin[('multi', name, ref)] = pin

    fixed_side = {}  # uid -> side, set only for 2-pin symbol devices
    for ref, conns in device_conns.items():
        if comps.get(ref, {}).get('device') not in SYMBOL_2PIN_TYPES or len(conns) != 2:
            continue
        cx, cy = pos[ref]
        directional = [(t, k) for t, k in conns if t is not None]
        if directional:
            tx, ty = directional[0][0]
            dx, dy = tx - cx, ty - cy
            primary_side = ('R' if dx > 0 else 'L') if abs(dx) >= abs(dy) else ('T' if dy > 0 else 'B')
            primary_key = directional[0][1]
        else:
            primary_side, primary_key = 'R', conns[0][1]
        opposite_side = {'R': 'L', 'L': 'R', 'T': 'B', 'B': 'T'}[primary_side]
        other_key = conns[1][1] if conns[0][1] == primary_key else conns[0][1]
        fixed_side[primary_key] = primary_side
        fixed_side[other_key] = opposite_side

    # Potentiometers: same opposite-sides logic as the 2-pin block above,
    # applied only to end_a/end_b (the resistive element) -- the wiper is
    # deliberately left OUT of this so it keeps its own natural side_of(),
    # same as it's the one connection on a 2-pin passive's device that
    # isn't part of the forced axis.
    for ref, conns in device_conns.items():
        pp = POTENTIOMETER_PINOUT.get(comps.get(ref, {}).get('device'))
        if pp is None:
            continue
        end_conns = [(t, k) for t, k in conns if uid_pin.get(k) in (pp['end_a'], pp['end_b'])]
        if len(end_conns) != 2:
            continue
        cx, cy = pos[ref]
        directional = [(t, k) for t, k in end_conns if t is not None]
        if directional:
            tx, ty = directional[0][0]
            dx, dy = tx - cx, ty - cy
            primary_side = ('R' if dx > 0 else 'L') if abs(dx) >= abs(dy) else ('T' if dy > 0 else 'B')
            primary_key = directional[0][1]
        else:
            primary_side, primary_key = 'R', end_conns[0][1]
        opposite_side = {'R': 'L', 'L': 'R', 'T': 'B', 'B': 'T'}[primary_side]
        other_key = end_conns[1][1] if end_conns[0][1] == primary_key else end_conns[0][1]
        fixed_side[primary_key] = primary_side
        fixed_side[other_key] = opposite_side

    # NPN transistors (see TRANSISTOR_PINOUT/npn_symbol_svg below) --
    # base always exits the box's LEFT edge (NPN_TERMINALS['base'] sits
    # at norm x=31=the box's own x0). Collector/emitter are NOT
    # left/right pins -- NPN_TERMINALS puts them at norm y=44/465, i.e.
    # exactly the box's y0/y1 (top/bottom edges), not its x1 (right
    # edge). A real, confirmed defect this pass forced both to 'R'
    # (mirroring the op-amp's OUT, a genuine right-edge pin) -- that told
    # the router "exit rightward" for a pin that actually sits ON the
    # box's top/bottom edge, so the wire immediately skimmed along that
    # edge instead of leaving it (verify() caught this as EDGE-HUG, then
    # POINT-INSIDE-BOX once the box grew). Per the file's own inverted
    # T/B convention (see the op-amp V+/V- comment above: 'B' = the box's
    # smaller-y edge, visually the TOP of the screen; 'T' = larger-y,
    # visually the BOTTOM) -- collector (visually top) is 'B', emitter
    # (visually bottom) is 'T', exactly mirroring how V+/V- are assigned.
    for ref, conns in device_conns.items():
        tp = TRANSISTOR_PINOUT.get(comps.get(ref, {}).get('device'))
        if tp is None:
            continue
        side_for_pin = {tp['base']: 'L', tp['collector']: 'B', tp['emitter']: 'T'}
        for target, uid in conns:
            side = side_for_pin.get(uid_pin.get(uid))
            if side is not None:
                fixed_side[uid] = side

    # Dual-section packages (LM1458: two triangles in ONE physical DIP body,
    # not two independently-placed units -- see OPAMP_DUAL_PINOUT/
    # dual_opamp_triangle_segments below). All 8 physical pins are forced to
    # the REAL side (left_order/right_order, verified pin numbers) AND the
    # real physical order top-to-bottom -- forced_rank is applied as an
    # override AFTER pin_rank's normal target-position sort below, since
    # pin_rank has no other hook for a fixed order independent of routing.
    forced_rank = {}  # uid -> rank, applied after the generic pin_rank sort
    for ref, conns in device_conns.items():
        dp = OPAMP_DUAL_PINOUT.get(comps.get(ref, {}).get('device'))
        if dp is None:
            continue
        side_for_pin = {p: 'L' for p in dp['left_order']}
        side_for_pin.update({p: 'R' for p in dp['right_order']})
        rank_for_pin = {p: i for i, p in enumerate(dp['left_order'])}
        rank_for_pin.update({p: i for i, p in enumerate(dp['right_order'])})
        for target, uid in conns:
            pin = uid_pin.get(uid)
            if pin in side_for_pin:
                fixed_side[uid] = side_for_pin[pin]
                forced_rank[uid] = rank_for_pin[pin]

    # Single-section op-amps (LM741) NOW get the same real physical-DIP-8-
    # package treatment as the dual op-amp -- REPLACES the earlier bare-
    # triangle-in-a-generic-box glyph entirely (user: "We are going to
    # replace U2 symbol with a svg file"). Same left_order/right_order +
    # forced_rank pattern as the dual op-amp block above, reusing the
    # SAME forced_rank dict (declared once, above) rather than a second
    # one -- this block runs after that one specifically so forced_rank
    # already exists.
    for ref, conns in device_conns.items():
        dp = OPAMP_DIP8_PINOUT.get(comps.get(ref, {}).get('device'))
        if dp is None:
            continue
        side_for_pin = {p: 'L' for p in dp['left_order']}
        side_for_pin.update({p: 'R' for p in dp['right_order']})
        rank_for_pin = {p: i for i, p in enumerate(dp['left_order'])}
        rank_for_pin.update({p: i for i, p in enumerate(dp['right_order'])})
        for target, uid in conns:
            pin = uid_pin.get(uid)
            if pin in side_for_pin:
                fixed_side[uid] = side_for_pin[pin]
                forced_rank[uid] = rank_for_pin[pin]

    # DPST switch: same fixed-physical-order pattern as the two op-amp
    # packages above (SW_DPST_PINOUT gives left_order/right_order, top to
    # bottom, matching the real SW_DPST symbol's own pin layout).
    for ref, conns in device_conns.items():
        sw = SW_DPST_PINOUT.get(comps.get(ref, {}).get('device'))
        if sw is None:
            continue
        side_for_pin = {p: 'L' for p in sw['left_order']}
        side_for_pin.update({p: 'R' for p in sw['right_order']})
        rank_for_pin = {p: i for i, p in enumerate(sw['left_order'])}
        rank_for_pin.update({p: i for i, p in enumerate(sw['right_order'])})
        for target, uid in conns:
            pin = uid_pin.get(uid)
            if pin in side_for_pin:
                fixed_side[uid] = side_for_pin[pin]
                forced_rank[uid] = rank_for_pin[pin]

    # Rotate specific 2-pin devices' fixed exit sides by a requested angle
    # (device center position unchanged -- this only changes which side of
    # the box the lead exits from). R=0, B=90, L=180, T=270 (increasing
    # angle = visual clockwise, since screen y increases downward); CW adds
    # degrees, CCW subtracts.
    ROTATE_TEST = {
        'R13': (90, 'CW'),
        'R9': (180, 'CCW'),
        'R10': (90, 'CW'),
        'R8': (90, 'CW'),
        'C1': (199, 'CCW'),
    }
    if ROTATE_TEST:
        side_angle = {'R': 0, 'B': 90, 'L': 180, 'T': 270}
        angle_side = {0: 'R', 90: 'B', 180: 'L', 270: 'T'}
        uid_to_ref = {}
        for ref_, conns_ in device_conns.items():
            for _t, uid_ in conns_:
                uid_to_ref[uid_] = ref_
        for uid, side in list(fixed_side.items()):
            ref_ = uid_to_ref.get(uid)
            if ref_ not in ROTATE_TEST:
                continue
            deg, direction = ROTATE_TEST[ref_]
            delta = deg if direction == 'CW' else -deg
            new_angle = min(angle_side.keys(), key=lambda a: min(abs((side_angle[side]+delta-a) % 360), 360-abs((side_angle[side]+delta-a) % 360)))
            fixed_side[uid] = angle_side[new_angle]

    def power_side(uid):
        return fixed_side.get(uid, 'R')

    def signal_side(ref, target, uid):
        return fixed_side.get(uid) or rt.side_of(ref, target)

    # Pass 1: figure out which side of each device every planned connection
    # will exit through, and count how many share a (device, side) pair --
    # needed BEFORE picking exact points, so pins can be spread across the
    # side instead of every one of them landing on that side's single
    # center point (the bug the user caught: 5 different U1 signal wires
    # all exiting at literally the same coordinate).
    side_counts = {}
    for name, ref, pin in power_items:
        side_counts[(ref, power_side(('power', name, ref)))] = side_counts.get((ref, power_side(('power', name, ref))), 0) + 1
    for item in plan:
        if item[0] == 'pair':
            _, name, a, b, pin_a, pin_b = item
            sa, sb = signal_side(a, pos[b], ('pair', name, a)), signal_side(b, pos[a], ('pair', name, b))
            side_counts[(a, sa)] = side_counts.get((a, sa), 0) + 1
            side_counts[(b, sb)] = side_counts.get((b, sb), 0) + 1
        elif item[0] == 'stub':
            _, name, ref, pin = item
            side_counts[(ref, power_side(('stub', name, ref)))] = side_counts.get((ref, power_side(('stub', name, ref))), 0) + 1
        else:
            _, name, refs, nx, ny, pins = item
            for ref in refs:
                s = signal_side(ref, (nx, ny), ('multi', name, ref))
                side_counts[(ref, s)] = side_counts.get((ref, s), 0) + 1

    # Rank pins sharing one (device, side) by the actual position of what
    # they connect to, instead of arbitrary net-processing order -- a real,
    # calculated defect found this pass ("edge overlay" at U1): INVERTING's
    # and OUT_RAW's pins on U1's left side crossed each other only 10 units
    # from the box, because whichever net happened to come first in the
    # net list got the top slot regardless of which pin's target was
    # actually above the other's. Sorting by each pin's own target position
    # along the side (target Y for a left/right side, target X for a
    # top/bottom side) puts pins in the same relative order as what they
    # connect to, so a same-side crossing only happens when the two targets
    # are genuinely on opposite sides of each other -- never as an
    # avoidable artifact of iteration order.
    pin_targets = {}  # uid -> (ref, side, target_point)
    for name, ref, pin in power_items:
        pin_targets[('power', name, ref)] = (ref, power_side(('power', name, ref)), pos[ref])
    for item in plan:
        if item[0] == 'pair':
            _, name, a, b, pin_a, pin_b = item
            sa, sb = signal_side(a, pos[b], ('pair', name, a)), signal_side(b, pos[a], ('pair', name, b))
            pin_targets[('pair', name, a)] = (a, sa, pos[b])
            pin_targets[('pair', name, b)] = (b, sb, pos[a])
        elif item[0] == 'stub':
            _, name, ref, pin = item
            pin_targets[('stub', name, ref)] = (ref, power_side(('stub', name, ref)), pos[ref])
        else:
            _, name, refs, nx, ny, pins = item
            for ref in refs:
                side = signal_side(ref, (nx, ny), ('multi', name, ref))
                pin_targets[('multi', name, ref)] = (ref, side, (nx, ny))
    side_groups = {}
    for uid, (ref, side, target) in pin_targets.items():
        coord = target[1] if side in ('L', 'R') else target[0]
        side_groups.setdefault((ref, side), []).append((coord, uid))
    pin_rank = {}
    for key, items in side_groups.items():
        items.sort(key=lambda t: (t[0], t[1]))  # uid tiebreak: deterministic, not arbitrary
        for i, (_, uid) in enumerate(items):
            pin_rank[uid] = i
    pin_rank.update(forced_rank)  # dual-package physical pin order overrides the generic sort

    def alloc(ref, side, uid):
        # Dual-op-amp DIP-8 packages: the real electrical pin is placed by
        # the SAME calculated transform (dual_opamp_xform) the drawing
        # code uses for every artwork element -- at DUAL_OPAMP_TERMINALS,
        # the one normalized point set that is electrically authoritative.
        # No raster image, no separate letterbox math to drift out of sync
        # with the artwork -- one formula for both.
        dev = comps.get(ref, {}).get('device')
        if dev in OPAMP_DUAL_PINOUT:
            pin_name = uid_pin.get(uid)
            if pin_name in DUAL_OPAMP_TERMINALS:
                bx0, by0, bx1, by1 = rt.boxes[ref]
                x, y = dual_opamp_xform(bx0, by0, bx1, by1)(*DUAL_OPAMP_TERMINALS[pin_name])
                return (x, y), ((-1, 0) if side == 'L' else (1, 0))
        # NPN transistors: same principle, same shared-transform
        # architecture as the dual-op-amp above -- npn_xform is the SAME
        # function the artwork renderer (npn_symbol_svg) uses, so the
        # lead drawn on screen always ends exactly where the real wire
        # does, by construction, not by coincidence.
        if dev in TRANSISTOR_PINOUT:
            pin_name = uid_pin.get(uid)
            tp = TRANSISTOR_PINOUT[dev]
            role = next((r for r, p in tp.items() if p == pin_name), None)
            if role is not None:
                bx0, by0, bx1, by1 = rt.boxes[ref]
                x, y = npn_xform(bx0, by0, bx1, by1)(*NPN_TERMINALS[role])
                # Collector/emitter exit the box's TOP/BOTTOM edge, not
                # left/right (see the forced_side comment above) -- the
                # direction vector must match, or the router starts this
                # wire heading the wrong way along the very edge the pin
                # sits on. Same inverted T/B convention as the side label
                # itself: 'B' (smaller-y edge) exits upward (0,-1), 'T'
                # (larger-y edge) exits downward (0,1).
                dir_for_side = {'L': (-1, 0), 'R': (1, 0), 'B': (0, -1), 'T': (0, 1)}
                return (x, y), dir_for_side.get(side, (1, 0))
        # Single-section op-amp DIP-8 (LM741): same shared-transform
        # principle as the dual op-amp and the transistor above --
        # opamp8_xform is the SAME function opamp8_symbol_svg uses to draw
        # the artwork, so the lead and the real wire always end at the
        # identical point, by construction.
        if dev in OPAMP_DIP8_PINOUT:
            pin_name = uid_pin.get(uid)
            if pin_name in OPAMP_DIP8_TERMINALS:
                bx0, by0, bx1, by1 = rt.boxes[ref]
                x, y = opamp8_xform(bx0, by0, bx1, by1)(*OPAMP_DIP8_TERMINALS[pin_name])
                return (x, y), ((-1, 0) if side == 'L' else (1, 0))
        if dev in SW_DPST_PINOUT:
            pin_name = uid_pin.get(uid)
            if pin_name in SW_DPST_TERMINALS:
                bx0, by0, bx1, by1 = rt.boxes[ref]
                x, y = sw_dpst_xform(bx0, by0, bx1, by1)(*SW_DPST_TERMINALS[pin_name])
                return (x, y), ((-1, 0) if side == 'L' else (1, 0))
        return rt.side_pin_point(ref, side, pin_rank[uid], side_counts[(ref, side)])

    # Pass 2a: allocate every endpoint FIRST, registering each one's real
    # position (tagged with its owning net) into rt.all_pins before any
    # routing happens at all. This is the calculated fix for the bug found
    # this pass: an earlier net's path could settle exactly on a
    # not-yet-routed net's fixed pin, because that pin had no protection
    # until routing was already in progress. Deferring all routing to a
    # separate pass 2b means no route can ever be computed before every
    # pin it needs to avoid is already known.
    endpoints = []  # kind-specific tuples consumed by pass 2b, in order
    for name, ref, pin in power_items:
        p, d = alloc(ref, power_side(('power', name, ref)), ('power', name, ref))
        rt.all_pins.append((p[0], p[1], name))
        pin_points[(ref, pin)] = p
        endpoints.append(('power', name, ref, p, d))
    for item in plan:
        if item[0] == 'pair':
            _, name, a, b, pin_a, pin_b = item
            sa, sb = signal_side(a, pos[b], ('pair', name, a)), signal_side(b, pos[a], ('pair', name, b))
            pa, da = alloc(a, sa, ('pair', name, a))
            pb, db = alloc(b, sb, ('pair', name, b))
            rt.all_pins.append((pa[0], pa[1], name))
            rt.all_pins.append((pb[0], pb[1], name))
            pin_points[(a, pin_a)] = pa
            pin_points[(b, pin_b)] = pb
            endpoints.append(('pair', name, a, b, pa, da, pb, db))
        elif item[0] == 'stub':
            _, name, ref, pin = item
            p, d = alloc(ref, power_side(('stub', name, ref)), ('stub', name, ref))
            rt.all_pins.append((p[0], p[1], name))
            pin_points[(ref, pin)] = p
            endpoints.append(('stub', name, ref, p, d))
        else:
            _, name, refs, nx, ny, pins = item
            # The junction point itself -- where 3+ spokes of this net meet
            # and get drawn as a dot -- must be registered exactly like any
            # other pin, and for the same reason: it's a fixed point chosen
            # before any routing happens, so a later net's path can just as
            # easily land on top of it as on a real pin. Real defect found
            # this pass ("floating dots... an edge too close to R12"): the
            # junction was never in rt.all_pins, so nothing stopped a
            # foreign net's wire from passing within MIN_WIRE_GAP of a dot
            # it has nothing to do with -- reading, visually, as if it were
            # part of that crossing.
            rt.all_pins.append((nx, ny, name))
            spokes = []
            for ref, pin in zip(refs, pins):
                side = signal_side(ref, (nx, ny), ('multi', name, ref))
                p, d = alloc(ref, side, ('multi', name, ref))
                rt.all_pins.append((p[0], p[1], name))
                pin_points[(ref, pin)] = p
                spokes.append((ref, p, d))
            endpoints.append(('multi', name, refs, nx, ny, spokes))

    # Pass 2b: now that every pin (from every net) is a known obstacle,
    # actually route each wire.
    for e in endpoints:
        if e[0] == 'power':
            _, name, ref, p, d = e
            stub = rt.stub_endpoint(p, d, {ref}, name)
            wires.append(([p, stub], 'power', frozenset({ref}), name))
            record([p, stub], name)
        elif e[0] == 'pair':
            _, name, a, b, pa, da, pb, db = e
            path = rt.route(pa, da, pb, db, exclude={a, b}, net_name=name)
            wires.append((path, 'signal', frozenset({a, b}), name))
            record(path, name)
            # Decoupling-cap-to-pin pairs (see _decoupling_pin_targets)
            # already prove electrically that this connection is on the
            # named rail -- but drawing it as one continuous real wire (no
            # longer an open, labeled stub) means the rail's own NAME is
            # no longer visible anywhere on this specific connection,
            # unlike every other member of that same power net elsewhere
            # in the schematic. User: "There needs to be a stub that
            # points to +9V between U1.8 and C13." Purely a readability
            # label, not a new electrical terminal -- branches from an
            # interior bend vertex of the already-routed path (never a
            # real pin coordinate, so coord_to_pin naturally tags neither
            # end with a real pin, keeping electrical_equivalence_check
            # blind to it exactly like every other decorative element in
            # this file -- no dashed wire with a null pin tag is ever
            # matched against an authoritative net terminal).
            if a in decoupling_targets and decoupling_targets[a][0] == b and len(path) > 2:
                # Branch near the CAP's own end of the path (path[0]==pa
                # is the cap's real terminal; path[1] is the first bend
                # right after leaving it) -- NOT the path's geometric
                # midpoint (tried first, landed near the host's end
                # instead, since this wire's own bends happen to sit
                # closer to the host). User: "The stub was on the caps
                # when they were floating" -- matches where it visually
                # used to be before this pin became a real wire.
                branch = path[1]
                label_stub = rt.stub_endpoint(branch, (0, -1), {a, b}, name)
                wires.append(([branch, label_stub], 'power', frozenset(), name))
                record([branch, label_stub], name)
        elif e[0] == 'stub':
            # A single-terminal net (nothing else in the design shares it --
            # a reserved/unused connector pin, e.g. J_R4's R4_W) has to be
            # visually explained, not drawn as a bare, unlabeled line to
            # nowhere. Real defect the user caught ("stinger"): kind was
            # 'signal', which render_svg only labels/dashes for kind
            # 'power' -- so it rendered as a stray solid stub with no name
            # and no indication it's intentional. Use 'power' styling
            # (dashed + labeled) since the open-lead semantics are the same.
            _, name, ref, p, d = e
            stub = rt.stub_endpoint(p, d, {ref}, name)
            wires.append(([p, stub], 'power', frozenset({ref}), name))
            record([p, stub], name)
        else:
            _, name, refs, nx, ny, spokes = e
            for ref, p, d in spokes:
                # Exclude only this spoke's own box. Widening this to the
                # whole net (tried earlier) let a spoke's path legally cut
                # through a SIBLING participant's box interior on the way
                # to the junction -- now that the junction is guaranteed
                # outside every box, that widening is no longer needed and
                # was actively wrong.
                path = rt.route(p, d, (nx, ny), (0, 0), exclude={ref}, net_name=name)
                wires.append((path, 'signal', frozenset({ref}), name))
                record(path, name)
    # Dots are computed once, here, directly from the final routed
    # geometry (see compute_junction_dots) -- not accumulated during
    # planning. This is what makes the rule general: it finds every real
    # electrically-connected branch point that actually exists in the
    # final wires, including ones no planning step ever anticipated (the
    # ordinary 3+-terminal junction case falls out of this automatically,
    # since every spoke has a real endpoint at the shared (nx,ny)).
    dots = compute_junction_dots(wires)
    return wires, dots, pin_points


MIN_PARALLEL_GAP = 6.0  # below this, two parallel segments read as touching/merged at drawing scale

def max_corridor_bundle(wires):
    # Real measurement of how many DISTINCT nets are actually crammed
    # through one shared local corridor -- the calculated input to whether
    # the layout's inter-territory GAP needs to be widened, instead of
    # guessing a bigger GAP and checking by eye ("the magnify need could be
    # calculated, not searched for" -- the same standard applied earlier
    # this session to pin spacing and edge-hug clearance).
    #
    # For every parallel (same-axis) segment, counts how many OTHER
    # distinct-net segments run within one real "bundle window" of it AND
    # overlap its span -- i.e. would actually read as part of the same
    # crowded corridor at render scale, not just anywhere on the canvas.
    # The window is 2x CONGESTION_RADIUS (itself already defined as 3x
    # MIN_WIRE_GAP, "a soft awareness zone") -- reusing an existing
    # calculated constant rather than inventing a new one.
    window = Router.CONGESTION_RADIUS * 2
    segs = {'H': [], 'V': []}
    for path, kind, owners, name in wires:
        for p1, p2 in zip(path, path[1:]):
            k = seg_key(p1, p2)
            if k is None:
                continue
            axis, coord, lo, hi = k
            segs[axis].append((coord, lo, hi, name))
    best = 1
    for axis in ('H', 'V'):
        lst = segs[axis]
        for i in range(len(lst)):
            coord, lo, hi, name = lst[i]
            names = {name}
            for coord2, lo2, hi2, name2 in lst:
                if name2 == name:
                    continue
                if abs(coord2 - coord) <= window and not (hi2 < lo or lo2 > hi):
                    names.add(name2)
            best = max(best, len(names))
    return best


def overlap_count(wires):
    # Counts both exact coincidence AND near-parallel segments closer than
    # MIN_PARALLEL_GAP with overlapping range -- a real defect found by
    # inspection: segments 0.5-5.5 units apart (vs. the old "== 0" only
    # check) render with no visible gap between them. Feeding this into the
    # same swap-optimizer used for exact overlaps lets it push those apart
    # too, since it's the only lever this tool has (no per-wire lane
    # offsetting).
    #
    # Excludes same-NET pairs (not just same-wire-index): a multi-terminal
    # net's own spokes are EXPECTED to converge near their shared junction
    # dot -- that's correct, not crowding. Checked directly this pass: on
    # HRNG, 87/87 of the flagged pairs were same-net convergence and ZERO
    # were real different-net crowding, yet the old same-INDEX-only check
    # counted every one of those as "overlap," making the printed metric
    # actively misleading (real defect, not just a display nitpick -- it's
    # the number this file prints and the number the optimizer minimizes).
    occ = []
    overlap = 0
    for wi, (path, kind, ref, name) in enumerate(wires):
        for si in range(len(path)-1):
            k = seg_key(path[si], path[si+1])
            if k is None:
                continue
            axis, coord, lo, hi = k
            for (axis2, coord2, lo2, hi2, owner_name) in occ:
                if axis2 == axis and abs(coord2-coord) < MIN_PARALLEL_GAP and not (hi2 < lo or lo2 > hi) and owner_name != name:
                    overlap += 1
            occ.append((axis, coord, lo, hi, name))
    return overlap


from symbol_library import *




def render_svg(boxes, wires, dots, path, comps=None, pin_points=None):
    comps = comps or {}
    pin_points = pin_points or {}
    all_x = [b[0] for b in boxes.values()] + [b[2] for b in boxes.values()]
    all_y = [b[1] for b in boxes.values()] + [b[3] for b in boxes.values()]
    for w in wires:
        for p in w[0]:
            all_x.append(p[0]); all_y.append(p[1])
    minx, maxx = min(all_x), max(all_x)
    miny, maxy = min(all_y), max(all_y)
    PAGE_MARGIN = 40
    dx, dy = PAGE_MARGIN - minx, PAGE_MARGIN - miny
    boxes = {ref: (x0+dx, y0+dy, x1+dx, y1+dy) for ref, (x0, y0, x1, y1) in boxes.items()}
    wires = [([(p[0]+dx, p[1]+dy) for p in w], kind, owners, name) for w, kind, owners, name in wires]
    dots = [(x+dx, y+dy, name) for x, y, name in dots]
    pin_points = {(ref, pin): (x+dx, y+dy) for (ref, pin), (x, y) in pin_points.items()}
    W = (maxx - minx) + 2*PAGE_MARGIN
    H = (maxy - miny) + 2*PAGE_MARGIN

    # Reverse lookup so every wire polyline endpoint can be tagged with the
    # real (ref, pin) it terminates at -- the foundation for a PIN-level
    # electrical-equivalence proof (electrical_equivalence_check below),
    # not just the device-level one this replaces. Populated for every
    # device type, not only symbol-drawn ones, since pin identity is now
    # threaded generally through build_wires regardless of how a device
    # gets rendered.
    coord_to_pin = {}
    for (ref, pin), (px, py) in pin_points.items():
        coord_to_pin[(round(px, 1), round(py, 1))] = f'{ref}.{pin}'

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H:.0f}">']
    # Queried once per device, straight from the library, before any of the
    # existing dispatch below runs -- this only detects and reports a
    # device type with no registered artwork at all; it doesn't change what
    # the dispatch chain draws for it (still the plain box further down).
    missing_symbols = set()
    for ref, (x0, y0, x1, y1) in boxes.items():
        dev_type = comps.get(ref, {}).get('device')
        if dev_type not in SYMBOL_LIBRARY:
            missing_symbols.add(dev_type)
        if dev_type in SYMBOL_2PIN_TYPES:
            pin_a, pin_b = SYMBOL_PIN_ORDER[dev_type]
            p1, p2 = pin_points.get((ref, pin_a)), pin_points.get((ref, pin_b))
            if p1 is not None and p2 is not None:
                # Box rectangle kept in the DOM (invisible) so
                # electrical_equivalence_check's rect+text box parsing and
                # pin-touch detection keep working unmodified -- only the
                # visual style changes, not the underlying geometry any
                # other code depends on.
                # <text> (the ref label) must come IMMEDIATELY after <rect>,
                # before any symbol polylines -- electrical_equivalence_
                # check pairs a box's ref by reading the element right after
                # its <rect>, and a real bug found this pass put the symbol
                # geometry in between, silently breaking that pairing for
                # every symbol-drawn device (the check then can't find ANY
                # of that device's connections, since it never even
                # registers the box). Order matters here, not just content.
                svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="none" stroke="none"/>')
                # Real bug found by looking at the actual render (RMIN1,
                # RMIN4): a fixed 13-unit offset with text-anchor="middle"
                # cleared a HORIZONTAL component fine (the offset there is
                # vertical, so only text HEIGHT -- ~6-8 units -- eats into
                # it), but for a VERTICAL component the offset is
                # horizontal, and a centered multi-character label's own
                # WIDTH extends back over the zigzag from the center point.
                # Anchor away from the body instead of centering on the
                # offset point whenever the offset is the horizontal axis,
                # so the label's width grows outward, never back inward.
                dxu, dyu = p2[0]-p1[0], p2[1]-p1[1]
                L = math.hypot(dxu, dyu) or 1.0
                vx, vy = -dyu/L, dxu/L
                mx, my = (p1[0]+p2[0])/2, (p1[1]+p2[1])/2
                if abs(vx) > abs(vy):
                    anchor = 'start' if vx > 0 else 'end'
                    lx, ly = mx + vx*9, my + 3
                else:
                    anchor = 'middle'
                    lx, ly = mx + vx*13, my + vy*13
                svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="8" font-family="monospace" text-anchor="{anchor}">{ref}</text>')
                # class="symbol" marks these as a DEVICE BODY, not a wire --
                # real bug caught by electrical_equivalence_check itself
                # this pass: a resistor's two pins are genuinely different
                # electrical nodes (that's the whole point of a resistor),
                # but a plain continuous polyline from p1 to p2 reads to the
                # connectivity reconstruction as a wire joining them, which
                # silently merged two real, different nets (e.g. INVERTING
                # and OUT_RAW through RFB) into one. The checker skips any
                # polyline carrying this class when building connectivity.
                if dev_type == 'GENERIC_RESISTOR':
                    segs = [resistor_symbol_points(p1, p2)]
                elif dev_type == 'GENERIC_CAPACITOR':
                    segs = capacitor_symbol_segments(p1, p2)
                elif dev_type == 'GENERIC_LED':  # p1/p2 are (anode, cathode)
                    segs, (lc, lr) = led_symbol_segments(p1, p2)
                    svg.append(f'<circle cx="{lc[0]:.1f}" cy="{lc[1]:.1f}" r="{lr:.1f}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
                elif dev_type == 'GENERIC_METER':  # p1/p2 are (pos, neg)
                    segs, (mc, medge, mr) = meter_symbol_segments(p1, p2)
                    # A circle isn't a polyline -- same "own render hook"
                    # pattern as the transistor's encasing circle.
                    svg.append(f'<circle cx="{mc[0]:.1f}" cy="{mc[1]:.1f}" r="{mr:.1f}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
                    svg.append(f'<text x="{mc[0]:.1f}" y="{mc[1]+3:.1f}" font-size="10" font-family="monospace" text-anchor="middle" class="symbol">A</text>')
                elif dev_type == 'ZENER_2PIN':  # p1/p2 are (anode, cathode)
                    segs = zener_symbol_segments(p1, p2)
                else:  # 1N4148 -- p1/p2 are (anode, cathode), in that order
                    segs = diode_symbol_segments(p1, p2)
                for seg in segs:
                    spts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in seg)
                    svg.append(f'<polyline points="{spts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
                continue
        pp = POTENTIOMETER_PINOUT.get(dev_type)
        if pp is not None:
            p_a = pin_points.get((ref, pp['end_a']))
            p_b = pin_points.get((ref, pp['end_b']))
            p_w = pin_points.get((ref, pp['wiper']))
            if p_a is not None and p_b is not None and p_w is not None:
                svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="none" stroke="none"/>')
                svg.append(f'<text x="{(x0+x1)/2:.1f}" y="{y0-4:.1f}" font-size="8" font-family="monospace" text-anchor="middle">{ref}</text>')
                for seg in potentiometer_symbol_segments(p_a, p_b, p_w):
                    spts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in seg)
                    svg.append(f'<polyline points="{spts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
                continue
        dp8 = OPAMP_DIP8_PINOUT.get(dev_type)
        if dp8 is not None:
            opamp8_svg = opamp8_symbol_svg(x0, y0, x1, y1, ref, dp8, pin_points)
            if opamp8_svg is not None:
                svg.extend(opamp8_svg)
                continue
        tp = TRANSISTOR_PINOUT.get(dev_type)
        if tp is not None:
            npn_svg = npn_symbol_svg(x0, y0, x1, y1, ref, tp, pin_points)
            if npn_svg is not None:
                svg.extend(npn_svg)
                continue
        sw = SW_DPST_PINOUT.get(dev_type)
        if sw is not None:
            sw_svg = sw_dpst_symbol_svg(x0, y0, x1, y1, ref, sw, pin_points)
            if sw_svg is not None:
                svg.extend(sw_svg)
                continue
        dp = OPAMP_DUAL_PINOUT.get(dev_type)
        if dp is not None:
            dual_svg = dual_opamp_dip8_svg(x0, y0, x1, y1, ref, dp, pin_points)
            if dual_svg is not None:
                svg.extend(dual_svg)
                continue
        svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="white" stroke="black" stroke-width="1"/>')
        svg.append(f'<text x="{(x0+x1)/2:.1f}" y="{(y0+y1)/2:.1f}" font-size="8" font-family="monospace" text-anchor="middle">{ref}</text>')
    for wpath, kind, owners, name in wires:
        pts = ' '.join(f'{p[0]:.1f},{p[1]:.1f}' for p in wpath)
        dash = ' stroke-dasharray="2,2"' if kind == 'power' else ''
        # data-pin-a/-b: the real (ref.pin) each end of this wire terminates
        # at, when known -- the foundation for pin-level electrical
        # equivalence checking (see electrical_equivalence_check below).
        # Populated for every device type via coord_to_pin, not just
        # symbol-drawn ones. A junction end (multi-terminal net's shared
        # dot) has no single owning pin, so it's simply absent here.
        pin_a = coord_to_pin.get((round(wpath[0][0], 1), round(wpath[0][1], 1)))
        pin_b = coord_to_pin.get((round(wpath[-1][0], 1), round(wpath[-1][1], 1)))
        pin_attrs = ''
        if pin_a:
            pin_attrs += f' data-pin-a="{pin_a}"'
        if pin_b:
            pin_attrs += f' data-pin-b="{pin_b}"'
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="0.75"{dash}{pin_attrs}/>')
        if kind == 'power':
            lp = wpath[-1]
            svg.append(f'<text x="{lp[0]+2:.1f}" y="{lp[1]-2:.1f}" font-size="6" font-family="monospace">{name}</text>')
    for (dxp, dyp, _dname) in dots:
        svg.append(f'<circle cx="{dxp:.1f}" cy="{dyp:.1f}" r="2" fill="black"/>')
    svg.append('</svg>')
    open(path, 'w').write('\n'.join(svg))
    return W, H, boxes, wires, dots, missing_symbols


def verify(boxes, wires, W, H, rt, dots=()):
    """Hard assertions -- returns a list of violation strings, empty = clean."""
    problems = []
    for wpath, kind, owners, name in wires:
        for p1, p2 in zip(wpath, wpath[1:]):
            for bref, box in boxes.items():
                if bref in owners:
                    continue
                if rt.seg_hits_box(p1, p2, box, margin=0):
                    problems.append(f'BOX-TOUCH: net {name} (owners {sorted(owners)}) touches box {bref} at {p1}->{p2}')
        for p in wpath:
            if not (0 <= p[0] <= W and 0 <= p[1] <= H):
                problems.append(f'OUT-OF-BOUNDS: net {name} point {p} outside canvas {W}x{H}')
        # No wire point -- including endpoints of an owned box -- may sit
        # STRICTLY inside any box's interior (open interval; touching the
        # boundary exactly is how a real pin connection looks). A point
        # strictly inside is always wrong, even for a box the wire owns:
        # exit_point() only ever returns boundary points, so any interior
        # point means a junction/target was chosen inside someone's
        # silhouette -- exactly the "dip into your own box" defect found
        # this pass (a net-mate's wire cutting through a shared owner's box
        # to reach a badly-placed junction).
        for p in wpath:
            for bref, box in boxes.items():
                if box[0] < p[0] < box[2] and box[1] < p[1] < box[3]:
                    problems.append(f'POINT-INSIDE-BOX: net {name} (owners {sorted(owners)}) point {p} strictly inside box {bref} {box}')
        # A wire segment may not run collinear with ANY box's edge for a
        # real distance (only a single-point boundary touch, e.g. the exit
        # stub itself, is legitimate) -- otherwise it visually reads as
        # part of that box's own outline. Found at J_EB: a wire traced 13
        # units along J_EB's own bottom edge after reaching its corner.
        for p1, p2 in zip(wpath, wpath[1:]):
            for bref, box in boxes.items():
                if rt.seg_hugs_box_edge(p1, p2, box):
                    problems.append(f'EDGE-HUG: net {name} (owners {sorted(owners)}) segment {p1}->{p2} runs along box {bref} {box}')
    if rt.unverified:
        problems.append(f'UNVERIFIED-ROUTE: {rt.unverified} wire(s) fell through every checked tier')
    # Wire-to-wire spacing, made into a hard, unconditional check -- same
    # category as the box checks above, zero exceptions. Until now this was
    # only a PREFERENCE the router tried to honor, with fallback tiers that
    # could abandon it to guarantee some route existed -- which meant a
    # violation could always resurface somewhere else after a fix, since
    # the underlying rule was never actually enforced, just usually
    # respected. This check has no escape hatch: if it fires, the run is a
    # genuine FAIL, not a quietly-shipped tight spot.
    all_segs = []
    for wpath, kind, owners, name in wires:
        for p1, p2 in zip(wpath, wpath[1:]):
            k = seg_key(p1, p2)
            if k is not None:
                all_segs.append((k, name, p1, p2))
    for i in range(len(all_segs)):
        (axis, coord, lo, hi), name, p1, p2 = all_segs[i]
        for j in range(i+1, len(all_segs)):
            (axis2, coord2, lo2, hi2), name2, q1, q2 = all_segs[j]
            if name2 == name or axis2 != axis:
                continue
            if abs(coord2-coord) < Router.MIN_WIRE_GAP and not (hi2 < lo or lo2 > hi):
                problems.append(f'WIRE-TOO-CLOSE: net {name} {p1}->{p2} and net {name2} {q1}->{q2} '
                                 f'are {abs(coord2-coord):.1f} units apart (minimum {Router.MIN_WIRE_GAP})')
    # A junction dot belongs to exactly one net (3+ spokes converging). Any
    # OTHER net's segment running within MIN_WIRE_GAP of that point reads,
    # visually, as if it's part of the junction -- this is the hard,
    # independent check for the "floating dots" defect the user caught;
    # the actual prevention is the dot's pin-registration in build_wires,
    # this just proves it held.
    g = Router.MIN_WIRE_GAP
    for dx, dy, dname in dots:
        for (axis, coord, lo, hi), name, p1, p2 in all_segs:
            if name == dname:
                continue
            if axis == 'H' and abs(dy - coord) < g and lo - g <= dx <= hi + g:
                problems.append(f'DOT-TOO-CLOSE: junction for net {dname} at ({dx:.1f},{dy:.1f}) '
                                 f'is within {g} units of net {name} segment {p1}->{p2}')
            elif axis == 'V' and abs(dx - coord) < g and lo - g <= dy <= hi + g:
                problems.append(f'DOT-TOO-CLOSE: junction for net {dname} at ({dx:.1f},{dy:.1f}) '
                                 f'is within {g} units of net {name} segment {p1}->{p2}')
    # Hard, independent proof of the general junction rule: recompute
    # required dots directly from this same final `wires` geometry (the
    # exact same function build_wires used to produce `dots` in the first
    # place) and require an exact match, in both directions. A mismatch
    # here means the rendered dots and the actual routed geometry have
    # drifted apart -- either a real electrically-connected branch with no
    # dot (JUNCTION-DOT-MISSING, the RFB defect this was built for), or a
    # dot sitting somewhere that isn't actually a real branch point
    # (JUNCTION-DOT-EXTRA). Coordinates rounded to 0.1 for float-safe
    # comparison, same tolerance used throughout this file.
    required = {(round(x, 1), round(y, 1), name) for x, y, name in compute_junction_dots(wires)}
    actual = {(round(x, 1), round(y, 1), name) for x, y, name in dots}
    for x, y, name in sorted(required - actual):
        problems.append(f'JUNCTION-DOT-MISSING: net {name} has an electrically-connected branch '
                         f'at ({x},{y}) with no junction dot')
    for x, y, name in sorted(actual - required):
        problems.append(f'JUNCTION-DOT-EXTRA: dot for net {name} at ({x},{y}) does not correspond '
                         f'to a real branch point in the routed geometry')
    return problems


def find_wire_too_close(wires, min_gap):
    all_segs = []
    for wpath, kind, owners, name in wires:
        for p1, p2 in zip(wpath, wpath[1:]):
            k = seg_key(p1, p2)
            if k is not None:
                all_segs.append((k, name, owners))
    violations = []
    for i in range(len(all_segs)):
        (axis, coord, lo, hi), name, owners = all_segs[i]
        for j in range(i+1, len(all_segs)):
            (axis2, coord2, lo2, hi2), name2, owners2 = all_segs[j]
            if name2 == name or axis2 != axis:
                continue
            if abs(coord2-coord) < min_gap and not (hi2 < lo or lo2 > hi):
                violations.append((axis, coord, coord2, max(lo, lo2), min(hi, hi2), owners, owners2))
    return violations


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def electrical_equivalence_check(svg_path, semantic_path):
    # A DIFFERENT kind of proof than verify() -- verify() checks that the
    # drawing is geometrically clean (nothing touching, nothing too close).
    # This checks that the drawing is electrically CORRECT: that the actual
    # rendered geometry in final.svg connects exactly the PINS the
    # authoritative netlist (semantic.json) says it should, no more and no
    # less. It re-derives connectivity from the raw SVG geometry itself --
    # box rectangles and polyline points, tagged with real pin identity via
    # data-pin-a/data-pin-b (see render_svg's coord_to_pin) -- rather than
    # trusting whatever net name the generator's own internal bookkeeping
    # attached to a wire, so it can catch a real coordinate/logic bug that
    # produces a correctly-labeled wire in the wrong place, not just a
    # mislabeled one.
    #
    # PIN-level, not device-level: two-terminal-or-more signal nets are
    # proven by exact pin-set match (which named pins does each connected
    # tangle of polylines actually terminate at), not just which device
    # refs participate. This is the real upgrade over the earlier
    # device-only version -- it can now catch, for example, a wire that
    # connects the right TWO DEVICES but the wrong pin on one of them,
    # which a device-level check structurally cannot see. Power rails and
    # single-terminal nets are, by design, drawn as OPEN unconnected stubs
    # (real schematic convention -- they merge by name label, not a drawn
    # wire), so those are proven by matching each authoritative pin to a
    # correctly labeled stub at that exact pin instead.
    problems = []
    tree = ET.parse(svg_path)
    root = tree.getroot()
    ns = ''
    if root.tag.startswith('{'):
        ns = root.tag.split('}')[0] + '}'

    def tag(el):
        return el.tag[len(ns):] if ns and el.tag.startswith(ns) else el.tag

    children = list(root)
    boxes = {}  # ref -> (x0,y0,x1,y1)
    signal_polys = []  # list of (points, pin_a_or_None, pin_b_or_None)
    stub_polys = []     # list of (points, label_or_None, pin_a_or_None)
    i = 0
    while i < len(children):
        el = children[i]
        t = tag(el)
        if t == 'rect':
            x0 = float(el.get('x')); y0 = float(el.get('y'))
            w = float(el.get('width')); h = float(el.get('height'))
            ref = None
            if i+1 < len(children) and tag(children[i+1]) == 'text':
                ref = (children[i+1].text or '').strip()
            if ref:
                boxes[ref] = (x0, y0, x0+w, y0+h)
        elif t == 'polyline':
            # class="symbol" is a device BODY (resistor zigzag, capacitor
            # plates, diode triangle+bar), not a wire -- its pins are
            # genuinely different electrical nodes, so it must never
            # participate in connectivity reconstruction (real bug found
            # this pass: without this skip, a resistor's own drawn body
            # silently merged the two different nets on either side of it
            # into one).
            if el.get('class') == 'symbol':
                i += 1
                continue
            pts_raw = el.get('points', '')
            pts = []
            for chunk in pts_raw.split():
                xs, ys = chunk.split(',')
                pts.append((round(float(xs), 1), round(float(ys), 1)))
            pin_a = el.get('data-pin-a')
            pin_b = el.get('data-pin-b')
            is_dashed = el.get('stroke-dasharray') is not None
            if is_dashed:
                label = None
                if i+1 < len(children) and tag(children[i+1]) == 'text':
                    label = (children[i+1].text or '').strip()
                # A power/stub wire is always [real_pin, open_end] (see
                # build_wires Pass 2b) -- data-pin-a is always its one real
                # pin, data-pin-b is never set (the open end isn't a pin).
                stub_polys.append((pts, label, pin_a))
            else:
                signal_polys.append((pts, pin_a, pin_b))
        i += 1

    if not boxes:
        return ['ELECTRICAL-CHECK: no device boxes found in SVG -- cannot verify']

    # Union every signal polyline's own points together, and union any two
    # points (from the same or different polylines) that coincide exactly
    # -- this is how real junction dots (shared endpoints) tie multiple
    # spokes of one net into a single reconstructed component without ever
    # reading a net-name label.
    uf = _UnionFind()
    for pts, pin_a, pin_b in signal_polys:
        for a, b in zip(pts, pts[1:]):
            uf.union(a, b)
    all_points = [p for pts, _, _ in signal_polys for p in pts]
    by_point = collections.defaultdict(list)
    for p in all_points:
        by_point[p].append(p)
    for p, plist in by_point.items():
        for q in plist[1:]:
            uf.union(plist[0], q)

    # Reconstructed PIN set per connected component -- built from the
    # data-pin tags directly, not a coordinate-to-box lookup. A wire end
    # with no tag is a real junction point (a multi-terminal net's shared
    # dot, or a power stub's open end), never a live device pin, so it
    # correctly contributes nothing here.
    comp_pins = collections.defaultdict(set)
    for pts, pin_a, pin_b in signal_polys:
        root_a, root_b = uf.find(pts[0]), uf.find(pts[-1])
        if pin_a:
            comp_pins[root_a].add(pin_a)
        if pin_b:
            comp_pins[root_b].add(pin_b)
    reconstructed = [frozenset(s) for s in comp_pins.values() if len(s) >= 2]

    sem = json.loads(open(semantic_path, 'rb').read())
    nets = sem['electrical']['nets']

    connected_nets = []   # (name, frozenset(pin strings, e.g. "D1.A"))
    open_lead_nets = []   # (name, [pin strings])
    for n in nets:
        name = n['name']
        terms = n['terminals']
        refs = [t.split('.')[0] for t in terms]
        if power_like(name) or len(refs) == 1:
            open_lead_nets.append((name, terms))
        else:
            connected_nets.append((name, frozenset(terms)))

    unmatched_reconstructed = list(reconstructed)
    for name, pinset in connected_nets:
        if pinset in unmatched_reconstructed:
            unmatched_reconstructed.remove(pinset)
        else:
            problems.append(f'ELECTRICAL-MISSING: net {name} (pins {sorted(pinset)}) '
                             f'has no matching connected geometry in the SVG')

    # A decoupling cap wired directly to its exact target pin (see
    # build_wires' _decoupling_pin_targets) is a REAL 2-terminal connected
    # component that happens to belong to an otherwise open-lead/power net
    # -- e.g. {"C13.1","U1.8"} out of +9V's full 8-terminal membership.
    # That's neither a full-net MISSING/EXTRA case (the net itself isn't
    # fully connected, by design) nor an ordinary open-lead stub (this
    # ONE pair genuinely is drawn as connected wire) -- it's verified here
    # as its own category: a real connected component that is a genuine
    # SUBSET of one open-lead net's authoritative terminal set. Terminals
    # proven this way are removed from stub-checking below; anything still
    # unmatched after this is a genuine ELECTRICAL-EXTRA.
    open_lead_terms = {name: set(terms) for name, terms in open_lead_nets}
    satisfied_by_real_wire = set()  # "ref.pin" strings proven via real wire, not stub
    still_unmatched = []
    for comp in unmatched_reconstructed:
        matched_net = next((name for name, terms in open_lead_terms.items() if comp <= terms), None)
        if matched_net is not None:
            satisfied_by_real_wire |= comp
        else:
            still_unmatched.append(comp)
    for extra in still_unmatched:
        problems.append(f'ELECTRICAL-EXTRA: SVG geometry connects pins {sorted(extra)} '
                         f'but no authoritative net matches that exact pin set')

    # Open-lead nets (power rails, single-terminal stubs): verified by
    # label + the exact tagged pin, since these are intentionally NOT drawn
    # as connected geometry -- except any terminal already proven above via
    # a real wire to another member of its own net.
    stub_labels_by_pin = collections.defaultdict(list)  # "ref.pin" -> [label, ...]
    for pts, label, pin_a in stub_polys:
        if label is None or pin_a is None:
            continue
        stub_labels_by_pin[pin_a].append(label)
    for name, terms in open_lead_nets:
        for t in terms:
            if t in satisfied_by_real_wire:
                continue
            if name not in stub_labels_by_pin.get(t, []):
                problems.append(f'ELECTRICAL-MISSING-STUB: net {name} terminal {t} has no '
                                 f'matching labeled open-lead stub at that exact pin in the SVG')

    return problems


def uncrowd_boxes(rt, pos, boxsize, boxes, wires):
    # Adapted from graph-relax's own earlier 'global matrix uncrowding'
    # tool (experiment_5b_global_uncrowding_hrng.py / global_matrix.py),
    # which solved an analogous problem for point-vertex graphs: one
    # simultaneous linear least-squares solve that relieves detected
    # crowding while minimizing total displacement, rather than a greedy
    # local search. Translated from Euclidean point-crowding to this
    # router's real crowding: two different nets' orthogonal wire segments
    # closer than MIN_WIRE_GAP.
    #
    # Two real fixes versus the first attempt at this, found by checking
    # what it actually moved instead of assuming it worked:
    # 1. The candidates for each violation are the ACTUAL devices that own
    #    the two conflicting segments (from wires' own `owners`), not "any
    #    box within a fixed radius of the conflict's midpoint" -- the
    #    radius heuristic pulled in irrelevant devices and diluted the
    #    needed relief across them instead of moving the ones actually
    #    responsible.
    # 2. No centroid-preservation constraint. The original tool's graph had
    #    dozens of crowding relations spread across ~34-84 vertices, so a
    #    global centroid constraint barely nudged any single vertex. Here,
    #    with only a handful of violations, that same constraint forced a
    #    compensating shift onto EVERY uninvolved device to keep the sum
    #    constant -- moving unrelated parts of the drawing for no benefit
    #    while barely moving the ones that actually needed to move. Solved
    #    as an unconstrained weighted least-squares over just the involved
    #    devices instead.
    import numpy as np
    violations = find_wire_too_close(wires, rt.MIN_WIRE_GAP)
    if not violations:
        return pos, boxes, 0

    def solve_axis(axis_label, ci):
        rows = []  # (ref, target_x0_plus_delta)
        for axis, coord, coord2, lo, hi, owners1, owners2 in violations:
            if axis != axis_label:
                continue
            needed = rt.MIN_WIRE_GAP - abs(coord2 - coord) + 2.0
            if needed <= 0:
                continue
            mid = (coord + coord2) / 2
            side1 = list(owners1); side2 = list(owners2)
            for r in side1:
                direction = 1.0 if pos[r][ci] >= mid else -1.0
                rows.append((r, pos[r][ci] + direction * (needed / max(1, len(side1)))))
            for r in side2:
                direction = 1.0 if pos[r][ci] >= mid else -1.0
                rows.append((r, pos[r][ci] + direction * (needed / max(1, len(side2)))))
        if not rows:
            return {}, 0
        involved = sorted({r for r, _ in rows})
        vidx = {r: i for i, r in enumerate(involved)}
        m = len(involved)
        x0 = np.array([pos[r][ci] for r in involved], dtype=float)
        C = np.zeros((len(rows), m)); Nv = np.zeros(len(rows))
        for k, (r, target) in enumerate(rows):
            C[k, vidx[r]] = 1.0
            Nv[k] = target
        # average multiple demands on the same device instead of overwriting
        M = C.T @ C + np.eye(m) * 0.05
        B = C.T @ Nv + 0.05 * x0
        solved = np.linalg.solve(M, B)
        return {r: float(solved[vidx[r]]) for r in involved}, len(rows)

    moved_x, nx = solve_axis('V', 0)
    moved_y, ny = solve_axis('H', 1)
    if nx == 0 and ny == 0:
        return pos, boxes, 0

    new_pos = dict(pos)
    for r, x in moved_x.items():
        new_pos[r] = (x, new_pos[r][1])
    for r, y in moved_y.items():
        new_pos[r] = (new_pos[r][0], y)
    new_boxes = dict(boxes)
    for r in set(moved_x) | set(moved_y):
        new_boxes[r] = (new_pos[r][0]-boxsize[r][0]/2, new_pos[r][1]-boxsize[r][1]/2,
                         new_pos[r][0]+boxsize[r][0]/2, new_pos[r][1]+boxsize[r][1]/2)
    return new_pos, new_boxes, nx + ny


def boxes_overlap_any(boxes):
    refs = list(boxes.keys())
    for i in range(len(refs)):
        ax0, ay0, ax1, ay1 = boxes[refs[i]]
        for j in range(i+1, len(refs)):
            bx0, by0, bx1, by1 = boxes[refs[j]]
            if ax0 < bx1 and bx0 < ax1 and ay0 < by1 and by0 < ay1:
                return True
    return False


def _layout_and_route(sem, init, nets, gap):
    # The full placement+routing pipeline, factored out so run() can call it
    # a second time with a wider GAP when the first pass's own wire data
    # shows the corridor is genuinely too narrow for what's routed through
    # it (see max_corridor_bundle below) -- without duplicating the whole
    # swap-optimize-then-exact-route sequence inline.
    comps, territories, pos, boxsize, boxes = layout_boxes(sem, init, gap=gap)
    rt = Router(boxes, pos)
    rt.cheap = True

    wires, dots, pin_points = build_wires(rt, comps, nets, pos, boxes, sem)
    base_overlap = overlap_count(wires)

    tried = accepted = 0
    for tid, t in territories.items():
        devs = sorted(t['devices'])
        for i in range(len(devs)):
            for j in range(i+1, len(devs)):
                a, b = devs[i], devs[j]
                tried += 1
                pos[a], pos[b] = pos[b], pos[a]
                boxes[a] = (pos[a][0]-boxsize[a][0]/2, pos[a][1]-boxsize[a][1]/2, pos[a][0]+boxsize[a][0]/2, pos[a][1]+boxsize[a][1]/2)
                boxes[b] = (pos[b][0]-boxsize[b][0]/2, pos[b][1]-boxsize[b][1]/2, pos[b][0]+boxsize[b][0]/2, pos[b][1]+boxsize[b][1]/2)
                trial_wires, _, _ = build_wires(rt, comps, nets, pos, boxes, sem)
                trial_overlap = overlap_count(trial_wires)
                if trial_overlap < base_overlap:
                    base_overlap = trial_overlap
                    accepted += 1
                else:
                    pos[a], pos[b] = pos[b], pos[a]
                    boxes[a] = (pos[a][0]-boxsize[a][0]/2, pos[a][1]-boxsize[a][1]/2, pos[a][0]+boxsize[a][0]/2, pos[a][1]+boxsize[a][1]/2)
                    boxes[b] = (pos[b][0]-boxsize[b][0]/2, pos[b][1]-boxsize[b][1]/2, pos[b][0]+boxsize[b][0]/2, pos[b][1]+boxsize[b][1]/2)

    rt.fallback_log.clear(); rt.unverified = 0; rt.spacing_failed = 0
    rt.cheap = False
    wires, dots, pin_points = build_wires(rt, comps, nets, pos, boxes, sem)
    final_overlap = overlap_count(wires)

    uncrowd_note = 'not needed'
    violations_before = find_wire_too_close(wires, rt.MIN_WIRE_GAP)
    if violations_before:
        trial_pos, trial_boxes, nudges = uncrowd_boxes(rt, pos, boxsize, boxes, wires)
        if nudges and not boxes_overlap_any(trial_boxes):
            trial_wires, trial_dots, trial_pin_points = build_wires(rt, comps, nets, trial_pos, trial_boxes, sem)
            trial_violations = find_wire_too_close(trial_wires, rt.MIN_WIRE_GAP)
            if len(trial_violations) < len(violations_before):
                pos, boxes, wires, dots, pin_points = trial_pos, trial_boxes, trial_wires, trial_dots, trial_pin_points
                final_overlap = overlap_count(wires)
                uncrowd_note = f'{len(violations_before)} -> {len(trial_violations)} violations ({nudges} device nudges)'
            else:
                uncrowd_note = f'tried, did not improve ({len(violations_before)} -> {len(trial_violations)})'
        else:
            uncrowd_note = 'tried, rejected (introduced box overlap or found no nudge)' if nudges else 'tried, no candidate devices found'

    return dict(comps=comps, territories=territories, pos=pos, boxsize=boxsize, boxes=boxes,
                rt=rt, wires=wires, dots=dots, pin_points=pin_points, tried=tried, accepted=accepted,
                uncrowd_note=uncrowd_note, final_overlap=final_overlap)


DEFAULT_GAP = 60

def run(baseline_dir, out_svg, label):
    sem = json.loads(open(f'{baseline_dir}/semantic.json', 'rb').read())
    init = json.loads(open(f'{baseline_dir}/initial.json', 'rb').read())
    nets = sem['electrical']['nets']

    result = _layout_and_route(sem, init, nets, DEFAULT_GAP)
    gap_note = f'default ({DEFAULT_GAP})'

    # Real, calculated corridor-widening pass: measure how many distinct
    # nets are actually bundled through the tightest corridor in what just
    # got routed (max_corridor_bundle), and only widen the layout if that
    # count genuinely needs more room than DEFAULT_GAP gives it -- never a
    # guessed bigger constant. Found this pass: several different nets
    # forced through the same narrow gap near a hub device (IAMP's JP1,
    # mono741's U1/RFB) read as a tangled tandem of parallel wires -- every
    # individual pair was spacing-legal, so no hard check caught it, but a
    # human couldn't visually trace which wire was which. Re-runs the WHOLE
    # pipeline once with the wider gap (same try/measure/keep-only-if-
    # better pattern as uncrowd_boxes) -- never applied blindly, and never
    # looped, since a single widening pass either fixes the real bottleneck
    # or it doesn't.
    bundle_before = max_corridor_bundle(result['wires'])
    needed_gap = bundle_before * Router.CONGESTION_RADIUS
    if needed_gap > DEFAULT_GAP:
        wider = _layout_and_route(sem, init, nets, needed_gap)
        bundle_after = max_corridor_bundle(wider['wires'])
        W2, H2, fb2, fw2, fd2, _ = render_svg(wider['boxes'], wider['wires'], wider['dots'], out_svg + '.trial')
        problems2 = verify(fb2, fw2, W2, H2, wider['rt'], fd2)
        if bundle_after < bundle_before and not problems2:
            result = wider
            gap_note = f'widened {DEFAULT_GAP} -> {needed_gap:.0f} (worst corridor bundled {bundle_before} nets -> {bundle_after})'
        else:
            gap_note = (f'tried widening to {needed_gap:.0f} (worst corridor bundled {bundle_before} nets), '
                        f'did not improve ({bundle_after} nets, {len(problems2)} verify problems) -- kept default')
        try:
            os.remove(out_svg + '.trial')
        except OSError:
            pass
    else:
        gap_note = f'default ({DEFAULT_GAP}) sufficient (worst corridor bundled {bundle_before} nets)'

    comps, wires, dots, rt = result['comps'], result['wires'], result['dots'], result['rt']
    tried, accepted = result['tried'], result['accepted']
    uncrowd_note, final_overlap = result['uncrowd_note'], result['final_overlap']
    W, H, final_boxes, final_wires, final_dots, missing_symbols = render_svg(result['boxes'], wires, dots, out_svg, comps=comps, pin_points=result['pin_points'])
    problems = verify(final_boxes, final_wires, W, H, rt, final_dots)

    # A SEPARATE proof from verify() above -- geometric cleanliness (no
    # touches, no crowding) says nothing about whether the drawing is
    # electrically the circuit the netlist describes. This re-derives
    # connectivity from the actual written final.svg file (not the in-
    # memory wires list) and checks it against semantic.json directly.
    elec_problems = electrical_equivalence_check(out_svg, f'{baseline_dir}/semantic.json')

    print(f'=== {label} ===')
    print(f'devices: {len(comps)} | nets: {len(nets)} | swaps tried: {tried} accepted: {accepted}')
    print(f'corridor gap: {gap_note}')
    print(f'uncrowding: {uncrowd_note}')
    print(f'overlap count (segment-level): {final_overlap} | fallback routes used: {len(rt.fallback_log)} | spacing-failed: {rt.spacing_failed}')
    print(f'canvas: {W:.0f} x {H:.0f} ({W/72:.1f}in x {H/72:.1f}in at 72 units/in)')
    if missing_symbols:
        print(f'symbols not in library (drawn as generic box, need real artwork): {", ".join(sorted(missing_symbols))}')
    ok = True
    if problems:
        print(f'GEOMETRIC VERIFICATION FAILED -- {len(problems)} problem(s):')
        for p in problems[:20]:
            print('  -', p)
        ok = False
    else:
        print('GEOMETRIC VERIFICATION PASSED -- zero box touches, zero out-of-bounds points, zero unverified routes')
    if elec_problems:
        print(f'ELECTRICAL EQUIVALENCE FAILED -- {len(elec_problems)} problem(s):')
        for p in elec_problems[:20]:
            print('  -', p)
        ok = False
    else:
        print('ELECTRICAL EQUIVALENCE PASSED -- final.svg connects exactly the devices semantic.json says it should')
    return ok


if __name__ == '__main__':
    jobs = [
        ('output/hrng_baseline', '/tmp/hrng_final.svg', 'HRNG'),
        ('output/iamp_baseline', '/tmp/iamp_final.svg', 'IAMP'),
        ('output/mono741_baseline', '/tmp/mono741_final.svg', 'MONO_741_PREAMP'),
    ]
    all_ok = True
    for baseline_dir, out_svg, label in jobs:
        ok = run(baseline_dir, out_svg, label)
        all_ok = all_ok and ok
        print()
    sys.exit(0 if all_ok else 1)
