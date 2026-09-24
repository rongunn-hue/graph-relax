#!/usr/bin/env python3
"""General-purpose, standalone translator: any .circuit netlist file -> a KiCad .kicad_pcb (footprints + real
net connectivity), with no .kicad_sch involved anywhere. Parses the .circuit file directly - no dependency on
any project-specific pre-built semantic.json.

Must run under KiCad's own bundled Python (this is where the real `pcbnew` module and real footprint libraries
live), e.g.:
    flatpak run --command=python3 org.kicad.KiCad xlate_circuit_to_pcb.py <circuit> <output.kicad_pcb> [options]

What this script does NOT do, on purpose (see project discussion): decide PCB placement (every footprint lands
on a neutral alphabetical filler grid - zero information from any schematic placer), or touch KiCad's schematic
editor/IPC at all.

Device-type library (FOOTPRINTS/PAD_MAP below) covers whatever device types have been seen so far. Running this
on a .circuit file with a new device type fails LOUDLY before touching the board (see check_coverage()) with a
clear list of what's missing, rather than guessing or silently skipping pads.
"""
import argparse
import json
import re
import sys
import uuid

import pcbnew

LIBDIR_DEFAULT = "/app/extensions/Library/footprints"

# device type -> (footprint library .pretty dir, footprint name)
# GENERIC_CAPACITOR is one uniform non-polarized THT ceramic-disc footprint by default - appropriate for AC
# coupling (low leakage matters) and IC-local bypass (low ESR matters). A specific capacitor that should be a
# real polarized electrolytic instead (bulk reservoir duty, larger values) is a per-board judgment call, not a
# device-type fact - see --overrides, not this table.
FOOTPRINTS = {
    "GENERIC_RESISTOR":    ("Resistor_THT.pretty",               "R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal"),
    "GENERIC_CAPACITOR":   ("Capacitor_THT.pretty",               "C_Disc_D5.0mm_W2.5mm_P5.00mm"),
    "ZENER_2PIN":          ("Diode_THT.pretty",                   "D_DO-35_SOD27_P7.62mm_Horizontal"),
    "1N4148":              ("Diode_THT.pretty",                   "D_DO-35_SOD27_P7.62mm_Horizontal"),
    "GENERIC_CONNECTOR_2": ("Connector_PinHeader_2.54mm.pretty",  "PinHeader_1x02_P2.54mm_Vertical"),
    "GENERIC_CONNECTOR_3": ("Connector_PinHeader_2.54mm.pretty",  "PinHeader_1x03_P2.54mm_Vertical"),
    "P2N2222A":            ("Package_TO_SOT_THT.pretty",          "TO-92_Inline"),
    "WS78L05":             ("Package_TO_SOT_THT.pretty",          "TO-92_Inline"),
    "LM1458":              ("Package_DIP.pretty",                 "DIP-8_W7.62mm"),
    "LM741":               ("Package_DIP.pretty",                 "DIP-8_W7.62mm"),
    "PICAXE_08M2":         ("Package_DIP.pretty",                 "DIP-8_W7.62mm"),
}

# device type -> {.circuit pin name: real footprint pad number}. Identity everywhere the .circuit file's own
# pin label already IS the real physical pin number. Two verified-against-real-datasheet exceptions baked in
# here as real physical facts about these device types (not board-specific judgment calls, so they belong in
# this shared library, not in --overrides):
#   - P2N2222A (TO-92): 1=C,2=B,3=E, confirmed against the ON Semi datasheet's own pinout diagram.
#   - WS78L05 (TO-92): 1=OUT,2=GND,3=IN, confirmed against the Wing Shing datasheet's own TO-92 lead diagram.
#   - Any diode/zener on this footprint: the real D_DO-35_SOD27_P7.62mm_Horizontal.kicad_mod file's own "K"
#     silkscreen label sits at pad 1's position, i.e. pad "1" = cathode - so A->pad2, K->pad1, not identity.
PAD_MAP = {
    "GENERIC_RESISTOR":    {"1": "1", "2": "2"},
    "GENERIC_CAPACITOR":   {"1": "1", "2": "2"},
    "ZENER_2PIN":          {"A": "2", "K": "1"},
    "1N4148":              {"A": "2", "K": "1"},
    "GENERIC_CONNECTOR_2": {"1": "1", "2": "2"},
    "GENERIC_CONNECTOR_3": {"1": "1", "2": "2", "3": "3"},
    "P2N2222A":            {"1": "1", "2": "2", "3": "3"},
    "WS78L05":             {"1": "1", "2": "2", "3": "3"},
    "LM1458":              {str(n): str(n) for n in range(1, 9)},
    "LM741":               {str(n): str(n) for n in range(1, 9)},
    "PICAXE_08M2":         {str(n): str(n) for n in range(1, 9)},
}

UUID_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "circuit-to-pcb-xlate")


def parse_circuit(path):
    """Direct parse of the .circuit text format - no dependency on any pre-built semantic.json. Only 3 line
    kinds matter for PCB generation: `component REF TYPE`, `net NAME TERM...`, `nc TERM`. Everything else
    (circuit title, placement hints) is real data for OTHER purposes (schematic placement, provenance) but
    irrelevant to what a PCB needs, so it's simply skipped here, not misparsed. Net names are sometimes quoted
    (e.g. `net "+9V" ...`) - quotes are stripped, not left embedded in the net name."""
    components, nets, nc = {}, [], set()
    for raw in open(path).read().splitlines():
        parts = raw.split()
        if not parts:
            continue
        kw = parts[0]
        if kw == "component" and len(parts) == 3:
            components[parts[1]] = parts[2]
        elif kw == "net" and len(parts) >= 3:
            nets.append({"name": parts[1].strip('"'), "terminals": parts[2:]})
        elif kw == "nc" and len(parts) == 2:
            nc.add(parts[1])
    return components, nets, nc


def parse_values(path):
    """Real per-ref component values/descriptions (e.g. "10k", "100nF", "LM1458", "MORSE KEY"), straight from
    the .circuit file's own `source_metadata '{...}'` line when present - never invented. Returns {} if the
    file has no source_metadata (older/simpler .circuit files), so callers must fall back sanely."""
    for raw in open(path).read().splitlines():
        if raw.startswith("source_metadata "):
            payload = raw[len("source_metadata "):].strip()
            if payload[:1] == "'" and payload[-1:] == "'":
                payload = payload[1:-1]
            meta = json.loads(payload)
            parts = meta.get("normalized", {}).get("parts", [])
            return {p["ref"]: p["value"] for p in parts if "ref" in p and "value" in p}
    return {}


def check_coverage(components, footprints, pad_map):
    """Fail loudly, before touching the board, if any device type in this .circuit file has no footprint/pad
    mapping yet - never guess a footprint or silently drop a pad for an unrecognized device type."""
    missing = sorted({dev for dev in components.values() if dev not in footprints or dev not in pad_map})
    if missing:
        sys.exit(
            "No footprint/pad mapping for device type(s): " + ", ".join(missing) + "\n"
            "Add entries to FOOTPRINTS and PAD_MAP in this script (verify against a real datasheet/footprint "
            "file, same as the existing entries) before running this .circuit file through the translator."
        )


def _match_paren(text, open_pos):
    """Index of the ')' matching the '(' at open_pos, respecting quoted strings (with backslash escapes)."""
    depth = 0
    in_str = False
    i = open_pos
    while i < len(text):
        ch = text[i]
        if in_str:
            if ch == "\\":
                i += 2
                continue
            if ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("unbalanced parens")


def canonicalize_pcb_text(text):
    """board.Save() writes footprints in an order that depends on pcbnew's own internal container (observed:
    genuinely different every run, unrelated to our own board.Add() insertion order or the .circuit file's
    content) - not something this script controls, so it's fixed here instead: pull out every top-level
    (footprint ...) block by its REAL "Reference" property (not by position - position isn't trustworthy, see
    above), sort those blocks by reference, and splice them back in that fixed order. Every uuid inside a block
    is also rewritten deterministically, keyed by (that block's real reference, its ordinal position within the
    block) - stable because the same library footprint always generates the same internal item structure/order
    every time, confirmed empirically. Net effect: same .circuit input -> byte-identical .kicad_pcb output,
    regardless of whatever internal order/randomness pcbnew's own Save() used this particular run."""
    starts = [m.start() for m in re.finditer(r'^\t\(footprint "', text, re.MULTILINE)]
    if not starts:
        return text
    ends = [_match_paren(text, s + 1) for s in starts]  # s points at the tab; s+1 is the opening '('
    spans = []
    for s, e in zip(starts, ends):
        e2 = e + 1
        if e2 < len(text) and text[e2] == "\n":
            e2 += 1
        spans.append((s, e2))
    preamble = text[:spans[0][0]]
    tail = text[spans[-1][1]:]
    tagged = []
    for s, e in spans:
        block = text[s:e]
        m = re.search(r'\(property "Reference" "([^"]*)"', block)
        ref = m.group(1) if m else ""
        tagged.append((ref, block))
    tagged.sort(key=lambda t: t[0])

    out_blocks = []
    for ref, block in tagged:
        counter = 0

        def repl(m, ref=ref):
            nonlocal counter
            key = f"{ref}#{counter}"
            counter += 1
            return f'(uuid "{uuid.uuid5(UUID_NAMESPACE, key)}")'

        block = re.sub(r'\(uuid "[0-9a-fA-F-]+"\)', repl, block)
        out_blocks.append(block)
    return preamble + "".join(out_blocks) + tail


def build_board(components, nets, nc, overrides, libdir, out_path, ncols, values=None):
    values = values or {}
    term_net = {t: n["name"] for n in nets for t in n["terminals"]}

    board = pcbnew.NewBoard(out_path)
    netinfo = {}
    for n in nets:
        ni = pcbnew.NETINFO_ITEM(board, n["name"])
        board.Add(ni)
        netinfo[n["name"]] = ni

    unmapped_pads = []
    step_mm = 25
    for i, ref in enumerate(sorted(components)):
        dev = components[ref]
        libsub, fpname = overrides.get(ref, FOOTPRINTS[dev])
        fp = pcbnew.FootprintLoad(f"{libdir}/{libsub}", fpname)
        if fp is None:
            sys.exit(f"Could not load footprint {libsub}/{fpname} for {ref} ({dev}) from {libdir}")
        fp.SetReference(ref)
        fp.SetValue(values.get(ref, dev))  # real component value/description when the .circuit file has one
        value_field = fp.Value()
        value_field.SetLayer(pcbnew.F_SilkS)  # library default is F.Fab (fab notes, not printed) - move onto
        value_field.SetVisible(True)          # the real silkscreen so it's actually printed on the board
        col, row = i % ncols, i // ncols
        fp.SetPosition(pcbnew.VECTOR2I(pcbnew.FromMM(20 + col * step_mm), pcbnew.FromMM(20 + row * step_mm)))
        pad_map = PAD_MAP[dev]
        inv_pad_map = {v: k for k, v in pad_map.items()}
        for pad in fp.Pads():
            padnum = pad.GetNumber()
            circuit_pin = inv_pad_map.get(padnum)
            if circuit_pin is None:
                unmapped_pads.append((ref, padnum))
                continue
            term = f"{ref}.{circuit_pin}"
            if term in nc:
                continue  # leave NC pads on no net, on purpose
            net_name = term_net.get(term)
            if net_name is None:
                unmapped_pads.append((ref, padnum, "no net found for", term))
                continue
            pad.SetNet(netinfo[net_name])
        board.Add(fp)

    board.Save(out_path)
    text = canonicalize_pcb_text(open(out_path).read())
    open(out_path, "w").write(text)
    return unmapped_pads


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("circuit", help="path to the .circuit netlist file")
    ap.add_argument("output", help="path to write the .kicad_pcb file to")
    ap.add_argument("--overrides", help="optional JSON file: {\"REF\": [\"lib.pretty\", \"footprint_name\"]} "
                                         "for board-specific footprint choices (e.g. one capacitor that should "
                                         "be a real electrolytic instead of the device type's default)")
    ap.add_argument("--libdir", default=LIBDIR_DEFAULT, help="KiCad footprint library root (.pretty dirs live here)")
    ap.add_argument("--cols", type=int, default=6, help="placeholder placement grid width (components per row)")
    args = ap.parse_args()

    components, nets, nc = parse_circuit(args.circuit)
    if not components:
        sys.exit(f"No `component` lines found in {args.circuit} - is this really a .circuit file?")
    check_coverage(components, FOOTPRINTS, PAD_MAP)

    overrides = {}
    if args.overrides:
        raw = json.loads(open(args.overrides).read())
        overrides = {ref: tuple(v) for ref, v in raw.items()}

    values = parse_values(args.circuit)
    unmapped_pads = build_board(components, nets, nc, overrides, args.libdir, args.output, args.cols, values)

    print(f"saved {args.output}")
    print(f"components: {len(components)}  nets: {len(nets)}  NC pins left unconnected: {len(nc)}")
    if unmapped_pads:
        print("UNMAPPED PADS (real problem, needs fixing):", unmapped_pads)
    else:
        print("every pad on every footprint was matched to a net or intentionally left NC - none skipped")


if __name__ == "__main__":
    main()
