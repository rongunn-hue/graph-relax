#!/usr/bin/env python3
"""Third and final stage: standalone, host-side script (no flatpak/pcbnew needed - kicad-cli alone does
everything here). Takes a placed+routed .kicad_pcb and produces the actual files a fab house needs.

Usage:
    python3 export_fab.py <board.kicad_pcb> <output_dir>

Refuses to produce a fab package from a board that isn't DRC-clean - a real, freshly-run check each time, not
an assumption carried over from an earlier stage.
"""
import argparse
import json
import subprocess
import sys
import zipfile
from pathlib import Path

GERBER_LAYERS = "F.Cu,B.Cu,F.Mask,B.Mask,F.Silkscreen,B.Silkscreen,Edge.Cuts"


def run(cmd):
    print("+", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        sys.exit(f"command failed: {' '.join(cmd)}")
    if r.stdout.strip():
        print(r.stdout.strip())
    return r.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("board", help="placed+routed .kicad_pcb file")
    ap.add_argument("outdir", help="directory to write everything into")
    args = ap.parse_args()

    board = args.board
    outdir = Path(args.outdir)
    fabdir = outdir / "fab"
    fabdir.mkdir(parents=True, exist_ok=True)

    print("=== DRC (real, freshly run - not assumed clean) ===")
    drc_path = outdir / "drc.json"
    run(["kicad-cli", "pcb", "drc", "--format", "json", "--output", str(drc_path), board])
    drc = json.loads(drc_path.read_text())
    nviol, nunc = len(drc["violations"]), len(drc["unconnected_items"])
    print(f"violations: {nviol}  unconnected_items: {nunc}")
    if nviol or nunc:
        sys.exit("DRC is not clean - refusing to produce a fab package from a board with real problems")

    print("\n=== board statistics (explains hole counts - vias are real drill holes too, not just pins) ===")
    stats_path = outdir / "stats.json"
    run(["kicad-cli", "pcb", "export", "stats", "--format", "json", "--output", str(stats_path), board])
    stats = json.loads(stats_path.read_text())
    print(json.dumps(stats, indent=2)[:1500])

    print("\n=== gerbers ===")
    run(["kicad-cli", "pcb", "export", "gerbers", "--output", f"{fabdir}/", "--layers", GERBER_LAYERS, board])

    print("\n=== drill ===")
    run(["kicad-cli", "pcb", "export", "drill", "--output", f"{fabdir}/", "--format", "excellon",
         "--drill-origin", "absolute", board])

    print("\n=== flat top/bottom plots (2D overlay, for tracing/inspection - copper+silkscreen+mask+outline) ===")
    # PDF export, not SVG: kicad-cli's SVG->PNG path (via librsvg/ImageMagick) was found to genuinely mis-render
    # some glyphs (confirmed on real evidence: the digit "4" came out with a malformed extra stroke at 900dpi,
    # in a direct A/B render against KiCad's own PDF plot of the identical board, which was clean) - this is a
    # real bug in that specific rasterization path, not the font and not the board data. PDF (via kicad-cli's
    # own plotter, rasterized by poppler/pdftoppm) does not show this defect. --bg-color is a real kicad-cli
    # flag, so the green background no longer needs an ImageMagick post-process step either.
    BOARD_GREEN = "#0f4d2a"  # real soldermask-green background, not white - yellow silkscreen/copper read poorly on white
    top_pdf = outdir / "render_top.pdf"
    bottom_pdf = outdir / "render_bottom.pdf"
    top_png = outdir / "render_top.png"
    bottom_png = outdir / "render_bottom.png"
    run(["kicad-cli", "pcb", "export", "pdf", "--mode-single", "--output", str(top_pdf),
         "--layers", "F.Cu,F.Silkscreen,F.Mask,Edge.Cuts", "--bg-color", BOARD_GREEN, board])
    run(["kicad-cli", "pcb", "export", "pdf", "--mode-single", "--output", str(bottom_pdf),
         "--layers", "B.Cu,B.Silkscreen,B.Mask,Edge.Cuts", "--mirror", "--bg-color", BOARD_GREEN, board])
    run(["pdftoppm", "-png", "-r", "600", "-singlefile", str(top_pdf), str(outdir / "render_top")])
    run(["pdftoppm", "-png", "-r", "600", "-singlefile", str(bottom_pdf), str(outdir / "render_bottom")])
    # unlike the svg export, pdf export has no --fit-page-to-board - it always plots the full fixed page size
    # (e.g. A4), so a small board ends up as a tiny cluster in one corner of a mostly-blank page. Trim to real
    # content instead of leaving that page margin in the delivered image (safe here: no border/title block was
    # requested, so the only non-background content on the page is the board itself).
    run(["convert", str(top_png), "-trim", "+repage", "-bordercolor", BOARD_GREEN, "-border", "40", str(top_png)])
    run(["convert", str(bottom_png), "-trim", "+repage", "-bordercolor", BOARD_GREEN, "-border", "40", str(bottom_png)])

    print("\n=== 3D render (real raytraced view, angled isometric) ===")
    render_3d_png = outdir / "render_3d.png"
    run(["kicad-cli", "pcb", "render", "--rotate", "-45,0,45", "--perspective", "--width", "1600",
         "--height", "1200", "--quality", "high", "--background", "opaque", "--output", str(render_3d_png),
         board])

    print("\n=== zipping fab package ===")
    zip_path = outdir / "fab_package.zip"
    fab_files = sorted(fabdir.iterdir())
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in fab_files:
            zf.write(f, f.name)
    print(f"{zip_path} ({len(fab_files)} files: {', '.join(f.name for f in fab_files)})")

    print("\nDONE")
    print(f"fab package: {zip_path}")
    print(f"flat plots: {top_png}, {bottom_png}")
    print(f"3D render: {render_3d_png}")
    print(f"DRC report: {drc_path}")
    print(f"stats report: {stats_path}")


if __name__ == "__main__":
    main()
