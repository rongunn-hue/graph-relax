#!/usr/bin/env python3
"""Rendering-only views of the completed common-metric direct experiment."""
import argparse, html, json, math, xml.etree.ElementTree as ET
from pathlib import Path

SVGNS="{http://www.w3.org/2000/svg}"
STAGES=("start","polygon","nodes","final")

def saved_coordinates(svg, vertices):
    """Recover authoritative graph-node circles from an already-saved stage SVG."""
    names={x["name"]:x["id"] for x in vertices.values()};root=ET.parse(svg).getroot();children=list(root);out={}
    for i,e in enumerate(children[:-1]):
        if e.tag==SVGNS+"circle" and e.get("fill")=="#1565c0":
            t=children[i+1]
            if t.tag==SVGNS+"text" and (t.text or "") in names:
                out[names[t.text]]=(float(e.get("cx")),float(e.get("cy")))
    if len(out)!=len(vertices):raise AssertionError(f"{svg}: recovered {len(out)} graph coordinates")
    return out

def bounds(co):
    xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];return min(xs),min(ys),max(xs),max(ys)

def graph_elements(vertices,edges,co,label_size=13):
    z=[]
    for a,b in edges:z.append(f'<line x1="{co[a][0]:.9f}" y1="{co[a][1]:.9f}" x2="{co[b][0]:.9f}" y2="{co[b][1]:.9f}" stroke="#777" stroke-width="1.35"/>')
    # Stable radial label offset away from the drawing centre reduces collisions with incident stars.
    cx=sum(p[0] for p in co.values())/len(co);cy=sum(p[1] for p in co.values())/len(co)
    for v in sorted(vertices):
        x,y=co[v];dx,dy=x-cx,y-cy;L=math.hypot(dx,dy) or 1;ox,oy=13*dx/L,13*dy/L
        anchor="end" if ox<0 else "start"
        z.append(f'<circle cx="{x:.9f}" cy="{y:.9f}" r="4.3" fill="#1565c0" stroke="white" stroke-width="1"/>')
        z.append(f'<text x="{x+ox:.9f}" y="{y+oy:.9f}" font-size="{label_size}" font-weight="600" text-anchor="{anchor}" dominant-baseline="middle" fill="#111">{html.escape(vertices[v]["name"])}</text>')
    return z

def clean_svg(vertices,edges,co,title,path):
    x0,y0,x1,y1=bounds(co);pad=65;w=x1-x0+2*pad;h=y1-y0+2*pad
    z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0-pad} {y0-pad} {w} {h}" width="1500" height="1100">',f'<rect x="{x0-pad}" y="{y0-pad}" width="{w}" height="{h}" fill="white"/>',f'<text x="{x0-pad+12}" y="{y0-pad+25}" font-size="18" font-weight="bold">{html.escape(title)}</text>']
    z+=graph_elements(vertices,edges,co);z.append('</svg>');path.write_text('\n'.join(z)+'\n')

def event_rows(report):
    prefixes=(("polygon","P"),("node_node","NN"),("node_edge","NE"),("unrelated_edge_edge","EE"),("shared_endpoint_overlay","SE"));rows=[]
    for key,prefix in prefixes:
        for i,e in enumerate(report["initial_map"][key],1):
            x=dict(e);x["display_id"]=f"{prefix}{i:02d}";rows.append(x)
    return rows

def detail(e):
    a=e["affected"]
    if e["type"]=="NODE_EDGE": geom=f'{a["vertex"]} vs {a["edge"][0]}–{a["edge"][1]}'
    elif e["type"]=="NODE_NODE":geom=f'{a[0]} ↔ {a[1]}'
    elif e["type"]=="EDGE_EDGE":geom=f'{a[0][0]}–{a[0][1]} vs {a[1][0]}–{a[1][1]}'
    elif e["type"]=="SHARED_OVERLAY":geom=f'{a[0][0]}–{a[0][1]} / {a[1][0]}–{a[1][1]}'
    else:geom=str(a)
    return f'{e["display_id"]} {e["type"]}  {geom}  d={e["d"]:.3f}  D={e["D"]:.3f}  d/D={e["ratio"]:.4f}'

def diagnostic_svg(vertices,edges,co,report,path):
    rows=event_rows(report);x0,y0,x1,y1=bounds(co);pad=65;gw=x1-x0+2*pad;gh=max(y1-y0+2*pad,950);cols=3;per=math.ceil(len(rows)/cols);colw=690;tablew=cols*colw;X=x0-pad;Y=y0-pad
    colors={"POLYGON":"#ef6c00","NODE_NODE":"#7b1fa2","NODE_EDGE":"#00838f","EDGE_EDGE":"#d32f2f","SHARED_OVERLAY":"#2e7d32"}
    z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X} {Y} {gw+tablew} {gh}" width="2600" height="1200">',f'<rect x="{X}" y="{Y}" width="{gw+tablew}" height="{gh}" fill="white"/>',f'<text x="{X+12}" y="{Y+24}" font-size="18" font-weight="bold">COMMON-METRIC DIRECT — START DIAGNOSTIC</text>']
    z+=graph_elements(vertices,edges,co,10)
    for e in rows:
        c=colors[e["type"]];eid=e["display_id"]
        if e["type"]=="NODE_NODE":
            a,b=e["affected"];p,q=co[a],co[b];z.append(f'<line x1="{p[0]}" y1="{p[1]}" x2="{q[0]}" y2="{q[1]}" stroke="{c}" stroke-width="1" stroke-dasharray="3 3"/>')
        elif e["type"]=="NODE_EDGE":
            p=co[e["affected"]["vertex"]];q=e["closest_point"];z.append(f'<line x1="{p[0]}" y1="{p[1]}" x2="{q[0]}" y2="{q[1]}" stroke="{c}" stroke-width="1"/>');z.append(f'<circle cx="{q[0]}" cy="{q[1]}" r="2.5" fill="{c}"/>')
        elif e["type"]=="EDGE_EDGE":
            p,q=e["closest_points"];z.append(f'<line x1="{p[0]}" y1="{p[1]}" x2="{q[0]}" y2="{q[1]}" stroke="{c}" stroke-width="1"/>');z.extend(f'<circle cx="{r[0]}" cy="{r[1]}" r="2.5" fill="{c}"/>' for r in (p,q))
        elif e["type"]=="SHARED_OVERLAY":
            v=e["shared_vertex"]
            for a in e["nonshared_endpoints"]:
                p,q=co[v],co[a];z.append(f'<line x1="{p[0]}" y1="{p[1]}" x2="{p[0]+.28*(q[0]-p[0])}" y2="{p[1]+.28*(q[1]-p[1])}" stroke="{c}" stroke-width="3" opacity=".75"/>')
        elif e["type"]=="POLYGON":
            pts=' '.join(f'{co[v][0]},{co[v][1]}' for v in e["boundary"]);z.append(f'<polygon points="{pts}" fill="none" stroke="{c}" stroke-width="1" stroke-dasharray="4 3" opacity=".7"/>')
        px,py=e["location"];z.append(f'<text x="{px+4}" y="{py-4}" font-size="7" fill="{c}" stroke="white" stroke-width="1" paint-order="stroke">{eid}</text>')
    tx=X+gw+12;z.append(f'<text x="{tx}" y="{Y+22}" font-size="12" font-weight="bold">EVENT TABLE ({len(rows)} annotations)</text>');z.append(f'<text x="{tx}" y="{Y+39}" font-size="10">P polygon · NN node-node · NE node-edge · EE unrelated edges · SE shared endpoint</text>')
    for n,e in enumerate(rows):
        col=n//per;row=n%per;xx=tx+col*colw;yy=Y+58+row*15
        z.append(f'<text x="{xx}" y="{yy}" font-size="8.4" fill="{colors[e["type"]]}">{html.escape(detail(e))}</text>')
    # Explicit audit callouts, kept beside the graph rather than across it.
    z.append(f'<text x="{tx}" y="{Y+gh-38}" font-size="10" font-weight="bold">TPO: own TPO–EOUT excluded; principal unrelated edge event is NE{next(i for i,e in enumerate(report["initial_map"]["node_edge"],1) if e["affected"]["vertex"]=="component:TPO"):02d} vs J_R4–EOUT.</text>')
    z.append(f'<text x="{tx}" y="{Y+gh-20}" font-size="10" font-weight="bold">TPB: NODE_B–TPB / NODE_B–R3 separation is measured away from legitimate endpoint NODE_B.</text>')
    z.append(f'<!-- EVENT_ANNOTATION_COUNT={len(rows)} --></svg>');path.write_text('\n'.join(z)+'\n');return len(rows)

def run(outdir):
    report=json.loads((outdir/"iamp_common_metric_direct_report.json").read_text());gd=json.loads((outdir/"iamp_graph.json").read_text());vertices={x["id"]:x for x in gd["vertices"]};edges=[(x["source"],x["target"]) for x in gd["edges"]]
    maps={s:saved_coordinates(outdir/f"iamp_common_metric_direct_{s}.svg",vertices) for s in STAGES}
    direct=json.loads((outdir/"iamp_nearness_rotation_direct_report.json").read_text());expected_start={v:tuple(p) for v,p in direct["final_coordinates"].items()};expected_final={v:tuple(p) for v,p in report["final_coordinates"].items()}
    assert maps["start"]==expected_start and maps["final"]==expected_final
    files=[]
    for s in STAGES:
        p=outdir/f"iamp_common_metric_direct_{s}_clean.svg";clean_svg(vertices,edges,maps[s],f"IAMP COMMON-METRIC DIRECT — {s.upper()}",p);files.append(p)
    p=outdir/"iamp_common_metric_direct_start_diagnostic.svg";count=diagnostic_svg(vertices,edges,maps["start"],report,p);assert count==sum(len(v) for v in report["initial_map"].values());files.append(p)
    # Rendering-only invariant: recovered source coordinates remain byte-for-number identical.
    assert maps["start"]==expected_start and maps["final"]==expected_final
    return files,maps

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,default=Path("output/graph_first"));a=p.parse_args();files,_=run(a.output);print('\n'.join(map(str,files)))
if __name__=="__main__":main()
