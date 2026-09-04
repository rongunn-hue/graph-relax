#!/usr/bin/env python3
"""One frozen-cell Voronoi-centroid movement pass on the IAMP graph."""
import argparse, json, math, statistics
from pathlib import Path
import placement_geometry as geometry
import graph_first_experiment_global_node_density as density

def centroid(poly):
    twice=sum(p[0]*q[1]-q[0]*p[1] for p,q in zip(poly,poly[1:]+poly[:1]))
    if abs(twice)<1e-15: raise AssertionError("degenerate Voronoi cell")
    x=sum((p[0]+q[0])*(p[0]*q[1]-q[0]*p[1]) for p,q in zip(poly,poly[1:]+poly[:1]))/(3*twice)
    y=sum((p[1]+q[1])*(p[0]*q[1]-q[0]*p[1]) for p,q in zip(poly,poly[1:]+poly[:1]))/(3*twice)
    return (x,y)

def cells_in_fixed_box(co,box):
    xmin,xmax,ymin,ymax=box; rectangle=[(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)]; result={}
    for v in sorted(co):
        x,y=co[v]; poly=list(rectangle)
        for w in sorted(co):
            if w==v: continue
            u,z=co[w]
            poly=density.clip_halfplane(poly,2*(u-x),2*(z-y),u*u+z*z-x*x-y*y)
            if not poly: raise AssertionError(f"empty fixed-domain Voronoi cell {v}")
        result[v]=poly
    return result

def diagnostics(co,edges):
    d=geometry.graph_geometry_diagnostics(co,edges)
    return {"proper_crossings":d["proper_unrelated_edge_crossing_count"],
            "proper_crossing_list":d["proper_unrelated_edge_crossings"],
            "coincident_vertices":d["coincident_vertex_pair_count"],
            "coincident_vertex_pairs":d["coincident_vertex_pairs"],
            "vertex_on_unrelated_edge_events":d["vertex_on_unrelated_edge_interior_count"],
            "vertex_on_unrelated_edge_list":d["vertices_on_unrelated_edge_interiors"]}

def records(cellmap,target):
    areas={v:density.polygon_area(p) for v,p in cellmap.items()}
    order=sorted(areas,key=lambda v:(areas[v],v)); out=[]
    for rank,v in enumerate(order,1):
        a=areas[v]; out.append({"rank":rank,"vertex_id":v,"voronoi_area":a,"area_over_target":a/target,
          "crowding_ratio":target/a,"percent_difference_from_target":100*(a/target-1)})
    return out

def stats(rows):
    a=[r["voronoi_area"] for r in rows]
    return {"minimum_area":min(a),"maximum_area":max(a),"median_area":statistics.median(a),
      "mean_area":statistics.mean(a),"population_standard_deviation":statistics.pstdev(a),"max_min_area_ratio":max(a)/min(a)}

def panel(vertices,edges,co,cells,ratios,box,xoff,title,moves=None):
    xmin,xmax,ymin,ymax=box; z=[f'<g transform="translate({xoff:.17g} 0)">',f'<rect x="{xmin:.17g}" y="{ymin:.17g}" width="{xmax-xmin:.17g}" height="{ymax-ymin:.17g}" fill="white" stroke="#111" stroke-width="2"/>']
    for i,v in enumerate(sorted(cells)):
        z.append('<polygon points="'+' '.join(f'{x:.17g},{y:.17g}' for x,y in cells[v])+f'" fill="hsl({(i*137)%360} 55% 70%)" fill-opacity=".25" stroke="#78909c" stroke-width=".8"/>')
    if moves:
        for m in moves:
            a=m["starting_coordinate_at_attempt"]; b=m["frozen_cell_centroid"]
            if m["decision"]=="ACCEPTED":
                z.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#2e7d32" stroke-width="2" stroke-dasharray="5 4"/>')
            else:
                z.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#c62828" stroke-width="1.5" stroke-dasharray="2 4"/>')
                z.append(f'<path d="M {b[0]-5:.17g} {b[1]-5:.17g} L {b[0]+5:.17g} {b[1]+5:.17g} M {b[0]-5:.17g} {b[1]+5:.17g} L {b[0]+5:.17g} {b[1]-5:.17g}" stroke="#c62828" stroke-width="2"/>')
    for a,b in edges:z.append(f'<line x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#444" stroke-width="1.5"/>')
    for v in sorted(vertices):
        x,y=co[v]; fill="#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
        z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{fill}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]} CR={ratios[v]:.3f}</text>')
    z.append(f'<text x="{xmin+12:.17g}" y="{ymin+24:.17g}" font-family="sans-serif" font-size="18" font-weight="bold">{title}</text></g>'); return z

def render(vertices,edges,initial,final,initial_cells,final_cells,before,after,box,attempts,path):
    xmin,xmax,ymin,ymax=box; w=xmax-xmin; h=ymax-ymin; gap=40
    rb={r["vertex_id"]:r["crowding_ratio"] for r in before}; ra={r["vertex_id"]:r["crowding_ratio"] for r in after}
    z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{xmin:.17g} {ymin:.17g} {2*w+gap:.17g} {h:.17g}" width="1800" height="900">']
    z+=panel(vertices,edges,initial,initial_cells,rb,box,0,"BEFORE — frozen cells")
    z+=panel(vertices,edges,final,final_cells,ra,box,w+gap,"AFTER — recomputed cells",attempts)
    z.append('</svg>'); path.write_text('\n'.join(z)+'\n')

def run(graph_path,density_report_path,outdir):
    g=json.loads(graph_path.read_text()); base=json.loads(density_report_path.read_text())
    vertices={v["id"]:v for v in g["vertices"]}; edges=[(e["source"],e["target"]) for e in g["edges"]]
    initial={v:tuple(p) for v,p in json.loads(Path("output/graph_first/iamp_sequential_local_face_pulls_report.json").read_text())["stages"][1]["coordinates"].items()}
    # Prove the reconstructed cells and ranking are exactly the prior measurement.
    box,initial_cells=density.cells(initial); target=base["target_area"]; before=records(initial_cells,target)
    if box!=(base["domain"]["xmin"],base["domain"]["xmax"],base["domain"]["ymin"],base["domain"]["ymax"]): raise AssertionError("domain mismatch")
    if [r["vertex_id"] for r in before]!=[r["vertex_id"] for r in base["ranked_vertices"]]: raise AssertionError("ranking mismatch")
    co=dict(initial); attempts=[]
    for row in before:
        v=row["vertex_id"]; start=co[v]; c=centroid(initial_cells[v]); co[v]=c
        d=diagnostics(co,edges); reasons=[]
        if d["proper_crossings"]: reasons.append({"type":"proper_crossings","details":d["proper_crossing_list"]})
        if d["coincident_vertices"]: reasons.append({"type":"coincident_vertices","details":d["coincident_vertex_pairs"]})
        if d["vertex_on_unrelated_edge_events"]: reasons.append({"type":"vertex_on_unrelated_edge","details":d["vertex_on_unrelated_edge_list"]})
        decision="REJECTED" if reasons else "ACCEPTED"
        if reasons: co[v]=start
        attempts.append({"attempt_order":row["rank"],"vertex_id":v,"original_frozen_voronoi_area":row["voronoi_area"],
          "original_crowding_ratio":row["crowding_ratio"],"starting_coordinate_at_attempt":list(start),"frozen_cell_centroid":list(c),
          "proposed_displacement_distance":math.dist(start,c),"decision":decision,"rejection_reasons":reasons,
          "complete_edge_regeneration_count":len(edges)})
    final_diag=diagnostics(co,edges)
    if any(final_diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")): raise AssertionError("invalid final graph")
    final_cells=cells_in_fixed_box(co,box); after=records(final_cells,target)
    if abs(sum(r["voronoi_area"] for r in after)-base["domain"]["area"])>1e-7: raise AssertionError("final partition mismatch")
    bb={r["vertex_id"]:r for r in before}; aa={r["vertex_id"]:r for r in after}
    comparison=[{"vertex_id":v,"original_rank":bb[v]["rank"],"final_rank":aa[v]["rank"],"original_area":bb[v]["voronoi_area"],"final_area":aa[v]["voronoi_area"],"original_crowding_ratio":bb[v]["crowding_ratio"],"final_crowding_ratio":aa[v]["crowding_ratio"]} for v in sorted(vertices)]
    accepted=[a for a in attempts if a["decision"]=="ACCEPTED"]
    report={"schema":"graph-relax.voronoi-centroid-move.v1","source_density_report":str(density_report_path),"domain":base["domain"],"target_area":target,
      "vertex_count":len(vertices),"electrical_incidence_count":len(edges),"initial_proper_crossings":0,"frozen_initial_partition":True,
      "attempts":attempts,"attempt_count":len(attempts),"accepted_moves":len(accepted),"rejected_moves":len(attempts)-len(accepted),
      "total_vertex_displacement":sum(a["proposed_displacement_distance"] for a in accepted),"maximum_individual_displacement":max(a["proposed_displacement_distance"] for a in accepted),
      "final_coordinates":{v:list(co[v]) for v in sorted(co)},"final_validation":final_diag,"before_density":before,"after_density":after,
      "before_statistics":stats(before),"after_statistics":stats(after),"before_after_by_vertex":comparison,
      "final_partition_area_sum":sum(r["voronoi_area"] for r in after),"final_partition_area_error":sum(r["voronoi_area"] for r in after)-base["domain"]["area"]}
    outdir.mkdir(parents=True,exist_ok=True); render(vertices,edges,initial,co,initial_cells,final_cells,before,after,box,attempts,outdir/"iamp_voronoi_centroid_move.svg")
    density.render(vertices,edges,co,after,final_cells,box,outdir/"iamp_voronoi_centroid_after.svg")
    (outdir/"iamp_voronoi_centroid_move_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n'); return report

def main():
    p=argparse.ArgumentParser(); b=Path("output/graph_first"); p.add_argument("--graph",type=Path,default=b/"iamp_graph.json"); p.add_argument("--density-report",type=Path,default=b/"iamp_global_node_density_report.json"); p.add_argument("--output",type=Path,default=b); a=p.parse_args(); r=run(a.graph,a.density_report,a.output)
    print("accepted",r["accepted_moves"],"rejected",r["rejected_moves"],"total displacement",r["total_vertex_displacement"])
    for x in r["attempts"]: print(x["attempt_order"],x["vertex_id"],x["decision"],x["proposed_displacement_distance"],x["rejection_reasons"])
if __name__=="__main__": main()
