#!/usr/bin/env python3
"""Measure and, only when finite on both sides, apply half-space corridor opening."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_corridor_separation as corridor
import placement_geometry as geometry


def ray_segment_hit(origin, direction, first, second):
    """Return forward-ray distance and segment fraction, or None."""
    ex,ey=second[0]-first[0],second[1]-first[1]
    determinant=direction[0]*ey-direction[1]*ex
    scale=max(1.,math.dist(first,second))
    if abs(determinant)<=geometry.GEOMETRY_TOLERANCE*scale:return None
    ax,ay=first[0]-origin[0],first[1]-origin[1]
    distance=(ax*ey-ay*ex)/determinant
    fraction=(ax*direction[1]-ay*direction[0])/determinant
    if distance<=geometry.GEOMETRY_TOLERANCE or fraction < -geometry.GEOMETRY_TOLERANCE or fraction > 1+geometry.GEOMETRY_TOLERANCE:return None
    return distance,max(0.,min(1.,fraction))


def available_distance(vertices,direction,coordinates,edges,boundary_edges):
    """First unrelated graph boundary on each participating vertex's normal ray."""
    records=[]
    boundary={corridor.edge_key(*e) for e in boundary_edges}
    for vertex in vertices:
        hits=[]
        for edge in edges:
            key=corridor.edge_key(*edge)
            if key in boundary or vertex in edge:continue
            result=ray_segment_hit(coordinates[vertex],direction,coordinates[edge[0]],coordinates[edge[1]])
            if result:
                distance,fraction=result
                hits.append({"distance":distance,"limiting_edge":list(edge),"segment_parameter":fraction,
                             "intersection":[coordinates[vertex][0]+distance*direction[0],coordinates[vertex][1]+distance*direction[1]]})
        hits.sort(key=lambda x:(x["distance"],x["limiting_edge"]))
        records.append({"vertex":vertex,"ray_origin":list(coordinates[vertex]),"direction":list(direction),
                        "first_hit":hits[0] if hits else None,"status":"FINITE" if hits else "UNBOUNDED / NO LIMITING BOUNDARY"})
    finite=[r for r in records if r["first_hit"]]
    if len(finite)!=len(records):return None,None,records
    limiting=min(finite,key=lambda r:(r["first_hit"]["distance"],r["vertex"]))
    return limiting["first_hit"]["distance"],limiting,records


def render(base,fixed,report,output):
    vertices=base["vertices"];edges=base["edges"];before=base["before"];fixed_after=base["fixed_after"]
    xs=[p[0] for p in before.values()];ys=[p[1] for p in before.values()];margin=90
    minx,maxx,miny,maxy=min(xs),max(xs),min(ys),max(ys);pw=maxx-minx+2*margin;h=maxy-miny+2*margin
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx-margin:.17g} {miny-margin:.17g} {3*pw:.17g} {h:.17g}" width="2100" height="800">',f'<rect x="{minx-margin:.17g}" y="{miny-margin:.17g}" width="{3*pw:.17g}" height="{h:.17g}" fill="white"/>']
    panels=(("ORIGINAL",before),("FIXED delta=0.5",fixed_after),("HALF AVAILABLE — UNDEFINED; UNCHANGED",before))
    paths=(report["corridor"]["rmin1_side_boundary"]["ordered_vertices"],report["corridor"]["j_out_side_boundary"]["ordered_vertices"])
    start,end=(tuple(x) for x in report["centerline"]["endpoints"]);normal=tuple(report["centerline"]["unit_normal_toward_rmin1_side"])
    for index,(title,co) in enumerate(panels):
        shift=index*pw;parts.append(f'<g id="panel-{index}" transform="translate({shift:.17g} 0)">')
        for a,b in edges:parts.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        for path,color in zip(paths,("#00897b","#8e24aa")):parts.append('<polyline points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in path)+f'" fill="none" stroke="{color}" stroke-width="5" opacity=".7"/>')
        parts.append(f'<line x1="{start[0]:.17g}" y1="{start[1]:.17g}" x2="{end[0]:.17g}" y2="{end[1]:.17g}" stroke="#1565c0" stroke-width="2" stroke-dasharray="8 5"/>')
        for v in sorted(vertices):
            x,y=co[v];color="#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            parts.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-size="10">{vertices[v]["name"]}</text>')
        if index==2:
            for side,color in (("side_a_available_space","#00897b"),("side_b_available_space","#8e24aa")):
                for ray in report[side]["per_vertex_rays"]:
                    p=ray["ray_origin"];hit=ray["first_hit"]
                    if hit:q=hit["intersection"]
                    else:q=[p[0]+130*ray["direction"][0],p[1]+130*ray["direction"][1]]
                    parts.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{color}" stroke-width="3" stroke-dasharray="7 4"/>')
                    if hit:
                        e=hit["limiting_edge"];a,b=before[e[0]],before[e[1]]
                        parts.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#f9a825" stroke-width="7" opacity=".55"/>')
                    else:parts.append(f'<text x="{q[0]+5:.17g}" y="{q[1]:.17g}" font-size="11" fill="{color}">no finite boundary</text>')
        parts.append(f'<text x="{minx-margin+14:.17g}" y="{miny-margin+25:.17g}" font-size="17" font-weight="bold">{title}</text></g>')
    parts.append('</svg>');output.write_text('\n'.join(parts)+'\n')


def run(graph_path,fixed_report_path,output_dir):
    fixed=json.loads(fixed_report_path.read_text())
    g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]]
    before={v:tuple(p) for v,p in fixed["coordinates_before"].items()};fixed_after={v:tuple(p) for v,p in fixed["coordinates_after"].items()}
    if (len(vertices),len(edges),fixed["before_metrics"]["proper_crossings"],fixed["before_metrics"]["clearance_violations"])!=(34,44,0,6):raise AssertionError("baseline differs from established corridor experiment")
    n=tuple(fixed["centerline"]["unit_normal_toward_rmin1_side"]);minus=(-n[0],-n[1])
    ba=[tuple(x) for x in fixed["corridor"]["rmin1_side_boundary"]["ordered_edges"]];bb=[tuple(x) for x in fixed["corridor"]["j_out_side_boundary"]["ordered_edges"]]
    da,limit_a,rays_a=available_distance(fixed["rmin1_side_movable_vertices"],n,before,edges,ba)
    db,limit_b,rays_b=available_distance(fixed["j_out_side_movable_vertices"],minus,before,edges,bb)
    applicable=da is not None and db is not None
    # The rule defines a simultaneous two-sided operation. An unbounded side
    # has no finite half-distance, so no partial or substituted move is made.
    final=dict(before);delta_a=da/2 if da is not None else None;delta_b=db/2 if db is not None else None
    if applicable:
        for v in fixed["rmin1_side_movable_vertices"]:final[v]=(before[v][0]+delta_a*n[0],before[v][1]+delta_a*n[1])
        for v in fixed["j_out_side_movable_vertices"]:final[v]=(before[v][0]+delta_b*minus[0],before[v][1]+delta_b*minus[1])
    diagnostics,analysis,events,clusters,final_metrics=corridor.metrics(vertices,edges,final)
    segments=geometry.rebuild_straight_edge_segments(final,edges);reconnect=geometry.validate_rebuilt_edge_segments(final,edges,segments)
    report={"schema":"graph-relax.corridor-half-available.v1","status":"APPLIED" if applicable else "NOT APPLICABLE — AT LEAST ONE OUTWARD SIDE HAS NO FINITE LIMITING BOUNDARY",
      "source_coordinates":fixed["source_coordinates"],"vertex_count":len(vertices),"electrical_incidence_count":len(edges),"readability_clearance":geometry.READABILITY_CLEARANCE,
      "corridor":fixed["corridor"],"centerline":fixed["centerline"],"movement_eligibility":{"rmin1_side":fixed["rmin1_side_movable_vertices"],"j_out_side":fixed["j_out_side_movable_vertices"]},
      "side_a_available_space":{"direction":list(n),"D_A":da,"delta_A":delta_a,"limiting_record":limit_a,"per_vertex_rays":rays_a},
      "side_b_available_space":{"direction":list(minus),"D_B":db,"delta_B":delta_b,"limiting_record":limit_b,"per_vertex_rays":rays_b},
      "operation_applied":applicable,"coordinates_unchanged":final==before,"coordinates":{v:list(p) for v,p in sorted(final.items())},"edge_reconnection":reconnect,
      "baseline_crowding_events":fixed["before_crowding_events"],"baseline_crowding_clusters":fixed["before_crowding_clusters"],
      "baseline_metrics":fixed["before_metrics"],"fixed_delta_0_5_metrics":fixed["after_metrics"],"half_available_metrics":final_metrics if applicable else None,
      "baseline_boundary_separation":fixed["boundary_separation_before"],"fixed_delta_0_5_boundary_separation":fixed["boundary_separation_after"],"half_available_boundary_separation":corridor.boundary_metrics(fixed["corridor"]["rmin1_side_boundary"]["ordered_vertices"],fixed["corridor"]["j_out_side_boundary"]["ordered_vertices"],final) if applicable else None,
      "interpretation":"D. NO ADVANTAGE — RULE UNDEFINED ON UNBOUNDED SIDE B","reason":"The J_OUT-side movable structure casts into the unbounded exterior; no first unrelated graph/face boundary exists, so D_B and delta_B are not finite geometric quantities.",
      "comparison":{"original_corridor_geometry":fixed["before_metrics"],"fixed_delta_0_5":fixed["after_metrics"],"half_available_space":None},
      "graph_was_not_modified":not applicable}
    output_dir.mkdir(parents=True,exist_ok=True)
    base={"vertices":vertices,"edges":edges,"before":before,"fixed_after":fixed_after};render(base,fixed,report,output_dir/'iamp_corridor_half_available.svg')
    (output_dir/'iamp_corridor_half_available_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--fixed-report',type=Path,default=b/'iamp_corridor_separation_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.fixed_report,a.output)
    print(r['status']);print('D_A',r['side_a_available_space']['D_A'],'delta_A',r['side_a_available_space']['delta_A']);print('D_B',r['side_b_available_space']['D_B'],'delta_B',r['side_b_available_space']['delta_B']);print(r['reason'])

if __name__=='__main__':main()
