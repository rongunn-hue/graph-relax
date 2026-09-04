#!/usr/bin/env python3
"""Direct analytic IAMP uncrossing with a hard meaningful node spacing."""
import argparse,json,math
from itertools import combinations
from pathlib import Path
import graph_first_experiment_endpoint_rotation as loader
import graph_first_experiment_rotation_direct as direct
import placement_geometry as geometry

S_MIN=12.0
SANITATION_LIMIT=100000

def pair_distances(co):
 return sorted(({"vertices":[a,b],"distance":math.dist(co[a],co[b])} for a,b in combinations(sorted(co),2)),key=lambda x:(x["distance"],x["vertices"]))

def sanitize(co,minimum=S_MIN):
 p=dict(co);initial=pair_distances(p);corrections=[];nodes=sorted(p);pairs=list(combinations(nodes,2))
 for iteration in range(SANITATION_LIMIT):
  violation=None
  for pair_index,(a,b) in enumerate(pairs):
   d=math.dist(p[a],p[b])
   if d+geometry.GEOMETRY_TOLERANCE<minimum:violation=(pair_index,a,b,d);break
  if violation is None:break
  pair_index,a,b,d=violation;needed=minimum-d
  if d<=geometry.GEOMETRY_TOLERANCE:
   angle=2*math.pi*pair_index/len(pairs);ux,uy=math.cos(angle),math.sin(angle);method="stable-pair-order direction"
  else:ux=(p[b][0]-p[a][0])/d;uy=(p[b][1]-p[a][1])/d;method="joining line"
  half=needed/2;p[a]=(p[a][0]-half*ux,p[a][1]-half*uy);p[b]=(p[b][0]+half*ux,p[b][1]+half*uy)
  corrections.append({"iteration":iteration+1,"pair":[a,b],"distance_before":d,"required_separation_added":needed,"movement_each":needed/2,"direction":[ux,uy],"direction_method":method,"coordinates_after":{a:list(p[a]),b:list(p[b])}})
 else:raise AssertionError("minimum-spacing sanitation safety limit exceeded")
 final=pair_distances(p)
 if final[0]["distance"]+geometry.GEOMETRY_TOLERANCE<minimum:raise AssertionError("sanitation failed")
 return p,{"initial_minimum_node_distance":initial[0]["distance"],"initial_violating_pair_count":sum(x["distance"]<minimum for x in initial),"closest_initial_pairs":initial[:12],"corrections_made":len(corrections),"corrections":corrections,"iterations":len(corrections),"iteration_safety_limit":SANITATION_LIMIT,"post_sanitation_minimum_node_distance":final[0]["distance"],"closest_post_sanitation_pairs":final[:12]}

def panel(graph,co,title,xoff,view):
 minx,miny,w,h=view;z=[f'<g transform="translate({xoff:.17g} 0)"><rect x="{minx:.17g}" y="{miny:.17g}" width="{w:.17g}" height="{h:.17g}" fill="white" stroke="#bbb"/>',f'<text x="{minx+12:.17g}" y="{miny+25:.17g}" font-family="sans-serif" font-size="16" font-weight="bold">{title}</text>']
 for r in direct.rendered_edge_records(graph,co):
  a,b=r["segment_start"],r["segment_end"];z.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#666" stroke-width="1.5"/>')
 for v in sorted(graph):
  x,y=co[v];color="#1565c0" if graph.nodes[v]["type"]=="COMPONENT" else "#d84315";z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{graph.nodes[v]["name"]}</text>')
 z.append('</g>');return z

def render(graph,raw,sanitized,final,path):
 view=loader.viewbox_for([raw,sanitized,final]);minx,miny,w,h=view;gap=30;z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {3*w+2*gap:.17g} {h:.17g}" width="2100" height="800">'];z+=panel(graph,raw,"RAW PRE-UNCROSSING MDS",0,view);z+=panel(graph,sanitized,"POST S_MIN=12 SANITATION",w+gap,view);z+=panel(graph,final,"FINAL ANALYTIC UNCROSSING",2*(w+gap),view);z.append('</svg>');path.write_text('\n'.join(z)+'\n')

def run(graph_path,coordinate_path,mds_report_path,old_report_path,outdir):
 graph,raw,raw_diag,_=loader.load_start(graph_path,coordinate_path,mds_report_path);old=json.loads(old_report_path.read_text());sanitized,san=sanitize(raw);san_diag=geometry.graph_geometry_diagnostics(sanitized,graph.edges)
 final,final_diag,moves,states,history,counters=direct.direct_search(graph,sanitized,S_MIN)
 orientation_debug=[]
 for audit in counters["orientation_audits"]:
  for s in audit["solutions"]:
   spacing=s.get("spacing_intervals",[]);orientation_debug.append({"crossing_being_repaired":audit["crossing"],"moving_endpoint":s.get("moving"),"pivot":s.get("pivot"),"radius":s.get("radius"),"node_spacing_forbidden_interval_count":sum("interval" in x for x in spacing),"full_circle_infeasible_vertices":[x["vertex"] for x in spacing if x["status"]=="ENTIRE"],"crossing_forbidden_interval":s.get("target_interval"),"merged_forbidden_interval_set":s.get("union_segments_shifted"),"current_theta":s.get("theta_old"),"selected_nearest_legal_theta":s.get("theta_new"),"angular_displacement":s.get("signed_delta_radians"),"resulting_minimum_node_distance":s.get("resulting_minimum_node_distance"),"resulting_crossing_count":s.get("resulting_total_crossings"),"eligible":s.get("selection_status") in {"ELIGIBLE_STRICT_REDUCTION","SELECTED","ACCEPTED_GLOBAL_MOVE"},"reason":s.get("reason") or s.get("selection_status")})
 oldseq=[{"moving_vertex":m["moving_vertex"],"pivot_vertex":m["pivot_vertex"],"crossings_before":m["crossing_count_before"],"crossings_after":m["crossing_count_after"]} for m in old["moves"]];newseq=[{"moving_vertex":m["moving_vertex"],"pivot_vertex":m["pivot_vertex"],"crossings_before":m["crossing_count_before"],"crossings_after":m["crossing_count_after"]} for m in moves]
 report={"schema":"graph-relax.direct-uncross-min-spacing.v1","pre_uncrossing_coordinate_source":str(coordinate_path),"source_description":"authoritative spaced nearness-MDS coordinates immediately before direct analytic repairs","S_MIN":S_MIN,"S_MIN_rationale":"single first meaningful spacing test at the current coordinate and rendering scale","vertex_count":graph.number_of_nodes(),"electrical_edge_count":graph.number_of_edges(),"spacing_scope":{"angular_exclusion":"all vertices except moving endpoint and invariant fixed pivot","global_audit":"every distinct graph-vertex pair, including adjacent and moving-endpoint/pivot pairs"},"raw_start":{"crossing_count":raw_diag["proper_unrelated_edge_crossing_count"],"minimum_node_distance":raw_diag["minimum_node_distance"],"closest_node_pairs":pair_distances(raw)[:12]},"sanitation":{**san,"all_distinct_pairs_audited":len(pair_distances(sanitized)),"expected_distinct_pair_count":graph.number_of_nodes()*(graph.number_of_nodes()-1)//2,"crossing_count_after":san_diag["proper_unrelated_edge_crossing_count"],"coordinates":{v:list(sanitized[v]) for v in sorted(sanitized)}},"uncrossing":{"starting_crossing_count":san_diag["proper_unrelated_edge_crossing_count"],"final_crossing_count":final_diag["proper_unrelated_edge_crossing_count"],"accepted_move_count":len(moves),"moves":moves,"crossing_history":history,"orientation_debug":orientation_debug,"all_orientation_outcomes":counters["all_solution_outcomes"],"maximum_positions_per_orientation":max((x["positions_evaluated"] for x in counters["positions_per_choice"]),default=0),"maximum_orientations_per_crossing":max((x["calculated"] for x in counters["solutions_per_crossing"]),default=0),"moving_pivot_pairs_globally_valid_after_every_move":all(m["preserved_edge_length_before"]+geometry.GEOMETRY_TOLERANCE>=S_MIN for m in moves),"S_MIN_preserved_after_every_accepted_move":all(m["minimum_node_distance_after"]+geometry.GEOMETRY_TOLERANCE>=S_MIN for m in moves)},"final":{"coordinates":{v:list(final[v]) for v in sorted(final)},"minimum_node_distance":final_diag["minimum_node_distance"],"proper_crossings":final_diag["proper_unrelated_edge_crossing_count"],"coincidences":final_diag["coincident_vertex_pair_count"],"vertex_on_unrelated_edge_events":final_diag["vertex_on_unrelated_edge_interior_count"],"all_distinct_pair_count_audited":len(pair_distances(final)),"zero_crossings_achieved":final_diag["proper_unrelated_edge_crossing_count"]==0},"old_comparison":{"old_final_crossing_count":old["final_crossing_count"],"new_final_crossing_count":final_diag["proper_unrelated_edge_crossing_count"],"old_final_minimum_node_distance":old["final_geometry_diagnostics"]["minimum_node_distance"],"new_final_minimum_node_distance":final_diag["minimum_node_distance"],"old_accepted_sequence":oldseq,"new_accepted_sequence":newseq,"same_four_repairs":oldseq==newseq}}
 outdir.mkdir(parents=True,exist_ok=True);render(graph,raw,sanitized,final,outdir/"iamp_uncross_with_min_spacing.svg");(outdir/"iamp_uncross_with_min_spacing_report.json").write_text(json.dumps(report,indent=2,sort_keys=True,default=str)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--coordinates",type=Path,default=b/"iamp_nearness_mds.json");p.add_argument("--mds-report",type=Path,default=b/"iamp_nearness_mds_report.json");p.add_argument("--old-report",type=Path,default=b/"iamp_nearness_rotation_direct_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.coordinates,a.mds_report,a.old_report,a.output);print("raw",r["raw_start"]);print("sanitation",{k:v for k,v in r["sanitation"].items() if k not in {"corrections","coordinates"}});print("final",r["final"]);print("moves",r["old_comparison"])
if __name__=="__main__":main()
