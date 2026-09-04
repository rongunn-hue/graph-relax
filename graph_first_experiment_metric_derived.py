#!/usr/bin/env python3
"""Metric-derived polygon scaling plus aggregate-local direct crowding gate."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_common_metric as cm
import graph_first_experiment_common_metric_guard as guard
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as polygon
from render_iamp_common_metric_direct_clean import clean_svg

TOL=1e-7

def aggregate_gate(co,candidate,edges,faces,movers,target_key,target_before,target_after,exclude=()):
 ok,diag=cm.hard_valid(candidate,edges);result={'graph_valid':ok,'target_before':target_before,'target_after':target_after,'target_improved':target_after>target_before+TOL,'P_before':None,'P_after':None,'penalty_change':None,'diagnostics':diag}
 if not ok:return False,'REJECTED_GRAPH_INVALID',result
 snap=guard.affected_snapshot(co,edges,faces,movers);excluded=set(exclude);pb=pa=0.;terms=[]
 for k,b in snap.items():
  if k==target_key or k in excluded:continue
  a=guard.rel_value(k,candidate,faces);ra=a/b['D'];before=max(0.,1-b['r']);after=max(0.,1-ra);pb+=before;pa+=after
  if abs(after-before)>TOL:terms.append({'relationship':str(k),'d_before':b['d'],'d_after':a,'D':b['D'],'r_before':b['r'],'r_after':ra,'penalty_before':before,'penalty_after':after})
 result.update({'P_before':pb,'P_after':pa,'penalty_change':pa-pb,'changed_penalty_terms':terms})
 if not result['target_improved']:return False,'REJECTED_TARGET_NOT_IMPROVED',result
 if not pa<pb-TOL:return False,'REJECTED_AGGREGATE_CROWDING',result
 return True,'ACCEPTED',result

def roots_for_fixed(pole,v,a,D):
 q=(v[0]-pole[0],v[1]-pole[1]);b=(pole[0]-a[0],pole[1]-a[1]);A=q[0]*q[0]+q[1]*q[1];B=2*(b[0]*q[0]+b[1]*q[1]);C=b[0]*b[0]+b[1]*b[1]-D*D;disc=B*B-4*A*C
 if A<=TOL or disc<0:return []
 root=math.sqrt(max(0,disc));return sorted(x for x in ((-B-root)/(2*A),(-B+root)/(2*A)) if x>=1-TOL)

def metric_polygon_candidate(co,edges,faces,fid):
 f=faces[fid];w=f['ordered_boundary_walk'];floating=set(f['free_vertices']);fixed=set(w)-floating;pole=face_poles.safe_pole([co[v] for v in w])[0];metrics=polygon.metrics(w,co);D=metrics['average_nearest_boundary_node_distance'];nearest=[];seen=set()
 for v in sorted(set(w)):
  u=min((u for u in set(w) if u!=v),key=lambda u:(math.dist(co[v],co[u]),u));pair=tuple(sorted((v,u)))
  if pair in seen:continue
  seen.add(pair);d=math.dist(co[pair[0]],co[pair[1]])
  if d>=D-TOL:continue
  a,b=pair;kind=('MOVABLE' if a in floating else 'FIXED')+'_'+('MOVABLE' if b in floating else 'FIXED');s=None;reason=None
  if a in floating and b in floating:s=D/d
  elif a not in floating and b not in floating:reason='FIXED_FIXED_UNCORRECTABLE'
  else:
   v=a if a in floating else b;anchor=b if v==a else a;roots=roots_for_fixed(pole,co[v],co[anchor],D)
   if roots:s=roots[0]
   else:reason='NO_REAL_ROOT_AT_S_GTE_1'
  nearest.append({'pair':list(pair),'kind':kind,'d':d,'D':D,'r':d/D,'s_pair':s,'uncorrectable_reason':reason})
 valid_scales=[x['s_pair'] for x in nearest if x['s_pair'] is not None];s=max(valid_scales) if valid_scales else None
 candidate=dict(co)
 if s is not None:
  for v in floating:candidate[v]=(pole[0]+s*(co[v][0]-pole[0]),pole[1]+s*(co[v][1]-pole[1]))
 before_area=abs(polygon.area(w,co));after_area=abs(polygon.area(w,candidate));after_metrics=polygon.metrics(w,candidate);exclude=[('PNN',fid,v) for v in set(w)]
 if s is None:accepted=False;decision='REJECTED_NO_CORRECTABLE_DEFICIENT_RELATIONSHIP';gate={}
 else:
  accepted,decision,gate=aggregate_gate(co,candidate,edges,faces,floating,None,metrics['minimum_nearest_boundary_node_distance']/D,after_metrics['minimum_nearest_boundary_node_distance']/D,exclude)
  if accepted and (after_area<=before_area+TOL or after_metrics['minimum_nearest_boundary_node_distance']<=metrics['minimum_nearest_boundary_node_distance']+TOL):accepted=False;decision='REJECTED_TARGET_NOT_IMPROVED'
 return (candidate if accepted else co),{'polygon_id':fid,'pole':list(pole),'floating_vertices':sorted(floating),'fixed_vertices':sorted(fixed),'deficient_NN_relationships':nearest,'s_polygon':s,'area_before':before_area,'area_after_candidate':after_area,'minimum_NN_before':metrics['minimum_nearest_boundary_node_distance'],'minimum_NN_after_candidate':after_metrics['minimum_nearest_boundary_node_distance'],'average_NN_before':D,'average_NN_after_candidate':after_metrics['average_nearest_boundary_node_distance'],'decision':'ACCEPTED' if accepted else decision,'gate':gate,'candidate_count':1 if s is not None else 0}

def polygon_pass(co,edges,faces,event):
 target=event['face_id'];rank=polygon.ranking(list(faces.values()),co);_,_,paths,order=polygon.dependencies(list(faces.values()),rank,[target]);ops=[]
 for fid in order:co,o=metric_polygon_candidate(co,edges,faces,fid);ops.append(o)
 return co,{'target':target,'dependency_paths':paths,'processing_order':order,'operations':ops}

def summarize(ops):
 candidates=[]
 for o in ops:candidates+=o.get('candidate_results',[o])
 return {'attempted':len(candidates),'accepted':sum(x.get('decision')=='ACCEPTED' for x in candidates),'rejected_target_not_improved':sum(x.get('decision')=='REJECTED_TARGET_NOT_IMPROVED' for x in candidates),'rejected_aggregate_crowding':sum(x.get('decision')=='REJECTED_AGGREGATE_CROWDING' for x in candidates),'rejected_graph_validity':sum(x.get('decision')=='REJECTED_GRAPH_INVALID' for x in candidates),'skipped_no_longer_crowded':sum(x.get('decision')=='SKIPPED_NO_LONGER_CROWDED' for x in candidates),'other_rejections':sum(str(x.get('decision','')).startswith('REJECTED') and x.get('decision') not in ('REJECTED_TARGET_NOT_IMPROVED','REJECTED_AGGREGATE_CROWDING','REJECTED_GRAPH_INVALID') for x in candidates)}

def global_penalty(m):return sum(max(0.,1-x['ratio']) for rows in m.values() for x in rows)

def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());V={x['id']:x for x in gd['vertices']};E=[(x['source'],x['target']) for x in gd['edges']];pairs=guard.edgeset(E);start={v:tuple(p) for v,p in dr['final_coordinates'].items()};ok,sd=cm.hard_valid(start,E)
 if not ok or (len(V),len(E))!=(34,44):raise AssertionError('invalid authoritative start')
 fl,darts=inventory.enumerate_faces(V,E,start);inventory.add_classification_and_adjacency(fl,E,darts);faces={f['face_id']:f for f in fl};initial=cm.corrected_map(fl,start,E)
 # Reuse exact direct node/edge procedures; only their gate is replaced.
 original_gate=guard.gate;guard.gate=aggregate_gate
 try:
  co,pstage=polygon_pass(dict(start),E,faces,initial['polygon'][0]);polygon_co=dict(co);after_polygon=cm.corrected_map(fl,co,E);co,nops=guard.node_pass(co,E,faces,after_polygon);node_co=dict(co);after_nodes=cm.corrected_map(fl,co,E);co,eops=guard.edge_pass(co,E,faces,after_nodes)
 finally:guard.gate=original_gate
 final=cm.corrected_map(fl,co,E);ok,fd=cm.hard_valid(co,E);assert ok and pairs==guard.edgeset(E)
 stages={'start':start,'polygon':polygon_co,'nodes':node_co,'final':co};outdir.mkdir(parents=True,exist_ok=True);svgs=[]
 for name,c in stages.items():p=outdir/f'iamp_metric_derived_{name}.svg';clean_svg(V,E,c,f'IAMP METRIC-DERIVED — {name.upper()}',p);svgs.append(p.name)
 jeb={k:math.dist(c['component:J_EB'],c['component:C2']) for k,c in stages.items()};audit=[]
 for stage,ops in [('node',nops),('edge',eops)]:
  for o in ops:
   for c in o.get('candidate_results',[o]):
    mover=c.get('mover',o.get('selected_mover'))
    if mover in ('component:J_EB','component:C2') or ('component:J_EB' in str(o.get('offending_geometry',o)) or 'component:C2' in str(o.get('offending_geometry',o))):
     g=c.get('gate',o.get('gate',{}));audit.append({'stage':stage,'target_event':o.get('event_id'),'mover':mover,'target_r_before':o.get('r_before'),'target_r_after':c.get('r_after',o.get('r_after_candidate')),'P_before':g.get('P_before'),'P_after':g.get('P_after'),'decision':c.get('decision',o.get('decision'))})
 report={'schema':'graph-relax.metric-derived-aggregate.v1','source_coordinates':f'{direct_path}#final_coordinates','starting_validation':{'vertices':34,'edges':44,'endpoint_pairs_identical':pairs==guard.edgeset(E),'proper_crossings':sd['proper_unrelated_edge_crossing_count'],'coincidences':sd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':sd['vertex_on_unrelated_edge_interior_count']},'initial_map':initial,'start_counts':{k:len(v) for k,v in initial.items()},'start_min_ratio':{k:min((x['ratio'] for x in v),default=None) for k,v in initial.items()},'start_detector_penalty':global_penalty(initial),'polygon_stage':{**pstage,'summary':summarize(pstage['operations']),'map_after':after_polygon},'node_stage':{'operations':nops,'summary':summarize(nops),'map_after':after_nodes},'edge_stage':{'operations':eops,'summary':summarize(eops)},'J_EB_C2_distances':jeb,'J_EB_C2_candidate_audit':audit,'final_map':final,'final_counts':{k:len(v) for k,v in final.items()},'final_min_ratio':{k:min((x['ratio'] for x in v),default=None) for k,v in final.items()},'final_detector_penalty':global_penalty(final),'final_validation':{'vertices':34,'edges':44,'endpoint_pairs_identical':pairs==guard.edgeset(E),'proper_crossings':fd['proper_unrelated_edge_crossing_count'],'coincidences':fd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':fd['vertex_on_unrelated_edge_interior_count']},'svg_files':svgs,'stage_coordinates':{k:{v:list(p) for v,p in sorted(c.items())} for k,c in stages.items()}}
 (outdir/'iamp_metric_derived_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({k:r[k] for k in ('start_counts','start_min_ratio','start_detector_penalty','J_EB_C2_distances','final_counts','final_min_ratio','final_detector_penalty','final_validation')},indent=2));print('polygon',r['polygon_stage']['summary']);print('node',r['node_stage']['summary']);print('edge',r['edge_stage']['summary'])
if __name__=='__main__':main()
