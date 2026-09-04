#!/usr/bin/env python3
"""Direct common-metric IAMP pass with local crowding non-regression gates."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_common_metric as cm
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as polygon
import placement_geometry as geometry
from render_iamp_common_metric_direct_clean import clean_svg

TOL=1e-7
def edge(a,b):return tuple(sorted((a,b)))
def edgeset(edges):return sorted(edge(a,b) for a,b in edges)
def overlay_distance(co,v,a,b):return min(geometry.point_segment_distance(co[a],co[v],co[b]),geometry.point_segment_distance(co[b],co[v],co[a]))

def rel_value(key,co,faces):
 t=key[0]
 if t=='NN':return math.dist(co[key[1]],co[key[2]])
 if t=='NE':return geometry.point_segment_distance(co[key[1]],co[key[2][0]],co[key[2][1]])
 if t=='EE':p,q=cm.base.segcp(co[key[1][0]],co[key[1][1]],co[key[2][0]],co[key[2][1]]);return math.dist(p,q)
 if t=='SE':return overlay_distance(co,key[1],key[2],key[3])
 if t=='PNN':
  w=faces[key[1]]["ordered_boundary_walk"];v=key[2];return min(math.dist(co[v],co[q]) for q in set(w) if q!=v)
 raise AssertionError(key)

def affected_snapshot(co,edges,faces,movers):
 es=edgeset(edges);movers=set(movers);inc={e for e in es if set(e)&movers};keys=set()
 for m in movers:
  for q in co:
   if q!=m:keys.add(('NN',*sorted((m,q))))
  for e in es:
   if m not in e:keys.add(('NE',m,e))
 # Other nodes against an incident moving edge also change geometrically.
 for e in inc:
  for v in co:
   if v not in e:keys.add(('NE',v,e))
  for f in es:
   if not set(e)&set(f):keys.add(('EE',*sorted((e,f))))
 # All shared overlays containing any moved incident edge.
 adjacency={v:[] for v in co}
 for a,b in es:adjacency[a].append(b);adjacency[b].append(a)
 for v in adjacency:
  for i,a in enumerate(sorted(adjacency[v])):
   for b in sorted(adjacency[v])[i+1:]:
    if edge(v,a) in inc or edge(v,b) in inc:keys.add(('SE',v,a,b))
 for fid,f in faces.items():
  if f['is_bounded'] and f['simple_boundary'] and movers&set(f['ordered_boundary_walk']):
   for v in sorted(set(f['ordered_boundary_walk'])):keys.add(('PNN',fid,v))
 out={}
 for k in sorted(keys,key=str):
  if k[0]=='NN':participants=k[1:]
  elif k[0]=='NE':participants=[k[1],*k[2]]
  elif k[0]=='EE':participants=[*k[1],*k[2]]
  elif k[0]=='SE':participants=k[1:]
  else:participants=list(set(faces[k[1]]['ordered_boundary_walk']))
  D,_=cm.local_d(co,participants) if k[0]!='PNN' else (polygon.metrics(faces[k[1]]['ordered_boundary_walk'],co)['average_nearest_boundary_node_distance'],None)
  d=rel_value(k,co,faces);out[k]={'d':d,'D':D,'r':d/D,'crowded':d<D-TOL}
 return out

def gate(co,candidate,edges,faces,movers,target_key,target_before,target_after,exclude=()):
 ok,diag=cm.hard_valid(candidate,edges)
 result={'graph_valid':ok,'target_before':target_before,'target_after':target_after,'target_improved':target_after>target_before+TOL,'new_crowding':[],'existing_crowding_worsened':[],'diagnostics':diag}
 if not ok:return False,'REJECTED_GRAPH_INVALID',result
 snap=affected_snapshot(co,edges,faces,movers);excluded=set(exclude)
 for k,b in snap.items():
  if k==target_key or k in excluded:continue
  a=rel_value(k,candidate,faces);ra=a/b['D']
  row={'relationship':str(k),'d_before':b['d'],'d_after':a,'D':b['D'],'r_before':b['r'],'r_after':ra}
  if b['r']>=1-TOL and ra<1-TOL:result['new_crowding'].append(row)
  elif b['r']<1-TOL and ra<b['r']-TOL:result['existing_crowding_worsened'].append(row)
 if not result['target_improved']:return False,'REJECTED_TARGET_NOT_IMPROVED',result
 if result['new_crowding'] or result['existing_crowding_worsened']:return False,'REJECTED_NON_REGRESSION',result
 return True,'ACCEPTED',result

def current_target(e,co):
 if e['type']=='NODE_NODE':
  a,b=e['affected'];D,_=cm.local_d(co,[a,b]);return math.dist(co[a],co[b]),D,('NN',*sorted((a,b)))
 if e['type']=='NODE_EDGE':
  v=e['affected']['vertex'];a,b=e['affected']['edge'];D,_=cm.local_d(co,[v,a,b]);return geometry.point_segment_distance(co[v],co[a],co[b]),D,('NE',v,edge(a,b))
 if e['type']=='EDGE_EDGE':
  e1,e2=map(tuple,e['affected']);p,q=cm.base.segcp(co[e1[0]],co[e1[1]],co[e2[0]],co[e2[1]]);D,_=cm.local_d(co,[*e1,*e2]);return math.dist(p,q),D,('EE',*sorted((edge(*e1),edge(*e2))))
 v=e['shared_vertex'];a,b=e['nonshared_endpoints'];D,_=cm.local_d(co,[v,a,b]);return overlay_distance(co,v,a,b),D,('SE',v,a,b)

def polygon_pass(co,edges,faces,event):
 target=event['face_id'];rank=polygon.ranking(list(faces.values()),co);_,_,paths,order=polygon.dependencies(list(faces.values()),rank,[target]);ops=[]
 for fid in order:
  f=faces[fid];w=f['ordered_boundary_walk'];floating=f['free_vertices'];bm=polygon.metrics(w,co);ba=abs(polygon.area(w,co));dirs=polygon.directions(w,co,floating) if floating else {};t=bm['minimum_nearest_boundary_node_distance'];candidate=polygon.proposed(co,dirs,t);am=abs(polygon.area(w,candidate));mm=polygon.metrics(w,candidate)
  # Target is the polygon minimum NN; its individual polygon-NN relationships are target scope.
  excludes=[('PNN',fid,v) for v in set(w)];valid,decision,g=gate(co,candidate,edges,faces,floating,None,bm['minimum_nearest_boundary_node_distance']/bm['average_nearest_boundary_node_distance'],mm['minimum_nearest_boundary_node_distance']/bm['average_nearest_boundary_node_distance'],excludes)
  shape_ok=polygon.shell_simple(w,candidate) and am>ba+TOL and mm['minimum_nearest_boundary_node_distance']>bm['minimum_nearest_boundary_node_distance']+TOL and mm['average_nearest_boundary_node_distance']>bm['average_nearest_boundary_node_distance']+TOL
  if valid and shape_ok:co=candidate
  elif valid:valid=False;decision='REJECTED_TARGET_NOT_IMPROVED'
  ops.append({'face_id':fid,'candidate_count':1 if floating else 0,'movers':floating,'target_displacement':t,'decision':decision,'area_before':ba,'area_after_candidate':am,'metrics_before':bm,'metrics_after_candidate':mm,'gate':g})
 return co,{'target':target,'dependency_paths':paths,'processing_order':order,'operations':ops}

def node_pass(co,edges,faces,m):
 ops=[]
 for e in sorted(m['node_node']+m['node_edge'],key=lambda x:(x['ratio'],x['event_id'])):
  d,D,key=current_target(e,co)
  if d>=D-TOL:ops.append({'event_id':e['event_id'],'decision':'SKIPPED_NO_LONGER_CROWDED','current_d':d,'current_D':D});continue
  if e['type']=='NODE_NODE':
   a,b=e['affected'];deg={z:sum(z in x for x in edges) for z in (a,b)};mover=sorted((a,b),key=lambda z:(deg[z],z))[0];fixed=b if mover==a else a;u=cm.base.unit(co[mover],co[fixed],mover);target=(co[fixed][0]+D*u[0],co[fixed][1]+D*u[1]);score=lambda c:math.dist(c[mover],c[fixed])
  else:
   mover=e['affected']['vertex'];a,b=e['affected']['edge'];q=cm.base.cp(co[mover],co[a],co[b]);u=cm.base.unit(co[mover],q,mover);target=(q[0]+D*u[0],q[1]+D*u[1]);score=lambda c:geometry.point_segment_distance(c[mover],c[a],c[b])
  before_coordinate=co[mover];cand=dict(co);cand[mover]=target;after=score(cand);accepted,decision,g=gate(co,cand,edges,faces,[mover],key,d/D,after/D)
  jeb_before=math.dist(co['component:J_EB'],co['component:C2']);jeb_after=math.dist(cand['component:J_EB'],cand['component:C2'])
  if accepted:co=cand
  ops.append({'event_id':e['event_id'],'event_type':e['type'],'selected_mover':mover,'calculated_target':list(target),'displacement':math.dist(before_coordinate,target),'d_before':d,'D':D,'r_before':d/D,'d_after_candidate':after,'r_after_candidate':after/D,'decision':decision,'gate':g,'J_EB_C2_before':jeb_before,'J_EB_C2_candidate':jeb_after})
 return co,ops

def edge_pass(co,edges,faces,m):
 ops=[]
 for e in sorted(m['unrelated_edge_edge']+m['shared_endpoint_overlay'],key=lambda x:(x['ratio'],x['event_id'])):
  d,D,key=current_target(e,co)
  if d>=D-TOL:ops.append({'event_id':e['event_id'],'decision':'SKIPPED_NO_LONGER_CROWDED','current_d':d,'current_D':D});continue
  candidates=[]
  if e['type']=='EDGE_EDGE':
   e1,e2=e['affected'];p,q=cm.base.segcp(co[e1[0]],co[e1[1]],co[e2[0]],co[e2[1]])
   def score(c):x,y=cm.base.segcp(c[e1[0]],c[e1[1]],c[e2[0]],c[e2[1]]);return math.dist(x,y)
   for mover in e1+e2:
    ref=q if mover in e1 else p;u=cm.base.unit(co[mover],ref,mover);c=dict(co);c[mover]=(ref[0]+D*u[0],ref[1]+D*u[1]);candidates.append((mover,c,score(c)))
  else:
   v=e['shared_vertex'];a,b=e['nonshared_endpoints'];score=lambda c:overlay_distance(c,v,a,b)
   for mover,other in ((a,b),(b,a)):
    q=cm.base.cp(co[mover],co[v],co[other]);u=cm.base.unit(co[mover],q,mover);c=dict(co);c[mover]=(q[0]+D*u[0],q[1]+D*u[1]);candidates.append((mover,c,score(c)))
  results=[]
  for mover,c,after in candidates:
   accepted,decision,g=gate(co,c,edges,faces,[mover],key,d/D,after/D);results.append({'mover':mover,'target':list(c[mover]),'displacement':math.dist(co[mover],c[mover]),'d_after':after,'r_after':after/D,'decision':decision,'gate':g,'coordinates':c,'J_EB_C2_before':math.dist(co['component:J_EB'],co['component:C2']),'J_EB_C2_candidate':math.dist(c['component:J_EB'],c['component:C2'])})
  legal=[x for x in results if x['decision']=='ACCEPTED']
  if legal:
   chosen=max(legal,key=lambda x:(x['r_after'],x['mover']));co=chosen['coordinates'];decision='ACCEPTED';selected=chosen['mover']
  else:decision='REJECTED_ALL_CANDIDATES';selected=None
  for x in results:x.pop('coordinates')
  ops.append({'event_id':e['event_id'],'event_type':e['type'],'d_before':d,'D':D,'r_before':d/D,'candidate_count':len(results),'candidate_results':results,'selected_mover':selected,'decision':decision})
 return co,ops

def counts(ops):
 flat=[]
 for o in ops:
  if 'candidate_results' in o:flat+=o['candidate_results']
  else:flat.append(o)
 return {'attempted':len(flat),'accepted':sum(x.get('decision')=='ACCEPTED' for x in flat),'rejected_graph_validity':sum(x.get('decision')=='REJECTED_GRAPH_INVALID' for x in flat),'rejected_non_regression':sum(x.get('decision')=='REJECTED_NON_REGRESSION' for x in flat),'rejected_other':sum(str(x.get('decision','')).startswith('REJECTED') and x.get('decision') not in ('REJECTED_GRAPH_INVALID','REJECTED_NON_REGRESSION') for x in flat)}

def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());V={x['id']:x for x in gd['vertices']};E=[(x['source'],x['target']) for x in gd['edges']];start={v:tuple(p) for v,p in dr['final_coordinates'].items()};start_pairs=edgeset(E);ok,sd=cm.hard_valid(start,E)
 if not ok or (len(V),len(E))!=(34,44):raise AssertionError('invalid authoritative start')
 faces_list,darts=inventory.enumerate_faces(V,E,start);inventory.add_classification_and_adjacency(faces_list,E,darts);faces={f['face_id']:f for f in faces_list};initial=cm.corrected_map(faces_list,start,E);co,pstage=polygon_pass(dict(start),E,faces,initial['polygon'][0]);polygon_co=dict(co);after_polygon=cm.corrected_map(faces_list,co,E);co,nops=node_pass(co,E,faces,after_polygon);node_co=dict(co);after_nodes=cm.corrected_map(faces_list,co,E);co,eops=edge_pass(co,E,faces,after_nodes);final=cm.corrected_map(faces_list,co,E);ok,fd=cm.hard_valid(co,E)
 assert ok and start_pairs==edgeset(E)
 outdir.mkdir(parents=True,exist_ok=True);stage={'start':start,'polygon':polygon_co,'nodes':node_co,'final':co};svgs=[]
 for name,c in stage.items():p=outdir/f'iamp_common_metric_guard_{name}.svg';clean_svg(V,E,c,f'IAMP COMMON-METRIC GUARD — {name.upper()}',p);svgs.append(p.name)
 jeb={k:math.dist(c['component:J_EB'],c['component:C2']) for k,c in stage.items()};related=[]
 for section,ops in [('node',nops),('edge',eops)]:
  for o in ops:
   rs=o.get('candidate_results',[o])
   for x in rs:
    if x.get('decision')=='REJECTED_NON_REGRESSION' and ('component:J_EB' in str(o) or 'component:C2' in str(o) or abs(x.get('J_EB_C2_candidate',jeb['final'])-x.get('J_EB_C2_before',jeb['final']))>TOL):related.append({'stage':section,'target_event':o.get('event_id'),'candidate':x})
 report={'schema':'graph-relax.common-metric-guard.v1','source_coordinates':f'{direct_path}#final_coordinates','connectivity_endpoint_pairs':[list(x) for x in start_pairs],'starting_validation':{'vertices':34,'edges':44,'proper_crossings':sd['proper_unrelated_edge_crossing_count'],'coincidences':sd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':sd['vertex_on_unrelated_edge_interior_count']},'initial_event_counts':{k:len(v) for k,v in initial.items()},'initial_map':initial,'polygon_stage':{**pstage,'summary':counts(pstage['operations']),'map_after':after_polygon},'node_stage':{'operations':nops,'summary':counts(nops),'map_after':after_nodes},'edge_stage':{'operations':eops,'summary':counts(eops)},'J_EB_C2_distances':jeb,'J_EB_C2_related_non_regression_rejections':related,'final_map':final,'final_remaining_counts':{k:len(v) for k,v in final.items()},'final_min_ratio':{k:min((x['ratio'] for x in v),default=None) for k,v in final.items()},'final_validation':{'vertices':34,'edges':44,'endpoint_pairs_identical':start_pairs==edgeset(E),'proper_crossings':fd['proper_unrelated_edge_crossing_count'],'coincidences':fd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':fd['vertex_on_unrelated_edge_interior_count']},'svg_files':svgs,'stage_coordinates':{s:{v:list(p) for v,p in sorted(c.items())} for s,c in stage.items()}}
 (outdir/'iamp_common_metric_guard_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({'initial_event_counts':r['initial_event_counts'],'polygon':r['polygon_stage']['summary'],'node':r['node_stage']['summary'],'edge':r['edge_stage']['summary'],'J_EB_C2':r['J_EB_C2_distances'],'final_counts':r['final_remaining_counts'],'final_min_ratio':r['final_min_ratio'],'final_validation':r['final_validation']},indent=2))
if __name__=='__main__':main()
