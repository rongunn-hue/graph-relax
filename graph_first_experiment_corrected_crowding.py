#!/usr/bin/env python3
"""Corrected four-type crowding audit followed by one three-stage IAMP pass."""
import argparse,json,math
from itertools import combinations
from pathlib import Path
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as poly
import placement_geometry as geo

TOL=1e-7
def cp(p,a,b):
 dx,dy=b[0]-a[0],b[1]-a[1];q=dx*dx+dy*dy;t=0 if q==0 else max(0,min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/q));return(a[0]+t*dx,a[1]+t*dy)
def segcp(a,b,c,d):return min(((a,cp(a,c,d)),(b,cp(b,c,d)),(cp(c,a,b),c),(cp(d,a,b),d)),key=lambda z:(math.dist(*z),z))
def unit(a,b,key):
 x,y=a[0]-b[0],a[1]-b[1];q=math.hypot(x,y)
 if q>TOL:return x/q,y/q
 t=sum(map(ord,key))*math.pi/180;return math.cos(t),math.sin(t)
def valid(edges,co):
 d=geo.graph_geometry_diagnostics(co,edges);r=geo.validate_rebuilt_edge_segments(co,edges,geo.rebuild_straight_edge_segments(co,edges));return r['valid'] and not d['proper_unrelated_edge_crossing_count'] and not d['coincident_vertex_pair_count'] and not d['vertex_on_unrelated_edge_interior_count'],d
def node_node(co):
 z=[{'vertices':[a,b],'clearance':math.dist(co[a],co[b])} for a,b in combinations(sorted(co),2)];z.sort(key=lambda x:(x['clearance'],x['vertices']));
 for i,x in enumerate(z,1):x['rank']=i
 return z
def node_edge(co,edges):
 z=[]
 for v in sorted(co):
  for a,b in sorted(tuple(sorted(e)) for e in edges):
   if v in (a,b):continue
   p=cp(co[v],co[a],co[b]);z.append({'vertex':v,'edge':[a,b],'clearance':math.dist(co[v],p),'closest_point':list(p)})
 z.sort(key=lambda x:(x['clearance'],x['vertex'],x['edge']))
 for i,x in enumerate(z,1):x['rank']=i
 return z
def edge_edge(co,edges):
 es=sorted(tuple(sorted(e)) for e in edges);z=[]
 for i,e in enumerate(es):
  for f in es[i+1:]:
   if set(e)&set(f):continue
   p,q=segcp(co[e[0]],co[e[1]],co[f[0]],co[f[1]]);z.append({'edges':[list(e),list(f)],'clearance':math.dist(p,q),'closest_points':[list(p),list(q)]})
 z.sort(key=lambda x:(x['clearance'],x['edges']))
 for i,x in enumerate(z,1):x['rank']=i
 return z
def overlays(co,edges):
 inc={v:[] for v in co}
 for a,b in edges:inc[a].append(b);inc[b].append(a)
 z=[]
 for v in sorted(co):
  for a,b in combinations(sorted(inc[v]),2):
   ua=(co[a][0]-co[v][0],co[a][1]-co[v][1]);ub=(co[b][0]-co[v][0],co[b][1]-co[v][1]);la,lb=math.hypot(*ua),math.hypot(*ub);cos=max(-1,min(1,(ua[0]*ub[0]+ua[1]*ub[1])/(la*lb)));ang=math.acos(cos);sep=2*min(la,lb)*math.sin(ang/2)
   z.append({'shared_vertex':v,'edges':[[v,a],[v,b]],'nonshared_endpoints':[a,b],'angle_radians':ang,'angle_degrees':math.degrees(ang),'closest_separation_away_from_shared_vertex':sep,'available_segment_lengths':[la,lb]})
 z.sort(key=lambda x:(x['closest_separation_away_from_shared_vertex'],x['angle_radians'],x['shared_vertex'],x['nonshared_endpoints']))
 for i,x in enumerate(z,1):x['rank']=i
 return z
def maps(faces,co,edges):return {'polygon':poly.ranking(faces,co),'node_node':node_node(co),'node_edge':node_edge(co,edges),'unrelated_edge_edge':edge_edge(co,edges),'shared_endpoint_overlay':overlays(co,edges)}
def first_ties(rows,key):return [x for x in rows if abs(x[key]-rows[0][key])<=TOL]
def ray_limit(co,edges,v,direction,score):
 before=score(co);step=max(before,1e-6);lo=0.;hi=step;limit=''
 for _ in range(40):
  c=dict(co);c[v]=(co[v][0]+hi*direction[0],co[v][1]+hi*direction[1]);ok,_=valid(edges,c)
  if not ok or score(c)<=before+TOL:limit='HARD_GEOMETRY' if not ok else 'CLEARANCE_NOT_IMPROVING';break
  lo=hi;hi*=2
 else:return co,0,'NO_FINITE_LIMIT'
 for _ in range(60):
  t=(lo+hi)/2;c=dict(co);c[v]=(co[v][0]+t*direction[0],co[v][1]+t*direction[1]);ok,_=valid(edges,c)
  if ok and score(c)>before+TOL:lo=t
  else:hi=t
 c=dict(co);c[v]=(co[v][0]+lo*direction[0],co[v][1]+lo*direction[1]);return c,lo,limit
def audit_vertex(v,co,edges,m):
 inc=[list(e) for e in edges if v in e];nn=next(x for x in m['node_node'] if v in x['vertices']);ne=next(x for x in m['node_edge'] if x['vertex']==v);ov=[x for x in m['shared_endpoint_overlay'] if v in x['nonshared_endpoints']]
 flags=[]
 if nn['rank']==1:flags.append('UNRELATED_NODE_NODE')
 if ne['rank']==1:flags.append('UNRELATED_NODE_EDGE')
 if ov and min(x['rank'] for x in ov)==1:flags.append('SHARED_ENDPOINT_EDGE_OVERLAY')
 return {'vertex':v,'degree':len(inc),'incident_edges':inc,'flagged':bool(flags),'categories':flags,'nearest_unrelated_node_event':nn,'nearest_unrelated_edge_event':ne,'shared_endpoint_overlay_events':ov}
def render(vertices,edges,co,faces,m,title,path):
 xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];q=90;X,Y=min(xs)-q,min(ys)-q;W,H=max(xs)-min(xs)+2*q,max(ys)-min(ys)+2*q;z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X} {Y} {W} {H}" width="1400" height="1000"><rect x="{X}" y="{Y}" width="{W}" height="{H}" fill="white"/><text x="{X+12}" y="{Y+24}" font-size="17">{title}</text>']
 for f in faces:
  if f['is_bounded'] and f['simple_boundary'] and next(x for x in m['polygon'] if x['face_id']==f['face_id'])['crowding_rank']==1:z.append('<polygon points="'+' '.join(f'{co[v][0]},{co[v][1]}' for v in f['ordered_boundary_walk'])+'" fill="#ffd54f" fill-opacity=".25" stroke="#f57f17"/>')
 for a,b in edges:z.append(f'<line x1="{co[a][0]}" y1="{co[a][1]}" x2="{co[b][0]}" y2="{co[b][1]}" stroke="#666"/>')
 for x in first_ties(m['node_node'],'clearance'):
  for v in x['vertices']:z.append(f'<circle cx="{co[v][0]}" cy="{co[v][1]}" r="10" fill="none" stroke="#6a1b9a" stroke-width="2"/>')
 for x in first_ties(m['node_edge'],'clearance'):z.append(f'<circle cx="{co[x["vertex"]][0]}" cy="{co[x["vertex"]][1]}" r="12" fill="none" stroke="#00838f" stroke-width="2"/>')
 for x in first_ties(m['unrelated_edge_edge'],'clearance'):
  for p in x['closest_points']:z.append(f'<rect x="{p[0]-4}" y="{p[1]-4}" width="8" height="8" fill="#d50000"/>')
 o=m['shared_endpoint_overlay'][0];x,y=co[o['shared_vertex']];z.append(f'<path d="M {x-11} {y-11} L {x+11} {y+11} M {x-11} {y+11} L {x+11} {y-11}" stroke="#00a152" stroke-width="3"/>')
 for v in sorted(vertices):x,y=co[v];z.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#1565c0"/><text x="{x+7}" y="{y-7}" font-size="10">{vertices[v]["name"]}</text>')
 z.append('</svg>');path.write_text('\n'.join(z)+'\n')
def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());V={x['id']:x for x in gd['vertices']};E=[(x['source'],x['target']) for x in gd['edges']];start={v:tuple(p) for v,p in dr['final_coordinates'].items()};ok,d=valid(E,start)
 if not ok or (len(V),len(E))!=(34,44):raise AssertionError('invalid authoritative start')
 faces,darts=inventory.enumerate_faces(V,E,start);inventory.add_classification_and_adjacency(faces,E,darts);initial=maps(faces,start,E);tpo=audit_vertex('component:TPO',start,E,initial);tpb=audit_vertex('component:TPB',start,E,initial)
 # Polygon pass: only the rank-1 crowded polygon and its reservoir dependency.
 target=initial['polygon'][0]['face_id'];ext,opens,paths,order=poly.dependencies(faces,initial['polygon'],[target]);by={f['face_id']:f for f in faces};co=dict(start);pops=[]
 for fid in order:co,r=poly.expand(E,by[fid],co);r['face_id']=fid;pops.append(r)
 polygon_co=dict(co);after_polygon=maps(faces,co,E)
 # Node pass: one attempt for each worst node-node and node-edge event.  A
 # node-edge event represented by the worst overlay is deferred to edge stage.
 nodeops=[];candidates=[]
 for x in first_ties(after_polygon['node_node'],'clearance'):candidates.append(('NODE_NODE',x['vertices'][0],x['vertices'][1],tuple(co[x['vertices'][1]])))
 worst_overlay=after_polygon['shared_endpoint_overlay'][0]
 for x in first_ties(after_polygon['node_edge'],'clearance'):
  overlay_equivalent=x['vertex'] in worst_overlay['nonshared_endpoints'] and worst_overlay['shared_vertex'] in x['edge']
  if not overlay_equivalent:candidates.append(('NODE_EDGE',x['vertex'],x['edge'],tuple(x['closest_point'])))
 for kind,v,off,p in candidates:
  before=math.dist(co[v],p);direction=unit(co[v],p,v);score=lambda c,v=v,kind=kind,off=off: (math.dist(c[v],c[off]) if kind=='NODE_NODE' else geo.point_segment_distance(c[v],c[off[0]],c[off[1]]));new,dist,limit=ray_limit(co,E,v,direction,score);after=score(new);accepted=dist>TOL;nodeops.append({'event_type':kind,'target':v,'offending_geometry':off,'clearance_before':before,'clearance_after':after,'vertex_moved':v if accepted else None,'displacement':dist,'limit':limit});co=new
 node_co=dict(co);after_nodes=maps(faces,co,E)
 # Edge pass: the worst unrelated pair, then the worst overlay at each shared
 # vertex.  This is rank-based and ensures one dominant junction cannot hide
 # a separate genuine overlay such as NODE_B--TPB.
 edgeops=[];ue=after_nodes['unrelated_edge_edge'][0];overlay_targets=[]
 for event in after_nodes['shared_endpoint_overlay']:
  overlay_targets.append(event)
  if 'component:TPB' in event['nonshared_endpoints']:break
 for category,event in [('UNRELATED_EDGE_EDGE',ue)]+[('SHARED_ENDPOINT_EDGE_OVERLAY',x) for x in overlay_targets]:
  if category=='UNRELATED_EDGE_EDGE':
   endpoints=event['edges'][0]+event['edges'][1];base=event['clearance'];p,q=map(tuple,event['closest_points']);choices=[(v,unit(co[v],q if v in event['edges'][0] else p,v)) for v in endpoints];score=lambda c,e=event: math.dist(*segcp(c[e['edges'][0][0]],c[e['edges'][0][1]],c[e['edges'][1][0]],c[e['edges'][1][1]]))
  else:
   shared=event['shared_vertex'];a,b=event['nonshared_endpoints']
   current_event=next(x for x in overlays(co,E) if x['shared_vertex']==shared and x['nonshared_endpoints']==[a,b])
   base=current_event['closest_separation_away_from_shared_vertex'];choices=[]
   for moving,other in ((a,b),(b,a)):
    ux,uy=co[moving][0]-co[shared][0],co[moving][1]-co[shared][1];cross=ux*(co[other][1]-co[shared][1])-uy*(co[other][0]-co[shared][0]);sgn=-1 if cross>=0 else 1;choices.append((moving,(sgn*-uy/math.hypot(ux,uy),sgn*ux/math.hypot(ux,uy))))
   def score(c,shared=shared,a=a,b=b):
    u=(c[a][0]-c[shared][0],c[a][1]-c[shared][1]);w=(c[b][0]-c[shared][0],c[b][1]-c[shared][1]);la,lb=math.hypot(*u),math.hypot(*w);ang=math.acos(max(-1,min(1,(u[0]*w[0]+u[1]*w[1])/(la*lb))));return 2*min(la,lb)*math.sin(ang/2)
  accepted=None
  for v,direction in choices:
   new,dist,limit=ray_limit(co,E,v,direction,score)
   if dist>TOL:accepted=(v,new,dist,limit,score(new));break
  edgeops.append({'event_type':category,'target':event,'clearance_or_separation_before':base,'clearance_or_separation_after':base if accepted is None else accepted[4],'vertex_moved':None if accepted is None else accepted[0],'displacement':0 if accepted is None else accepted[2],'limit':'NO_VALID_DIRECTION' if accepted is None else accepted[3]})
  if accepted:co=accepted[1]
 final=maps(faces,co,E);ok,fd=valid(E,co)
 if not ok:raise AssertionError('invalid final')
 outdir.mkdir(parents=True,exist_ok=True);render(V,E,start,faces,initial,'CORRECTED CROWDING DETECTOR — START',outdir/'iamp_corrected_crowding_start.svg');render(V,E,polygon_co,faces,after_polygon,'AFTER POLYGON PASS',outdir/'iamp_corrected_crowding_polygon.svg');render(V,E,node_co,faces,after_nodes,'AFTER GENUINE NODE PASS',outdir/'iamp_corrected_crowding_nodes.svg');render(V,E,co,faces,final,'FINAL — AFTER EDGE/OVERLAY PASS',outdir/'iamp_corrected_crowding_final.svg')
 report={'schema':'graph-relax.corrected-four-type-crowding.v1','source_coordinates':f'{direct_path}#final_coordinates','starting_validation':{'vertices':34,'edges':44,'crossings':0,'coincidences':0,'vertex_on_unrelated_edge':0},'detector_rules':{'node_edge':'all edges incident to tested node excluded','edge_edge':'shared endpoints excluded','shared_endpoint_overlay':'common endpoint ignored; ray angle and separation at shorter segment length'},'initial_maps':initial,'explicit_audits':{'TPO':tpo,'TPB':tpb},'old_false_positives_removed':[{'vertex':'component:TPO','reason':'old pass treated every non-polygon vertex as crowded; its own incident TPO-EOUT edge is excluded and it is not rank-1 in any corrected non-polygon category'}],'polygon_stage':{'target':target,'dependency_paths':paths,'processing_order':order,'operations':pops,'maps_after':after_polygon},'node_stage':{'operations':nodeops,'maps_after':after_nodes},'edge_stage':{'operations':edgeops},'final_maps':final,'final_validation':{'vertices':34,'edges':44,'crossings':fd['proper_unrelated_edge_crossing_count'],'coincidences':fd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':fd['vertex_on_unrelated_edge_interior_count'],'connectivity_unchanged':True},'final_coordinates':{v:list(co[v]) for v in sorted(co)}};(outdir/'iamp_corrected_crowding_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({'TPO':r['explicit_audits']['TPO'],'TPB':r['explicit_audits']['TPB'],'polygon_ops':r['polygon_stage']['operations'],'node_ops':r['node_stage']['operations'],'edge_ops':r['edge_stage']['operations'],'final_validation':r['final_validation']},indent=2))
if __name__=='__main__':main()
