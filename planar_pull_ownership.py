"""Generic face-incidence and exclusive Pull-Cell ownership primitives.

Face incidence is many-to-many planar topology.  Pull-Cell ownership is
movement authority and is zero-or-one per graph vertex.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping, Sequence


def boundary_edges_of_walk(walk: Sequence[str]) -> set[tuple[str,str]]:
    return {tuple(sorted((walk[i],walk[(i+1)%len(walk)]))) for i in range(len(walk))}


def derive_face_incidence_and_free_for_face(
        face_boundary_walks: Mapping[str,Sequence[str]],
        graph_edges: Iterable[tuple[str,str]]) -> tuple[dict[str,list[str]],dict[str,list[str]],dict[str,dict]]:
    edges=sorted(tuple(sorted(edge)) for edge in graph_edges)
    incident=defaultdict(set)
    for edge in edges:
        for vertex in edge:incident[vertex].add(edge)
    boundary_faces=defaultdict(list);free_faces=defaultdict(list);details={}
    for face_id in sorted(face_boundary_walks):
        walk=list(face_boundary_walks[face_id]);face_edges=boundary_edges_of_walk(walk)
        for vertex in dict.fromkeys(walk):
            boundary_faces[vertex].append(face_id)
            local={edge for edge in face_edges if vertex in edge}
            qualifies=len(incident[vertex])==2 and len(local)==2 and incident[vertex]==local
            if qualifies:free_faces[vertex].append(face_id)
            details[f"{vertex}|{face_id}"]={"vertex":vertex,"face_id":face_id,"boundary_incident":True,
                "free_for_face":qualifies,"degree":len(incident[vertex]),"incident_graph_edges":[list(e) for e in sorted(incident[vertex])],
                "incident_face_boundary_edges":[list(e) for e in sorted(local)]}
    return ({v:sorted(fs) for v,fs in sorted(boundary_faces.items())},
            {v:sorted(fs) for v,fs in sorted(free_faces.items())},details)


def assign_pull_ownership(boundary_faces_by_vertex,free_for_faces_by_vertex,pull_layer_by_face):
    owner={};cells={face:[] for face in sorted(pull_layer_by_face)}
    table=[]
    for vertex in sorted(boundary_faces_by_vertex):
        qualifying=sorted(free_for_faces_by_vertex.get(vertex,[]),key=lambda face:(pull_layer_by_face[face],face))
        chosen=qualifying[0] if qualifying else None
        if chosen:
            owner[vertex]=chosen;cells[chosen].append(vertex)
        table.append({"vertex":vertex,"all_boundary_incidental_selected_faces":boundary_faces_by_vertex[vertex],
                      "all_free_for_face_qualifying_faces":sorted(free_for_faces_by_vertex.get(vertex,[])),
                      "owner_face":chosen,"pull_layer":pull_layer_by_face[chosen] if chosen else None,
                      "pull_cell":f"cell:{chosen.split(':',1)[1]}" if chosen else None})
    members=[vertex for values in cells.values() for vertex in values]
    if len(members)!=len(set(members)):raise AssertionError("duplicate Pull Cell membership")
    if set(owner)!=set(free_for_faces_by_vertex):raise AssertionError("movable vertex lacks exactly one Owner Face")
    return owner,{face:sorted(vertices) for face,vertices in cells.items()},table


def calculate_owned_destinations(coordinates,pole,owner_face_by_vertex,pull_factor_by_face):
    destinations={};factor_by_vertex={}
    for vertex in sorted(owner_face_by_vertex):
        face=owner_face_by_vertex[vertex];factor=pull_factor_by_face[face];x,y=coordinates[vertex]
        if vertex in destinations:raise AssertionError("multiple destination coordinates")
        factor_by_vertex[vertex]=factor
        destinations[vertex]=(pole[0]+factor*(x-pole[0]),pole[1]+factor*(y-pole[1]))
    if set(destinations)!=set(owner_face_by_vertex) or set(factor_by_vertex)!=set(owner_face_by_vertex):raise AssertionError("incomplete owned movement")
    return destinations,factor_by_vertex


def plant_destinations(coordinates,destinations):
    planted=dict(coordinates);planted.update(destinations);return planted


def rebuild_edge_segments(coordinates,graph_edges):
    return [{"edge":[a,b],"start":list(coordinates[a]),"end":list(coordinates[b])} for a,b in graph_edges]


def assert_reconnected(segments,coordinates,graph_edges):
    expected=[tuple(edge) for edge in graph_edges]
    if [tuple(segment["edge"]) for segment in segments]!=expected:raise AssertionError("electrical endpoints changed")
    for segment in segments:
        a,b=segment["edge"]
        if segment["start"]!=list(coordinates[a]) or segment["end"]!=list(coordinates[b]):raise AssertionError("stale edge endpoint")

