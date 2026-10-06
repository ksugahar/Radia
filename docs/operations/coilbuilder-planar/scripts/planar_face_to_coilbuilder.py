"""Recognize constant-width planar line/circle ribbon faces without CAD guessing.

Input is an audited Netgen OCC Face with one hole, plus explicit electrical
and extrusion data. Does not execute arbitrary source files or infer turns.
"""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np

def unit(v):
    v=np.asarray(v,dtype=float);l=np.linalg.norm(v)
    if not np.isfinite(l) or l<1e-14:raise ValueError('degenerate direction')
    return v/l

def circle(points,normal):
    p0,p1,p2=np.asarray(points);a=p1-p0;b=p2-p0
    M=np.array([a,b,normal]);rhs=np.array([np.dot(a,a)/2,np.dot(b,b)/2,0.])
    if abs(np.linalg.det(M))<1e-16:raise ValueError('degenerate circular arc')
    c=p0+np.linalg.solve(M,rhs);return c,float(np.linalg.norm(p0-c))

def arc_angle(start,mid,end,center,normal):
    a=unit(start-center);b=unit(mid-center);c=unit(end-center)
    theta=lambda q:math.atan2(np.dot(normal,np.cross(a,q)),np.dot(a,q))%(2*math.pi)
    tm,te=theta(b),theta(c)
    return te if tm<te else te-2*math.pi

def _edges(face,normal,tol,scale):
    records=[]
    for wire_id,wire in enumerate(face.wires):
        for edge in wire.edges:
            lo,hi=edge.parameter_interval
            samples=np.array([list(edge.Value(float(lo+(hi-lo)*t))) for t in np.linspace(0,1,17)],dtype=float)*scale
            chord=samples[-1]-samples[0]
            isline=np.linalg.norm(chord)>tol and np.max(np.linalg.norm(np.cross(samples-samples[0],unit(chord)),axis=1))<=tol/10
            if isline:records.append({'kind':'line','wire':wire_id,'samples':samples});continue
            center,radius=circle(samples[[0,4,8]],normal)
            error=np.max(np.abs(np.linalg.norm(samples-center,axis=1)-radius))
            if error>tol/10:raise ValueError(f'unsupported boundary curve: circle-fit error {error:g} m')
            closed=np.linalg.norm(samples[-1]-samples[0])<=tol
            chunks=[samples[i:i+5] for i in (0,4,8,12)] if closed else [samples]
            for chunk in chunks:
                records.append({'kind':'arc','wire':wire_id,'samples':chunk,'center':center,'radius':radius})
    return records

def face_to_coilbuilder(face,*,current_A,turns,extrusion_vector_m,
                        length_scale_to_m=1.,section_width_m=None,
                        path_normal=None,geometry_tolerance_m=1e-8):
    """Return (CoilBuilder, certificate) for a closed planar ribbon extrusion.

    The signed current and winding normal determine polarity. The default
    normal is the extrusion direction. Geometry may be provided in CAD units
    via length_scale_to_m. The extrusion is a vector in metres, perpendicular
    to the input Face. Constant rectangular sections only. Two boundary wires
    are required; join a known symmetry half before calling this function.
    Coordinates within the declared tolerance are snapped and the measured
    adjustments are reported. Larger ambiguity or geometry changes raise.
    """
    from radia.coil_builder import CoilBuilder
    scale=float(length_scale_to_m);tol=float(geometry_tolerance_m)
    vector=np.asarray(extrusion_vector_m,dtype=float);height=float(np.linalg.norm(vector));normal=unit(vector if path_normal is None else path_normal)
    if scale<=0 or tol<=0 or not np.isfinite([scale,tol,height,current_A,turns]).all() or height<=0 or int(turns)!=turns or turns<1:raise ValueError('explicit finite SI geometry, signed current, and positive integer turns required')
    if abs(np.dot(normal,unit(vector)))<1-1e-10:raise ValueError('extrusion must be perpendicular to the coil plane')
    if len(face.wires)!=2:raise ValueError('expected one closed planar ribbon with one hole (two wires)')
    if len(face.faces)!=1:raise ValueError('expected exactly one planar face')
    face=face.faces[0]
    edges=_edges(face,normal,tol,scale)
    allpoints=np.concatenate([e['samples'] for e in edges]);plane_error=float(np.max(np.abs((allpoints-allpoints[0])@normal)))
    if plane_error>tol:raise ValueError(f'nonplanar face: {plane_error:g} m')
    used=set();primitives=[];widths=[];pairing_error=0.
    for i,a in enumerate(edges):
        if i in used:continue
        candidates=[]
        for j,b in enumerate(edges):
            if j in used or j==i or a['wire']==b['wire'] or a['kind']!=b['kind']:continue
            sa=a['samples'];sb=b['samples']
            for reverse in (False,True):
                points=sb[::-1] if reverse else sb
                if len(points)!=len(sa):continue
                delta=points-sa
                if a['kind']=='line':
                    tangent=unit(sa[-1]-sa[0]);along=np.max(np.abs(delta@tangent));constant=np.max(np.linalg.norm(delta-delta.mean(axis=0),axis=1));error=max(along,constant);width=float(np.linalg.norm(delta.mean(axis=0)))
                else:
                    width=abs(a['radius']-b['radius']);error=max(np.linalg.norm(a['center']-b['center']),np.max(np.linalg.norm((sa-a['center'])/a['radius']-(points-b['center'])/b['radius'],axis=1))*min(a['radius'],b['radius']))
                if section_width_m is not None:error=max(error,abs(width-section_width_m))
                if error<=tol and width>tol:candidates.append((float(error),j,points,width))
        if len(candidates)!=1:raise ValueError(f'edge {i} has {len(candidates)} compatible opposite edges; declare/refine geometry instead of guessing')
        error,j,points,width=candidates[0];used.update((i,j));mid=(a['samples']+points)/2+vector/2
        primitives.append({'kind':a['kind'],'samples':mid,'rail_a':a['samples']+vector/2,'rail_b':points+vector/2});widths.append(width);pairing_error=max(pairing_error,error)
    measured_width=float(np.mean(widths))
    width=measured_width if section_width_m is None else float(section_width_m)
    if max(widths)-min(widths)>2*tol:raise ValueError('variable-width face is not a constant rectangular winding pack')
    nodes=[];ends=[];snapping=0.
    for p in primitives:
        ids=[]
        for point in (p['samples'][0],p['samples'][-1]):
            matches=[k for k,n in enumerate(nodes) if np.linalg.norm(point-n)<=tol]
            if len(matches)>1:raise ValueError('ambiguous centerline junction')
            if matches:idx=matches[0];snapping=max(snapping,float(np.linalg.norm(point-nodes[idx])))
            else:idx=len(nodes);nodes.append(point.copy())
            ids.append(idx)
        ends.append(tuple(ids))
    links={k:[] for k in range(len(nodes))}
    for k,(a,b) in enumerate(ends):links[a].append(k);links[b].append(k)
    if any(len(v)!=2 for v in links.values()):raise ValueError('centerline is not a single degree-two closed chain')
    ordered=[];visited=set();start=0;at=start
    while len(visited)<len(primitives):
        choices=[k for k in links[at] if k not in visited]
        if not choices:raise ValueError('multiple disconnected centerline loops')
        k=choices[0];a,b=ends[k];rev=at==b;ordered.append((k,rev));visited.add(k);at=a if rev else b
    if at!=start:raise ValueError('centerline did not close')
    area=0.
    for k,rev in ordered:
        pts=primitives[k]['samples'];pts=pts[::-1] if rev else pts
        area+=float(np.sum(np.cross(pts[:-1],pts[1:])@normal))/2
    if area<0:ordered=[(k,not rev) for k,rev in ordered[::-1]]
    builder=CoilBuilder(current=float(current_A)*int(turns));deviation=0.;gap=0.;tangents=[];boundary_error=0.;volume=0.
    for k,rev in ordered:
        original=primitives[k]['samples'];original=original[::-1] if rev else original
        a,b=ends[k];a,b=(b,a) if rev else (a,b);entry,exit=nodes[a],nodes[b]
        if builder.segments:gap=max(gap,float(np.linalg.norm(entry-builder.segments[-1].end_pos)));entry=builder.segments[-1].end_pos.copy()
        if primitives[k]['kind']=='line':
            tangent=unit(exit-entry);frame=np.array([np.cross(tangent,normal),tangent,normal]);builder.set_start(entry,frame).set_cross_section(width,height).add_straight(float(np.linalg.norm(exit-entry)))
            converted=np.array([entry+(exit-entry)*t for t in np.linspace(0,1,len(original))]);tin=tout=tangent
            radial_vectors=np.tile(frame[0],(len(original),1));volume+=width*height*np.linalg.norm(exit-entry)
        else:
            center,radius=circle([entry,original[len(original)//2],exit],normal);theta=arc_angle(entry,original[len(original)//2],exit,center,normal);radial=unit(entry-center);sign=1 if theta>0 else -1;tangent=sign*np.cross(normal,radial)
            frame=np.array([radial,tangent,sign*normal]);builder.set_start(entry,frame).set_cross_section(width,height).add_arc(radius,float(abs(theta)*180/math.pi))
            angles=np.linspace(0,theta,len(original));converted=np.array([center+radius*(radial*np.cos(t)+np.cross(normal,radial)*np.sin(t)) for t in angles]);tin=tangent;tout=builder.segments[-1].end_orientation[1]
            radial_vectors=(converted-center)/radius;volume+=width*height*radius*abs(theta)
        deviation=max(deviation,float(np.max(np.linalg.norm(converted-original,axis=1))));tangents.append((tin,tout))
        rail_a=primitives[k]['rail_a'];rail_b=primitives[k]['rail_b']
        if rev:rail_a,rail_b=rail_a[::-1],rail_b[::-1]
        plus=converted+width/2*radial_vectors;minus=converted-width/2*radial_vectors
        direct=max(np.max(np.linalg.norm(plus-rail_a,axis=1)),np.max(np.linalg.norm(minus-rail_b,axis=1)))
        swapped=max(np.max(np.linalg.norm(plus-rail_b,axis=1)),np.max(np.linalg.norm(minus-rail_a,axis=1)))
        boundary_error=max(boundary_error,float(min(direct,swapped)))
    closure=float(np.linalg.norm(builder.segments[-1].end_pos-builder.segments[0].start_pos));tangent_error=max(float(np.linalg.norm(tangents[i][1]-tangents[(i+1)%len(tangents)][0])) for i in range(len(tangents)))
    if deviation>tol or boundary_error>tol or gap>tol or closure>1e-10:raise ValueError(f'conversion changes geometry beyond tolerance: deviation={deviation}, gap={gap}, closure={closure}')
    if tangent_error>max(1e-6,10*tol/min(width,height)):raise ValueError(f'non-tangent junctions are unsupported: {tangent_error:g}')
    input_volume=float(face.mass)*scale**2*height
    volume_error=abs(volume-input_volume)/input_volume
    if volume_error>4*tol/width:raise ValueError(f'analytic swept-volume mismatch {volume_error:g}')
    certificate={'passed':True,'boundary_sample_max_deviation_m':boundary_error,'original_extruded_volume_m3':input_volume,'generated_analytic_volume_m3':float(volume),'relative_analytic_volume_difference':float(volume_error),'schema':'radia.planar-coil-conversion.v1','method':'paired analytic line/circle boundaries; explicit extrusion and excitation','length_unit':'m','current_A':float(current_A),'turns':int(turns),'ampere_turns_A':float(current_A)*int(turns),'path_normal':normal.tolist(),'section_width_m':width,'measured_section_width_m':measured_width,'section_height_m':height,'plane_error_m':plane_error,'pairing_error_m':pairing_error,'node_snap_max_m':snapping,'centerline_sample_max_deviation_m':deviation,'join_gap_m':gap,'closure_error_m':closure,'tangent_join_error':tangent_error,'geometry_tolerance_m':tol,'segments':len(builder.segments),'source_model':'uniform current density rectangular CoilBuilder primitives; NOT automatically equal to resistive FEM current'}
    return builder,certificate

def write_builder_module(builder,path,certificate_path=None,certificate=None):
    """Save a portable build_coil(ampere_turns_A=None) module using public APIs."""
    path=Path(path);lines=['"""Generated and certified planar CoilBuilder; geometry in SI metres."""','from radia.coil_builder import CoilBuilder','','def build_coil(ampere_turns_A=None):',f'    coil = CoilBuilder(current={float(builder.current)!r} if ampere_turns_A is None else float(ampere_turns_A))']
    for seg in builder.segments:
        lines.append(f'    coil.set_start({seg.start_pos.tolist()!r}, orientation={seg.orientation.tolist()!r})')
        lines.append(f'    coil.set_cross_section({float(seg.width)!r}, {float(seg.height)!r})')
        if hasattr(seg,'radius'):lines.append(f'    coil.add_arc(radius={float(seg.radius)!r}, arc_angle={float(seg.arc_angle)!r})')
        else:lines.append(f'    coil.add_straight(length={float(seg.length)!r})')
    lines.append('    return coil');path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if certificate_path is not None:Path(certificate_path).write_text(json.dumps(certificate,indent=2),encoding='utf-8')

