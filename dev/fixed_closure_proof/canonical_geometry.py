"""A nonintersecting rational PL realization of the fixed cubic braid family.

The archived exterior routing is available only to reproduce its six contacts.
The new routing changes the exterior layout explicitly; no ambient equivalence
is asserted with the old, self-intersecting polygonal object.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import json


def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def zero(a):return all(x==0 for x in a)
def at(p,d,t):return tuple(x+t*y for x,y in zip(p,d))


def intersection(p,q,a,b):
    """Exact intersection: no points, one point, or endpoints of an overlap."""
    if any(max(min(p[k],q[k]),min(a[k],b[k]))>min(max(p[k],q[k]),max(a[k],b[k])) for k in range(3)):
        return ()
    d,e,o=sub(q,p),sub(b,a),sub(a,p)
    if zero(d) or zero(e):raise ValueError('Zero-length segment')
    normal=cross(d,e)
    if not zero(normal):
        if dot(o,normal)!=0:return ()
        norm=dot(normal,normal)
        t=dot(cross(o,e),normal)/norm
        s=dot(cross(o,d),normal)/norm
        return (at(p,d,t),) if 0<=t<=1 and 0<=s<=1 else ()
    if not zero(cross(o,d)):return ()
    k=next(i for i in range(3) if d[i])
    first,last=sorted(((a[k]-p[k])/d[k],(b[k]-p[k])/d[k]))
    low,high=max(F(0),first),min(F(1),last)
    if low>high:return ()
    start=at(p,d,low)
    return (start,) if low==high else (start,at(p,d,high))


def canonical_graph(word: str, *, legacy_closure: bool=False) -> dict:
    if not isinstance(word,str) or set(word)-set('AB'):
        raise ValueError('word must be a string over A,B, including the empty word')
    braid=[i for ch in word for i in ((0,0) if ch=='A' else (1,1))]
    s,h,d=F(0.90),F(0.32),F(0.45)
    lanes=(s,F(0),-s)
    n=F(max(1,len(braid)))
    def point(x,y,z=0):return (F(x),F(y),F(z))
    pos={'LT':point(0,s),'LMB':point(0,0),'LMC':point(-d,0),'LB':point(0,-s),
         'RT':point(n,s),'RMB':point(n,0),'RMC':point(n+d,0),'RB':point(n,-s)}
    paths=[[point(0,lane)] for lane in lanes]
    occupants=[0,1,2]
    for step,i in enumerate(braid):
        upper,lower=occupants[i],occupants[i+1]
        for lane,identity in enumerate(occupants):
            target=i+1 if identity==upper else i if identity==lower else lane
            height=h if identity==upper else -h if identity==lower else F(0)
            paths[identity].extend([point(F(step)+F(1,2),(lanes[lane]+lanes[target])/2,height),
                                    point(step+1,lanes[target])])
        occupants[i],occupants[i+1]=occupants[i+1],occupants[i]
    if not braid:
        for i in range(3):paths[i].append(point(n,lanes[i]))
    assert occupants==[0,1,2]
    edges=[]
    def add(u,v,points,role):edges.append({'u':u,'v':v,'points':points,'role':role})
    for u,v,path,role in zip(('LT','LMB','LB'),('RT','RMB','RB'),paths,('braid_top','braid_middle','braid_bottom')):
        add(u,v,path,role)
    for u,v,role in [('LMB','LMC','middle_left_spacer'),('RMB','RMC','middle_right_spacer'),
                     ('LT','LMB','top_middle_left'),('RT','RMB','top_middle_right'),
                     ('LMC','LB','middle_bottom_left'),('RMC','RB','middle_bottom_right')]:
        add(u,v,[pos[u],pos[v]],role)
    specs=[('RT','LT',s,3 if legacy_closure else 1,5 if legacy_closure else 3,'closure_top'),
           ('RMC','LMC',F(0),2,4,'closure_middle'),
           ('RB','LB',-s,1 if legacy_closure else 3,3 if legacy_closure else 5,'closure_bottom')]
    for u,v,y,margin,height,role in specs:
        add(u,v,[pos[u],point(n+margin,y),point(n+margin,height),
                 point(-margin,height),point(-margin,y),pos[v]],role)
    return {'word':word,'nodes':pos,'edges':edges,'legacy_closure':legacy_closure,
            'coordinates':'exact rationals','expected_braid_crossings':len(braid)}


def validate_geometry(graph: dict) -> dict:
    segments=[]
    for index,edge in enumerate(graph['edges']):
        points=edge['points']
        if points[0]!=graph['nodes'][edge['u']] or points[-1]!=graph['nodes'][edge['v']]:
            raise ValueError('Polyline endpoint mismatch')
        for j,(a,b) in enumerate(zip(points,points[1:])):
            if a==b:raise ValueError('Degenerate edge segment')
            ta=('node',edge['u']) if j==0 else ('poly',index,j)
            tb=('node',edge['v']) if j==len(points)-2 else ('poly',index,j+1)
            segments.append((index,j,a,b,{ta:a,tb:b}))
    contacts=[]
    for first,second in combinations(segments,2):
        i,j,p,q,tokens=first;k,l,a,b,other=second
        hits=intersection(p,q,a,b)
        allowed={tokens[t] for t in tokens.keys() & other.keys() if tokens[t]==other[t]}
        if hits and (len(hits)>1 or hits[0] not in allowed):
            contacts.append({'left_role':graph['edges'][i]['role'],'right_role':graph['edges'][k]['role'],
                             'points':[[str(v) for v in point] for point in hits]})
    unique_points=sorted({tuple(pt) for row in contacts for pt in row['points']})
    return {'valid_PL_embedding':not contacts,'segments':len(segments),'unintended_contacts':contacts,
            'distinct_unintended_points':[list(p) for p in unique_points],
            'coordinate_scope':'Exact supplied rational PL graph; no tolerance or sampling test.'}

if __name__=='__main__':
    from pathlib import Path
    records=[]
    for word in ('','A','B','AB','BA','AAB','ABA','AAABA','ABBAAB'):
        for legacy in (True,False):
            record={'word':word,'legacy_closure':legacy,**validate_geometry(canonical_graph(word,legacy_closure=legacy))}
            records.append(record)
            assert len(record['distinct_unintended_points'])==(6 if legacy else 0)
    out=Path(__file__).parent/'results';out.mkdir(exist_ok=True)
    (out/'geometry_audit.json').write_text(json.dumps({'source':'Archived outer-route control points; new rational braid realization',
        'new_spatial_realization_has_same_abstract_closure':True,'records':records},indent=2)+'\n')
    print('Nine old closures: six exact unintended contact points each.')
    print('Nine corrected rational PL graphs: no unintended segment intersections.')
