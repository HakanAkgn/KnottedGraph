"""Read-only stage-level TPMS extraction audit, independent of Yamada."""
import dataclasses
import hashlib
import json
from collections import Counter
from pathlib import Path
from time import perf_counter

import networkx as nx
import numpy as np

from knotted_graph.applications.phase_map_examples._tpms import tpms_families,sample_grid
from knotted_graph.applications.phase_maps import volume_topology
from knotted_graph.core import remove_leaf_nodes,simplify_edges,smooth_edges
from knotted_graph.extraction import skeletonize_volume
from knotted_graph.extraction._optimized import sparse_adjacency_exact_cropped
from knotted_graph.extraction._topology_optimized import _prepared_components,_trace_prepared,_diagnostic_summary,persistent_extract

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--out-dir',type=Path,required=True)
OUT=_parser.parse_args().out_dir
OUT.mkdir(parents=True,exist_ok=True)

def summary(g):
    components=list(nx.connected_components(g))
    return {'V':len(g),'E':g.number_of_edges(),'b0':len(components),
            'b1':g.number_of_edges()-len(g)+len(components),
            'genera':sorted(g.subgraph(c).number_of_edges()-len(c)+1 for c in components),
            'max_degree':max((d for _,d in g.degree()),default=0)}

def betti(mask):
    t=volume_topology(mask)
    return {'b0':t.connected_components,'b1':t.handle_rank,'b2':t.enclosed_voids,'voxels':t.interior_voxels,'euler':t.euler_characteristic}

def fingerprint(g):
    nodes=[(int(n),tuple(d['pos'])) for n,d in g.nodes(data=True)]
    edges=[(int(a),int(b),int(k),np.asarray(d['pts']).tolist()) for a,b,k,d in g.edges(keys=True,data=True)]
    return repr((nodes,edges))

rows=[]
source={}
import knotted_graph.extraction._tracing as tracing
import knotted_graph.extraction._topology_optimized as topology
import knotted_graph.core.embedding as embedding
for module in (tracing,topology,embedding):
    path=Path(module.__file__)
    source[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()

for dimension in (24,64,96):
    for family in tpms_families(dimension=dimension,thresholds=(0.,.105,.3)):
        _,x,y,z=sample_grid(family.span,dimension)
        domain=family.domain.sample(x,y,z)
        for lam in (0.,.25,.5,.75,1.):
            field=family.field_at(lam,x,y,z)
            for c in family.thresholds:
                started=perf_counter()
                mask=(field<=c)&(domain<=0)
                row={'dimension':dimension,'family':family.key,'lambda':lam,'c':c,'mask':betti(mask)}
                try:
                    skeleton=skeletonize_volume(mask)
                    row['skeleton']=betti(skeleton)
                    coords,adj=sparse_adjacency_exact_cropped(skeleton)
                    prepared=_prepared_components(coords,adj)
                    row['voxel_adjacency']={'V':len(coords),'E':sum(map(len,adj))//2,'b0':len(prepared),
                        'b1':sum(map(len,adj))//2-len(coords)+len(prepared)}
                    candidates=[]
                    graphs=[]
                    for hop in range(5):
                        graph=_trace_prepared(prepared,hop)
                        graphs.append(graph)
                        reduced,clean,fp,safe=_diagnostic_summary(graph,max_degree=None,anomaly_ratio=.15)
                        candidates.append({'hop':hop,**summary(graph),'clean':bool(clean),'one_hop_safe':bool(safe)})
                    row['candidates']=candidates
                    selected=persistent_extract(coords,adj)
                    row['selected']={**summary(selected),'matching_hops':[i for i,g in enumerate(graphs) if fingerprint(g)==fingerprint(selected)]}
                    leaf=remove_leaf_nodes(selected)
                    row['leaf_pruned']=summary(leaf)
                    try:
                        simplified=simplify_edges(leaf)
                        row['simplified']=summary(simplified)
                        smoothed=smooth_edges(simplified,epsilon=0)
                        row['smoothed']=summary(smoothed)
                    except Exception as exc:
                        row['postprocessing_error']=type(exc).__name__+':'+str(exc)
                except Exception as exc:
                    row['extraction_error']=type(exc).__name__+':'+str(exc)
                row['seconds']=perf_counter()-started
                rows.append(row)
        (OUT/'tpms_extraction_stage_diagnosis.json').write_text(json.dumps({'source_sha256':source,'records':rows},indent=2))
        completed=[r for r in rows if r['dimension']==dimension and r['family']==family.key]
        print(dimension,family.key,'rows',len(completed),'selected_bad_betti',sum('selected'in r and any(r['selected'][b]!=r['mask'][b] for b in ('b0','b1')) for r in completed),
              'post_errors',sum('postprocessing_error'in r for r in completed),flush=True)

print('DONE',len(rows),flush=True)
