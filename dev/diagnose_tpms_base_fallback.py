"""Examine zero-radius fallback only for empirically incorrect persistence picks."""
import json
from pathlib import Path
import networkx as nx
import numpy as np
from knotted_graph.applications.phase_map_examples._tpms import tpms_families,sample_grid
from knotted_graph.core import remove_leaf_nodes,simplify_edges,smooth_edges
from knotted_graph.extraction import skeletonize_volume,skeleton_image_to_graph
from knotted_graph.extraction._topology_optimized import _embedded_geometry_safe

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--out-dir',type=Path,required=True)
OUT=_parser.parse_args().out_dir
OUT.mkdir(parents=True,exist_ok=True)
rows=json.loads((OUT/'tpms_extraction_stage_diagnosis.json').read_text())['records']
selected=[r for r in rows if 'selected'in r and any(r['mask'][k]!=r['selected'][k] for k in ('b0','b1'))]
result=[]
for r in selected:
    family=next(f for f in tpms_families(dimension=r['dimension'],thresholds=(r['c'],)) if f.key==r['family'])
    _,x,y,z=sample_grid(family.span,r['dimension'])
    mask=(family.field_at(r['lambda'],x,y,z)<=r['c'])&(family.domain.sample(x,y,z)<=0)
    graph=skeleton_image_to_graph(skeletonize_volume(mask),adaptive_max_hops=0)
    row={k:r[k] for k in ('dimension','family','lambda','c','mask')}
    row['base_geometry_safe']=bool(_embedded_geometry_safe(graph))
    row['base_max_degree']=max((d for _,d in graph.degree()),default=0)
    try:
        graph=remove_leaf_nodes(graph)
        if graph.number_of_edges():
            graph=simplify_edges(graph)
            graph=smooth_edges(graph,epsilon=0)
        row['result']={'V':len(graph),'E':graph.number_of_edges(),'b0':nx.number_connected_components(graph),
                       'b1':graph.number_of_edges()-len(graph)+nx.number_connected_components(graph)}
        row['matches']=all(row['result'][k]==r['mask'][k] for k in ('b0','b1'))
    except Exception as exc:
        row['error']=type(exc).__name__+':'+str(exc)
    result.append(row)
(OUT/'tpms_zero_radius_fallback_diagnosis.json').write_text(json.dumps(result,indent=2))
print('cases',len(result),'safe_geometry',sum(r['base_geometry_safe'] for r in result),'passing_betti',sum(r.get('matches',False) for r in result),'errors',sum('error'in r for r in result))
print('base_max_degrees',sorted(set(r['base_max_degree'] for r in result)))
