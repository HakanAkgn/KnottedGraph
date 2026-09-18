"""Bounded, read-only audit of archived TPMS display polylines."""
import gzip
import json
from collections import Counter
from pathlib import Path
from time import perf_counter

import networkx as nx

from knotted_graph.core.embedded_contraction import certify_straight_edge_contraction

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--geometry-dir',type=Path,required=True)
_parser.add_argument('--out-dir',type=Path,required=True)
_args=_parser.parse_args()
SOURCE=_args.geometry_dir
OUTPUT=_args.out_dir
OUTPUT.mkdir(parents=True,exist_ok=True)
summaries = []
for path in sorted(SOURCE.glob('*.json.gz')):
    with gzip.open(path,'rt') as handle:
        record = json.load(handle)
    geometry = record['graph']
    edges = geometry['edges']
    summaries.append({
        'file':path.name,
        'nodes':len(geometry['nodes']),
        'edges':len(edges),
        'segments':sum(len(e['points_index'])-1 for e in edges),
        'straight_nonloop_edges':sum(e['u']!=e['v'] and len(e['points_index'])==2 for e in edges),
        'downsampled_edges':sum(e.get('original_point_count',len(e['points_index']))>len(e['points_index']) for e in edges),
    })

families = ['gyroid_to_diamond','gyroid_to_schwarz_p','schwarz_p_to_diamond']
selected = []
for family in families:
    options = [s for s in summaries if s['file'].startswith(family) and s['straight_nonloop_edges']]
    # At most four distinct archived geometries per family, spanning complexity.
    options.sort(key=lambda s:(s['segments'],s['file']))
    for index in sorted(set([0,len(options)//3,2*len(options)//3,len(options)-1])):
        if options:
            selected.append(options[index])

rows = []
for summary in selected:
    with gzip.open(SOURCE/summary['file'],'rt') as handle:
        geometry = json.load(handle)['graph']
    graph = nx.MultiGraph()
    for node in geometry['nodes']:
        graph.add_node(node['id'],pos=node['index_pos'])
    for edge in geometry['edges']:
        graph.add_edge(edge['u'],edge['v'],key=edge['key'],pts=edge['points_index'])
    candidate_edges = [(u,v,k) for u,v,k,d in graph.edges(keys=True,data=True) if u!=v and len(d['pts'])==2]
    for edge in candidate_edges[:3]:
        started = perf_counter()
        result = certify_straight_edge_contraction(graph,edge,max_segments=1500,max_predicates=100000)
        row = {'file':summary['file'],'edge':edge,'status':result.status,'reason':result.reason,
               'seconds':perf_counter()-started,'predicates':result.predicate_count}
        if result.status == 'certified':
            row['witness']=result.witness
        rows.append(row)
    print(summary['file'],dict(Counter(r['status']+':'+r['reason'] for r in rows if r['file']==summary['file'])),flush=True)
    (OUTPUT/'embedding_contraction_tpms_probe.json').write_text(json.dumps({'scope':'Archived, rounded display PL graphs only; no source-volume certificate',
        'geometry_files':len(summaries),'graphs_with_downsampled_edges':sum(s['downsampled_edges']>0 for s in summaries),
        'graphs_with_straight_nonloop_edges':sum(s['straight_nonloop_edges']>0 for s in summaries),
        'sample_selection':selected,'attempts':rows},indent=2))

(OUTPUT/'embedding_contraction_geometry_inventory.json').write_text(json.dumps(summaries,indent=2))
print('STATUS_COUNTS',dict(Counter(r['status']+':'+r['reason'] for r in rows)),flush=True)
