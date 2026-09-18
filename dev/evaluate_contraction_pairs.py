"""Small bounded TPMS comparison probe; exact descendants only, no conjectural merges."""
import gzip
import json
from collections import Counter
from pathlib import Path
from time import perf_counter

import networkx as nx

from knotted_graph.core.embedded_contraction import compare_by_certified_contractions, embedded_fingerprint

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--tpms-dir',type=Path,required=True)
_parser.add_argument('--out-dir',type=Path,required=True)
_args=_parser.parse_args()
ROOT=_args.tpms_dir
OUTPUT=_args.out_dir
OUTPUT.mkdir(parents=True,exist_ok=True)
source=json.loads((ROOT/'data/tpms_compact_c0_03_stable_up_to_contraction_phase_row_panel_c_evolution_source_data.json').read_text())

def load_graph(name):
    with gzip.open(ROOT/'geometry'/name,'rt') as handle:
        geometry=json.load(handle)['graph']
    graph=nx.MultiGraph()
    for node in geometry['nodes']:
        graph.add_node(node['id'],pos=node['index_pos'])
    for edge in geometry['edges']:
        graph.add_edge(edge['u'],edge['v'],key=edge['key'],pts=edge['points_index'])
    return graph

rows=[]
for family in source['transition_order']:
    grid=source['grids_original_phase_ids'][family]
    candidates=[col for col in range(20) if grid[7][col]==grid[7][col+1]]
    chosen=sorted(set([candidates[0],candidates[len(candidates)//2],candidates[-1]]))
    for col in chosen:
        a=f'{family}_lambda{col:03d}_c007.json.gz'
        b=f'{family}_lambda{col+1:03d}_c007.json.gz'
        left,right=load_graph(a),load_graph(b)
        started=perf_counter()
        result=compare_by_certified_contractions(left,right,max_depth=1,max_states=8,max_attempts=16,max_predicates=100000)
        row={'left':a,'right':b,'old_shared_class':grid[7][col],
             'abstract_isomorphic':nx.is_isomorphic(left,right),
             'status':result.status,'reason':result.reason,
             'states':result.explored_states,'seconds':perf_counter()-started,
             'witness':result.witness}
        rows.append(row)
        print(a,b,result.status,result.reason,flush=True)
        (OUTPUT/'embedding_contraction_pair_probe.json').write_text(json.dumps({'scope':'Nine neighboring pairs sharing an old stabilized contraction label; archived rounded PL geometry only',
            'limits':{'max_depth':1,'max_states':8,'max_attempts':16,'max_predicates_per_move':100000},'results':rows},indent=2))
print(dict(Counter(row['status']+':'+row['reason'] for row in rows)),flush=True)
