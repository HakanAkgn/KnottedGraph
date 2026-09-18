#!/usr/bin/env python3
"""Regression-test notebook06 functions and benchmark a bounded volume pilot."""
import argparse,ast,dataclasses,hashlib,json,time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import networkx as nx
import sympy as sp
from skimage.morphology import skeletonize
from knotted_graph.core import remove_leaf_nodes,simplify_edges,smooth_edges
from knotted_graph.extraction import skeleton_image_to_graph
from knotted_graph.applications.nodal.skeleton import NodalSkeleton
from knotted_graph.applications.nodal.deformation import NodalBlochPath
from knotted_graph.applications.nodal.models import hopf_link_bloch_vector,trefoil_bloch_vector,solomon_bloch_vector,unknot_bloch_vector,pq_torus_knot_bloch_vector
from knotted_graph.applications.phase_maps import volume_topology
def namespace(notebook):
 nb=json.loads(Path(notebook).read_text())
 ns=dict(dataclasses=dataclasses,np=np,nx=nx,sp=sp,NodalSkeleton=NodalSkeleton,remove_leaf_nodes=remove_leaf_nodes,simplify_edges=simplify_edges,smooth_edges=smooth_edges,skeleton_image_to_graph=skeleton_image_to_graph,A=sp.Symbol('A'),Counter=Counter,ALLOW_WINDOW_BOUNDARY_CONTACT=False,APPLY_DISPLAY_FILTER=False,ALLOW_ABSTRACT_CONTRACTION_DIAGNOSTIC=False,SMOOTHING_RETRIES=(4.,1.,.25,0.),PROJECTION_RETRY_SAMPLES=(1,),PHASE_LAMBDAS=np.array([0.,1.]),PHASE_GAMMAS=np.array([.3,.4]),STABLE_MIN_COMPONENT_CELLS=5,SKELETON_DIMENSION=64)
 exec(compile(''.join(nb['cells'][7]['source']),'notebook06-cell7','exec'),ns)
 keep={'record_to_json','record_from_json','gammas_from_records','records_to_grid','phase_ids_for_records','connected_label_components','label_component_stats','stable_partition_labels','stable_phase_ids_for_records','is_one_vertex_yamada_record','gamma_row_is_all_vertex','terminal_transition_gammas_from_records'}
 tree=ast.parse(''.join(nb['cells'][8]['source']));tree.body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in keep];exec(compile(tree,'notebook06-cell8-functions','exec'),ns)
 ns['stable_min_component_cells']=lambda *args:5
 return ns
def tests(ns):
 passed=[]
 def check(name,fn):fn();passed.append(name)
 def fake(mask):return SimpleNamespace(_interior_mask=mask,_skeleton_image=skeletonize(mask))
 def reject(mask):
  try:ns['core_candidates_from_skeleton'](fake(mask))
  except ns['UnsupportedPhaseVolume']:return
  raise AssertionError('unsupported mask accepted')
 empty=np.zeros((20,20,20),bool)
 check('empty volume stays unsupported',lambda:reject(empty))
 boundary=empty.copy();boundary[0:6,5:10,5:10]=True
 check('boundary contact stays unsupported',lambda:reject(boundary))
 shell=empty.copy();shell[3:17,3:17,3:17]=True;shell[7:13,7:13,7:13]=False
 check('positive b2 rejected',lambda:reject(shell))
 topo=ns['interior_topology_summary'](fake(shell));assert topo['enclosed_voids']==1 and topo['handle_rank']==0;passed.append('b2 included in Euler/Betti relation')
 wrong=nx.MultiGraph();wrong.add_node(0)
 disconnected=nx.MultiGraph();disconnected.add_node(0,pos=np.array([0.,0.,0.]));disconnected.add_node(1,pos=np.array([3.,0.,0.]));disconnected.add_edge(0,0,pts=np.array([[0.,0.,0.],[1.,1.,0.],[0.,0.,0.]]))
 cleaned=ns['drop_degenerate_edges'](disconnected);assert 1 in cleaned and nx.number_connected_components(cleaned)==2;passed.append('cleanup preserves valid isolated components')
 try:ns['_assert_graph_topology'](wrong,dict(components=1,handle_rank=1),'test')
 except ValueError:passed.append('graph cycle loss rejected')
 else:raise AssertionError('cycle mismatch accepted')
 def fail(*args,**kwargs):raise ValueError('projection failed')
 original=ns['_compute_yamada_audited'];ns['_compute_yamada_audited']=fail
 try:ns['evaluate_yamada_on_candidates']([('test',wrong)])
 except ValueError as exc:assert 'unavailable' in str(exc);passed.append('projection failure never abstract/point fallback')
 else:raise AssertionError('failed evaluation accepted')
 high=nx.MultiGraph();high.add_node(0);high.add_edge(0,0);high.add_edge(0,0)
 ns['_compute_yamada_audited']=lambda *args:(sp.Integer(7),dict(evaluation_kind='diagram-yamada',projection={'pd_code':'recorded'}))
 source,_,graph,_=ns['evaluate_yamada_on_candidates']([('test',high)]);assert source=='diagram-yamada' and graph.graph['yamada_evaluation_audit']['projection'];passed.append('high-valence diagram metadata retained')
 ns['_compute_yamada_audited']=original
 R=ns['HamiltonianPhaseRecord']
 def record(lam,gamma,signature,error=None):
  return R('t','t',lam,gamma,'unavailable' if error else 'spatial-yamada','test',1,0,1,0,(),signature,'' if error else '-1',error)
 bad=record(0,.3,'unavailable','failure')
 try:ns['phase_ids_for_records']([bad],lambdas=[0],gammas=[.3])
 except ValueError:passed.append('unavailable cells excluded from phase coloring')
 else:raise AssertionError('unavailable became phase')
 records=[record(x,g,'a' if i else 'b') for i,(g,x) in enumerate((g,x) for g in [.3,.4] for x in [0,1])]
 labels,*_=ns['phase_ids_for_records'](records,lambdas=[0,1],gammas=[.3,.4])
 res=ns['stable_phase_ids_for_records'](records,lambdas=[0,1],gammas=[.3,.4]);assert np.array_equal(res[2],labels) and not res[3].any();passed.append('raw labels unchanged with display filter off')
 assert not ns['is_one_vertex_yamada_record'](bad);passed.append('failure cannot trigger terminal vertex row')
 gammas,index=ns['terminal_transition_gammas_from_records'](records,lambdas=[0,1],gammas=[.3,.4]);assert index is None and np.array_equal(gammas,[.3,.4]);passed.append('no validated terminal preserves full requested energy range')
 return passed
def benchmark(ns):
 endpoints={'hopf_to_trefoil':(hopf_link_bloch_vector,trefoil_bloch_vector),'hopf_to_solomon':(hopf_link_bloch_vector,solomon_bloch_vector),'unknot_to_solomon':(unknot_bloch_vector,solomon_bloch_vector),'trefoil_to_cinquefoil':(trefoil_bloch_vector,lambda g:pq_torus_knot_bloch_vector(2,5,g))}
 cells=[('hopf_to_trefoil',0.,.3),('hopf_to_trefoil',.5,.8),('hopf_to_solomon',.7,1.1),('unknot_to_solomon',.5,1.3),('trefoil_to_cinquefoil',.25,.3)]
 rows=[]
 for key,lam,gamma in cells:
  t=time.perf_counter();s=NodalSkeleton(char=NodalBlochPath(*endpoints[key]).at(gamma,lam),dimension=64);topo=ns['interior_topology_summary'](s)
  try:
   candidates=ns['core_candidates_from_skeleton'](s);g=candidates[-1][1];status='betti-consistent-extraction';graphb=(nx.number_connected_components(g),g.number_of_edges()-g.number_of_nodes()+nx.number_connected_components(g));error=None
  except Exception as exc:status='unsupported' if isinstance(exc,ns['UnsupportedPhaseVolume']) else 'unavailable';graphb=None;error=str(exc)
  row=dict(transition=key,lam=lam,gamma=gamma,dimension=64,seconds=time.perf_counter()-t,status=status,volume=topo,graph_betti=graphb,error=error);rows.append(row);print(json.dumps(row),flush=True)
 return rows
def main():
 p=argparse.ArgumentParser();p.add_argument('--notebook',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);ns=namespace(a.notebook);passed=tests(ns);(a.out/'regression_checks.json').write_text(json.dumps({'passed':passed,'count':len(passed),'notebook_sha256':hashlib.sha256(a.notebook.read_bytes()).hexdigest()},indent=2)+'\n');print('TESTS PASS',len(passed),flush=True);rows=benchmark(ns);(a.out/'five_cell_timing.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
