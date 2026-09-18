#!/usr/bin/env python3
"""Fresh finite-window homology/extraction pilot; no Yamada phase classification."""
import argparse,concurrent.futures,hashlib,json,multiprocessing,os,time
from datetime import datetime,timezone
from pathlib import Path
import networkx as nx
from audit_hamiltonian_notebook import namespace
from knotted_graph.applications.nodal.skeleton import NodalSkeleton
from knotted_graph.applications.nodal.deformation import NodalBlochPath
from knotted_graph.applications.nodal.models import hopf_link_bloch_vector,trefoil_bloch_vector,solomon_bloch_vector,unknot_bloch_vector,pq_torus_knot_bloch_vector
_NS=None
_PATHS=None
def cinquefoil(g):return pq_torus_knot_bloch_vector(2,5,g)
def initialize(notebook):
 global _NS,_PATHS
 _NS=namespace(notebook)
 _PATHS={key:NodalBlochPath(start,end) for key,start,end in [
 ('hopf_to_trefoil',hopf_link_bloch_vector,trefoil_bloch_vector),
 ('hopf_to_solomon',hopf_link_bloch_vector,solomon_bloch_vector),
 ('unknot_to_trefoil',unknot_bloch_vector,trefoil_bloch_vector),
 ('unknot_to_solomon',unknot_bloch_vector,solomon_bloch_vector),
 ('trefoil_to_cinquefoil',trefoil_bloch_vector,cinquefoil)]}
def evaluate(item):
 t=time.perf_counter();key,li,gi=item;lam=li/10;gamma=.3+99*gi/980;topology=None
 try:
  s=NodalSkeleton(char=_PATHS[key].at(gamma,lam),dimension=64)
  topology=_NS['interior_topology_summary'](s)
  candidates=_NS['core_candidates_from_skeleton'](s)
  graph=candidates[-1][1]
  components=nx.number_connected_components(graph);cycle=graph.number_of_edges()-graph.number_of_nodes()+components
  status='betti-consistent-extraction';error=None
  stats=dict(nodes=graph.number_of_nodes(),edges=graph.number_of_edges(),components=components,cycle_rank=cycle,max_degree=max(dict(graph.degree()).values(),default=0),candidate_mode=candidates[-1][0])
 except Exception as exc:
  error=f'{type(exc).__name__}: {exc}';stats=None
  if isinstance(exc,_NS['UnsupportedPhaseVolume']):
   status='empty' if topology and not topology['interior_voxels'] else 'unsupported-cavity' if topology and topology['enclosed_voids'] else 'unsupported-boundary'
  else:status='unavailable-extraction'
 return dict(transition=key,lambda_index=li,gamma_index=gi,lam=lam,gamma=gamma,dimension=64,volume=topology,graph=stats,status=status,error=error,seconds=time.perf_counter()-t)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--notebook',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=int,default=600);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 specs=[('hopf_to_trefoil',16),('hopf_to_solomon',24),('unknot_to_trefoil',16),('unknot_to_solomon',24),('trefoil_to_cinquefoil',39)]
 plan=[(key,li,gi) for key,ng in specs for gi in range(ng) for li in range(11)]
 def save(name,data):
  temp=a.out/(name+'.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(a.out/name)
 save('plan.json',{'dimensions':64,'lambda_array':[i/10 for i in range(11)],'candidate_gammas':[.3+99*i/980 for i in range(50)],'families':dict(specs),'planned_cells':len(plan),'cells':plan})
 save('provenance.json',{'run_kind':'new_coarse_finite_window_volume_extraction_pilot','started_utc':datetime.now(timezone.utc).isoformat(),'notebook_sha256':hashlib.sha256(a.notebook.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_hashes':{str(f.relative_to(a.repo)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.repo/'src/knotted_graph').rglob('*.py')},'dimension':64,'workers':2,'yamada_evaluated':False,'display_filter':False,'convergence_established':False,'regular_neighborhood_or_isotopy_certified':False})
 start=time.monotonic();rows=[];status='complete'
 pool=concurrent.futures.ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(str(a.notebook.resolve()),))
 futures=[pool.submit(evaluate,item) for item in plan]
 try:
  for f in concurrent.futures.as_completed(futures,timeout=a.seconds):
   rows.append(f.result())
   if len(rows)%25==0 or len(rows)==len(plan):
    save('records.json',rows);print(f'PILOT {len(rows)}/{len(plan)} elapsed={time.monotonic()-start:.1f}s',flush=True)
 except concurrent.futures.TimeoutError:
  status='time-budget-exhausted'
  for f in futures:f.cancel()
 finally:
  if status!='complete':
   for proc in pool._processes.values():proc.terminate()
  pool.shutdown(wait=True,cancel_futures=True)
  save('records.json',rows)
 from collections import Counter
 summary={'status':status,'completed_cells':len(rows),'planned_cells':len(plan),'seconds':time.monotonic()-start,'status_counts':dict(Counter(r['status'] for r in rows)),'by_transition':{key:dict(Counter(r['status'] for r in rows if r['transition']==key)) for key,_ in specs},'yamada_evaluated':False,'convergence_established':False,'historical_7140_rerun':False}
 save('summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()

