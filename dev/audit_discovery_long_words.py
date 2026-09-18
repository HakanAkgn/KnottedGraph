#!/usr/bin/env python3
"""Bounded retrospective long-word checks against the exported exact operators."""
import argparse,csv,hashlib,itertools,json,os,signal,time
from pathlib import Path
from datetime import datetime,timezone
import networkx as nx
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from audit_discovery_artifacts import execute_cell,sha256,write_json,git,canonical
from knotted_graph.core.embedding import ensure_embedding
from knotted_graph.projection import PDCode
from knotted_graph.invariants.yamada.polynomial import Yamada
from knotted_graph.invariants.yamada import factorized_frontier
import knotted_graph
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True);p.add_argument('--artifact',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=int,default=240);a=p.parse_args();repo=a.repo.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 assert repo in Path(knotted_graph.__file__).resolve().parents
 nbpath=repo/'User_guide/applications/05_yamada_formula_discovery.ipynb';nb=json.loads(nbpath.read_text());formula=json.loads((a.artifact/'exact_transfer_matrices.json').read_text());storedhash=formula.pop('sha256_without_hash');assert hashlib.sha256(canonical(formula).encode()).hexdigest()==storedhash
 plans=json.loads((a.artifact/'word_plans.json').read_text())['long_word_plan'];Y=sp.Symbol('Y')
 def dm(name):return DomainMatrix.from_Matrix(sp.Matrix([[sp.sympify(x,locals={'Y':Y}) for x in row] for row in formula[name]])).to_field()
 TA=dm('T_A');F=TA.domain;TB=dm('T_B').convert_to(F);L=dm('L_transpose').convert_to(F);R=dm('R').convert_to(F);I=DomainMatrix.eye((15,15),F)
 projections={}
 for letter,T in [('A',TA),('B',TB)]:
  projectors=[]
  for j,e in enumerate([2,-2,-4]):
   P=I;den=F.one
   for other in [2,-2,-4]:
    if other!=e:P=P.matmul(T.sub(I.scalarmul(F.convert(Y**other))));den*=F.convert(Y**e-Y**other)
   projectors.append((e,P.scalarmul(F.one/den)))
  projections[letter]=projectors
 run_cache={}
 def run_matrix(letter,n):
  if n==1:return TA if letter=='A' else TB
  key=letter,n
  if key not in run_cache:
   terms=[P.scalarmul(F.convert(Y**(e*n))) for e,P in projections[letter]]
   run_cache[key]=terms[0].add(terms[1]).add(terms[2])
  return run_cache[key]
 ns=dict(np=np,nx=nx,sp=sp,json=json,hashlib=hashlib,csv=csv,os=os,itertools=itertools,time=time,ensure_embedding=ensure_embedding,PDCode=PDCode,Yamada=Yamada,A=sp.Symbol('A'),CURRENT_BRANCH=git(repo,'rev-parse','--abbrev-ref','HEAD'),GIT_COMMIT=git(repo,'rev-parse','HEAD'),CONSTRUCTOR_VERSION='2026-08-27.fixed-cubic-pure-braid-word.v2')
 for i in [68,69,71,72]:execute_cell(nb,i,ns)
 provenance={'kind':'new_retrospective_long_word_audit','started_utc':datetime.now(timezone.utc).isoformat(),'formula_sha256':storedhash,'historical_freeze_verified':False,'notebook_sha256':sha256(nbpath),'script_sha256':sha256(__file__),'commit':ns['GIT_COMMIT'],'source_status':git(repo,'status','--porcelain','--','src/knotted_graph'),'source_hashes':{str(f.relative_to(repo)):sha256(f) for f in (repo/'src/knotted_graph').rglob('*.py')},'native_binary':str(factorized_frontier._yamada_factorized_frontier.__file__),'native_binary_sha256':sha256(factorized_frontier._yamada_factorized_frontier.__file__)}
 write_json(out/'provenance.json',provenance);write_json(out/'word_plan.json',plans)
 def expired(*args):raise TimeoutError('bounded audit budget expired')
 signal.signal(signal.SIGALRM,expired);signal.alarm(a.seconds);rows=[];in_progress=None
 try:
  for idx,item in enumerate(plans,1):
   print(f'PREDICT {idx}/20 length={item["length"]} {item["design"]}',flush=True);in_progress=item;start=time.time();state=L
   for letter,group in itertools.groupby(item['word']):state=state.matmul(run_matrix(letter,sum(1 for _ in group)))
   pred=state.matmul(R).to_list()[0][0];pred_time=time.time()-start
   predicted_raw={int(mon[0]-pred.denom.monoms()[0][0]):int(c) for mon,c in pred.numer.terms()} if pred else {}
   # Independent exact coefficient decoding: denominator must be a unit monomial.
   if pred:
    assert len(pred.denom.terms())==1 and pred.denom.terms()[0][1]==1
   write_json(out/f'prediction_{idx:02d}.json',{'word':item['word'],'raw_coefficients':predicted_raw,'predicted_before_this_audit_direct_evaluation':True,'historical_freeze_verified':False})
   print(f'DIRECT {idx}/20 predicted in {pred_time:.2f}s',flush=True);actual=ns['evaluate_word_raw_direct'](item['word']);coeffs={int(k):int(v) for k,v in json.loads(actual['raw_coefficients_json']).items()};passed=coeffs==predicted_raw
   rows.append({**item,'raw_coefficient_identity':passed,'prediction_seconds':pred_time,'actual':actual})
   write_json(out/'records.json',rows);assert passed,item
  status='complete'
 except TimeoutError:
  status='time_budget_exhausted'
 finally:signal.alarm(0)
 result={'status':status,'completed_cases':len(rows),'planned_cases':len(plans),'all_completed':len(rows)==len(plans),'all_completed_raw_comparisons_pass':all(r['raw_coefficient_identity'] for r in rows),'unfinished_word':in_progress if len(rows)<len(plans) else None,'historical_freeze_verified':False,'all_family_identity_proved':False}
 write_json(out/'certificate.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()

