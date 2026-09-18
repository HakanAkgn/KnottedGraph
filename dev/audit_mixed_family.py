#!/usr/bin/env python3
"""Retrospective direct check of the notebook's mixed-theta word plan."""
import argparse,ast,hashlib,json,time
from pathlib import Path
from datetime import datetime,timezone
from audit_discovery_artifacts import execute_cell,sha256,write_json,git
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();repo=a.repo.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 nbpath=repo/'User_guide/applications/05_yamada_formula_discovery.ipynb';nb=json.loads(nbpath.read_text())
 ns={};tree=ast.parse(''.join(nb['cells'][41]['source']));tree.body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))];exec(compile(tree,'setup-imports','exec'),ns)
 assert repo in Path(ns['knotted_graph'].__file__).resolve().parents
 ns.update(ROOT=repo,RESULTS_DIR=out,PROJECTION_CACHE_DIR=out/'projection_cache',FULL_CSV=out/'records.csv',SHARE_CSV=out/'share.csv',A=ns['sp'].Symbol('A'),CURRENT_BRANCH=git(repo,'rev-parse','--abbrev-ref','HEAD'),GIT_COMMIT=git(repo,'rev-parse','HEAD'),SOURCE_CONSTRUCTOR_VERSION='2026-08-25.three-subcubic-theta-families.safe-return.v3',CONSTRUCTOR_VERSION='2026-08-27.mixed-theta-word-transfer-test.v1',RESULT_SCHEMA_VERSION=2)
 execute_cell(nb,43,ns);ns['REUSE_PROJECTION_CACHE']=False;ns['WRITE_PROJECTION_CACHE']=False
 for i in [47,48,49,50,52,54,55,56,58,60]:execute_cell(nb,i,ns)
 provenance={'kind':'new_retrospective_exact_audit','started_utc':datetime.now(timezone.utc).isoformat(),'commit':ns['GIT_COMMIT'],'source_status':git(repo,'status','--porcelain','--','src/knotted_graph'),'notebook_sha256':sha256(nbpath),'script_sha256':sha256(__file__),'selected_cells':[43,47,48,49,50,52,54,55,56,58,60],'source_hashes':{str(f.relative_to(repo)):sha256(f) for f in (repo/'src/knotted_graph').rglob('*.py')},'historical_freeze_verified':False}
 write_json(out/'provenance.json',provenance)
 plan=ns['build_word_plan']();write_json(out/'word_plan.json',plan)
 ns['run_mixed_preflights']()
 rows=[];started=time.time()
 for i,(word,split) in enumerate(plan.items(),1):
  try:
   row=ns['evaluate_mixed_word'](word,split=split)
  except Exception as exc:
   row={'word':word,'status':'error','error':f'{type(exc).__name__}: {exc}'}
  rows.append(row)
  if i%10==0 or i==len(plan):
   write_json(out/'records.json',rows);print(f'MIXED {i}/{len(plan)} elapsed={time.time()-started:.1f}s failures={sum(not r.get("master_formula_pass") for r in rows)}',flush=True)
 result={'completed_cases':len(rows),'planned_cases':len(plan),'all_completed':len(rows)==len(plan),'all_coefficient_comparisons_pass':all(r.get('master_formula_pass') for r in rows),'failures':[r['word'] for r in rows if not r.get('master_formula_pass')],'historical_freeze_verified':False,'all_family_law_proved':False}
 write_json(out/'certificate.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()

