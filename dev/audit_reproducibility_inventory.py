#!/usr/bin/env python3
"""Read notebook/retained figure inputs; export metadata without executing scans."""
import argparse
import ast
from collections import Counter
import csv
import hashlib
import json
import itertools
import math
from pathlib import Path
import random


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();repo=args.repo.resolve();out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    def notebook(name):
        p=repo/'User_guide/applications'/name
        return p,json.loads(p.read_text())
    dp,discovery=notebook('05_yamada_formula_discovery.ipynb')
    hp,hamiltonian=notebook('06_hamiltonian_yamada_phase_maps.ipynb')
    source=lambda nb,i: ''.join(nb['cells'][i].get('source',[]))
    controls={}
    for cell in hamiltonian['cells']:
        if cell['cell_type']!='code':continue
        try:tree=ast.parse(''.join(cell.get('source',[])))
        except SyntaxError:continue
        for node in tree.body:
            if isinstance(node,ast.Assign):
                for target in node.targets:
                    if isinstance(target,ast.Name) and target.id in {
                        'PHASE_LAMBDAS','PHASE_GAMMAS','SKELETON_DIMENSION',
                        'N_JOBS','SMOOTHING_RETRIES','PROJECTION_RETRY_SAMPLES',
                        'STABLE_MIN_COMPONENT_FRACTION'}:
                        controls[target.id]=ast.unparse(node.value)
    html=repo/'User_guide/applications/NewPhaseMapPlots/RealMaterials/html/07_hamiltonian_yamada_plotly_region_geometry_with_materials.html'
    text=html.read_text();marker='const payload = '
    payload,_=json.JSONDecoder().raw_decode(text.split(marker,1)[1])
    transitions=[]
    for tr in payload['transitions'][:5]:
        modes={}
        for name,mode in tr.get('modes',{}).items():
            modes[name]={k:v for k,v in mode.items() if isinstance(v,(str,int,float,bool)) or v is None}
        transitions.append({k:v for k,v in tr.items() if k not in {'modes','regions','colors','labels'}}|{'mode_metadata':modes})
    axes={'source':str(html.relative_to(repo)),'source_sha256':digest(html),
          'payload_version':payload.get('version'),'cache_key':payload.get('cache_key'),
          'notebook_controls':controls,'lambda_step_exact':'1/59',
          'candidate_gamma_step_exact':'99/980',
          'lambdas_exact':[{'numerator':j,'denominator':59} for j in range(60)],
          'candidate_gammas_exact':[{'numerator':294+99*j,'denominator':980} for j in range(50)],
          'embedded_lambdas_rounded':payload['lambdas'],
          'embedded_candidate_gammas_rounded':payload['candidate_gammas'],
          'analytic_transitions':transitions,
          'scope':'Retained embedded HTML data and current notebook controls; not a new phase-map rerun or certification of manuscript figure provenance.'}
    (out/'hamiltonian_sampling_inventory.json').write_text(json.dumps(axes,indent=2)+'\n')
    expected=['05a_homogeneous_theta_formula_discovery_full.csv',
              '05a_homogeneous_theta_formula_discovery_share.csv',
              '05a_homogeneous_theta_formula_discovery_heldout.csv',
              '05b_mixed_theta_transfer_words_full.csv',
              '05b_mixed_theta_count_formula_extreme_unseen.csv',
              '05c_symbolic_short_raw_training.csv','05c_symbolic_frozen_formula.json',
              '05c_symbolic_heldout_plan.csv','05c_symbolic_full_polynomial_predictions.csv',
              '05c_symbolic_heldout_actual.csv','05c_symbolic_coefficient_checks.csv',
              '05c_symbolic_heldout_summary.json']
    results=repo/'User_guide/applications/results'
    retained={p.name:str(p.relative_to(repo)) for p in results.rglob('*') if p.is_file()}
    formula_cells={str(i):{'source':source(discovery,i),'sha256':hashlib.sha256(source(discovery,i).encode()).hexdigest()} for i in [4,35,43,60,64,68,76,78,84,86]}
    mixed_ns=dict(itertools=itertools,math=math,random=random)
    names={'MAX_EXHAUSTIVE_WORD_LENGTH','STRESS_LENGTHS','STRESS_COMPOSITION_CLASSES_PER_LENGTH','STRESS_PERMUTATIONS_PER_CLASS','RNG_SEED','ALPHABET'}
    mixed_controls=ast.parse(source(discovery,43))
    mixed_controls.body=[n for n in mixed_controls.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in n.targets)]
    exec(compile(mixed_controls,'mixed-controls','exec'),mixed_ns)
    mixed_funcs=ast.parse(source(discovery,60))
    mixed_funcs.body=[n for n in mixed_funcs.body if isinstance(n,ast.FunctionDef) and n.name in {'all_words_of_length','all_count_signatures','words_for_counts','build_word_plan'}]
    exec(compile(mixed_funcs,'mixed-plan-functions','exec'),mixed_ns)
    mixed_plan=mixed_ns['build_word_plan']()
    bp=repo/'User_guide/benchmarks/results/03_knottedgraph_vs_topoly_scaling_rows.csv'
    with bp.open(newline='') as f:
        reader=csv.DictReader(f);rows=list(reader);columns=reader.fieldnames
    inventory={'discovery_notebook_sha256':digest(dp),'hamiltonian_notebook_sha256':digest(hp),
               'discovery_code_cells_with_retained_outputs':sum(bool(c.get('outputs')) for c in discovery['cells'] if c['cell_type']=='code'),
               'expected_discovery_artifacts':{n:retained.get(n) for n in expected},
               'notebook_candidate_and_plan_sources':formula_cells,
               'mixed_theta_planned_words':mixed_plan,'mixed_theta_planned_word_count':len(mixed_plan),
               'mixed_theta_plan_split_counts':dict(Counter(mixed_plan.values())),
               'historical_prompts_transcripts_or_freeze_verified':False,
               'benchmark_csv_sha256':digest(bp),'benchmark_columns':columns,
               'benchmark_row_count':len(rows),
               'benchmark_topoly_source_counts':dict(Counter(r.get('topoly_source','') for r in rows)),
               'benchmark_projection_source_counts':dict(Counter(r.get('projection_source','') for r in rows)),
               'historical_benchmark_cpu_model_recorded_in_csv':False,
               'caution':'Notebook candidate definitions and intended ordering are reproducible procedure, not evidence of historical discovery chronology.'}
    (out/'reproducibility_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps({'hamiltonian_controls':controls,'discovery_artifacts_found':sum(bool(v) for v in inventory['expected_discovery_artifacts'].values()),'benchmark_rows':len(rows)},indent=2))


if __name__=='__main__':main()
