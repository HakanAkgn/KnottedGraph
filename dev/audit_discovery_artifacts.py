#!/usr/bin/env python3
"""Export new, explicitly dated discovery-audit evidence from notebook constructors.

Run with the repository's uv environment. This does not recreate historical LLM
conversations or establish that a historical candidate preceded its holdout.
All outputs are written to --out; the inspected checkout is read only.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import itertools
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import networkx as nx
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def git(repo, *args):
    p = subprocess.run(["git", "-C", str(repo), *args], text=True, capture_output=True)
    return p.stdout.strip() if p.returncode == 0 else None


def cell_source(notebook, index):
    return "".join(notebook["cells"][index].get("source", []))


def execute_cell(notebook, index, ns, *, definitions_only=False):
    source = cell_source(notebook, index)
    tree = ast.parse(source)
    if definitions_only:
        tree.body = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    exec(compile(tree, f"notebook-cell-{index}", "exec"), ns)


def payload(matrix):
    return [[str(x) for x in row] for row in matrix.to_list()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--symbolic", action="store_true", help="reconstruct/check exact symbolic transfer operators")
    parser.add_argument("--heldout", type=int, default=0, help="number of deterministic long words to evaluate after new freeze (0..20)")
    args = parser.parse_args()
    repo, out = args.repo.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if not 0 <= args.heldout <= 20:
        parser.error("--heldout must be between zero and twenty")
    notebook_path = repo / "User_guide/applications/05_yamada_formula_discovery.ipynb"
    notebook = json.loads(notebook_path.read_text())
    import knotted_graph
    from knotted_graph.core.embedding import ensure_embedding
    from knotted_graph.projection import PDCode
    from knotted_graph.invariants.yamada.polynomial import Yamada
    from knotted_graph.invariants.yamada.factorized_frontier import native_factorized_available

    if repo not in Path(knotted_graph.__file__).resolve().parents:
        raise RuntimeError("Imported package is not the requested checkout")
    provenance = {
        "run_kind": "new_reproduction_audit_not_historical_discovery",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "repo": str(repo), "commit": git(repo, "rev-parse", "HEAD"),
        "branch": git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        "source_status": git(repo, "status", "--porcelain", "--", "src/knotted_graph"),
        "notebook_sha256": sha256(notebook_path),
        "audit_script_sha256": sha256(__file__),
        "notebook_cells_used": [68, 69, 71, 72, 86],
        "notebook_cell_hashes": {str(i): hashlib.sha256(cell_source(notebook,i).encode()).hexdigest() for i in [68,69,71,72,86]},
        "source_file_sha256": {str(p.relative_to(repo)): sha256(p) for p in sorted((repo/'src/knotted_graph').rglob('*')) if p.suffix in {'.py','.cpp','.hpp','.c','.h'}},
        "python": sys.version, "platform": platform.platform(),
        "machine": platform.machine(), "processor": platform.processor(),
        "packages": {p: importlib.metadata.version(p) for p in ['sympy','numpy','networkx']},
        "native_factorized_available": native_factorized_available(),
        "historical_prompts_transcripts_found": False,
        "historical_freeze_verified": False,
        "scope": "Exact finite-data algebraic checks, not a proof of the all-word skein realization.",
    }
    write_json(out/'provenance.json', provenance)
    ns = dict(np=np, nx=nx, sp=sp, json=json, hashlib=hashlib, csv=csv, os=os,
              itertools=itertools, time=time, ensure_embedding=ensure_embedding,
              PDCode=PDCode, Yamada=Yamada, A=sp.Symbol('A'),
              CURRENT_BRANCH=provenance['branch'], GIT_COMMIT=provenance['commit'],
              CONSTRUCTOR_VERSION='2026-08-27.fixed-cubic-pure-braid-word.v2')
    for index in [68, 69, 71, 72]:
        execute_cell(notebook,index,ns)
    execute_cell(notebook,86,ns,definitions_only=True)
    basis = ns['TRANSFER_BASIS']
    required = sorted({w for u in basis for v in basis for w in (u+v,u+'A'+v,u+'B'+v)},key=lambda w:(len(w),w))
    exhaustive = [''.join(w) for n in range(8) for w in itertools.product('AB',repeat=n)]
    data_words = sorted(set(required) | set(exhaustive), key=lambda w:(len(w),w))
    plan = [{'length':m,'design':design,'word':word} for m in ns['HELDOUT_LENGTHS'] for design,word in ns['heldout_words'](m)]
    write_json(out/'word_plans.json', {'basis':basis,'required_short_words':required,'exhaustive_words_through_length_7':exhaustive,'all_direct_short_words':data_words,'long_word_plan':plan,
               'short_words_count':len(required),'maximum_short_length':max(map(len,required)),
               'historical_execution_or_freeze_claim':False})
    rows = []
    checkpoint = out/'short_exact_records.json'
    checkpoint_key = hashlib.sha256(canonical({k:provenance[k] for k in ['notebook_sha256','source_file_sha256','packages']}).encode()).hexdigest()
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        if saved['checkpoint_key'] == checkpoint_key:
            rows=saved['records']
    existing={r['word']:r for r in rows}
    for i,word in enumerate(data_words,1):
        if word not in existing:
            existing[word]=ns['evaluate_word_raw_direct'](word)
            write_json(checkpoint,{'checkpoint_key':checkpoint_key,'records':[existing[w] for w in data_words if w in existing]})
        if i%25==0 or i==len(data_words):
            print(f'SHORT {i}/{len(data_words)}',flush=True)
    Y=sp.Symbol('Y')
    coeffs={w:{int(k):int(v) for k,v in json.loads(existing[w]['raw_coefficients_json']).items()} for w in data_words}
    count_classes={}
    for w in exhaustive:
        count_classes.setdefault((w.count('A'),w.count('B')),[]).append(w)
    nonreversal_collisions=[]
    reversal_mismatches=[]
    for words in count_classes.values():
        for i,w in enumerate(words):
            if coeffs[w] != coeffs[w[::-1]]: reversal_mismatches.append(w)
            for v in words[i+1:]:
                if coeffs[w] == coeffs[v] and v != w[::-1]:
                    nonreversal_collisions.append([w,v])
    write_json(out/'exhaustive_short_word_checks.json',{
        'words':exhaustive,'same_count_raw_nonreversal_collisions':nonreversal_collisions,
        'raw_reversal_mismatches':reversal_mismatches,
        'all_raw_reversals_match':not reversal_mismatches,
        'raw_same_count_collisions_only_reversal':not nonreversal_collisions})
    def raw(w):
        return sp.Add(*[sp.Integer(c)*Y**p for p,c in coeffs[w].items()])
    matrices={name:sp.Matrix([[raw(u+mid+v) for v in basis] for u in basis]) for name,mid in [('H',''),('H_A','A'),('H_B','B')]}
    write_json(out/'hankel_input_matrices.json',{'variable':'Y','basis':basis,**{k:[[str(x) for x in row] for row in m.tolist()] for k,m in matrices.items()}})
    print('Checking exact rank certificate at Y=2',flush=True)
    H2=matrices['H'].subs(Y,2)
    det2=H2.det(method='domain-ge')
    assert det2!=0
    certificate={'H_shape':[15,15],'H_rank_over_QY':15,'rank_certificate':{'method':'nonzero exact determinant at Y=2','specialization':'2','determinant':str(det2)},
                 'short_exact_words':len(data_words),'hankel_training_words':len(required),
                 'exhaustive_words_through_length_7':len(exhaustive),'all_raw_reversals_match':not reversal_mismatches,
                 'raw_same_count_collisions_only_reversal':not nonreversal_collisions,
                 'maximum_short_word_length':max(map(len,required)),
                 'symbolic_transfer_checked':False,'all_word_realization_proved':False,'historical_discovery_provenance_verified':False}
    write_json(out/'certificate.json',certificate)
    if not args.symbolic:
        return
    H=DomainMatrix.from_Matrix(matrices['H']).to_field(); F=H.domain
    HA=DomainMatrix.from_Matrix(matrices['H_A']).convert_to(F)
    HB=DomainMatrix.from_Matrix(matrices['H_B']).convert_to(F)
    print('Inverting H over '+str(F),flush=True)
    inv=H.inv(); I=DomainMatrix.eye((15,15),F)
    assert H.matmul(inv).sub(I).is_zero_matrix
    TA=inv.matmul(HA); TB=inv.matmul(HB)
    assert H.matmul(TA).sub(HA).is_zero_matrix
    assert H.matmul(TB).sub(HB).is_zero_matrix
    L=H.extract([0],list(range(15)));R=I.extract(list(range(15)),[0])
    print('Checking cubic identities',flush=True)
    for name,T in [('A',TA),('B',TB)]:
        product=I
        for exponent in [2,-2,-4]:
            product=product.matmul(T.sub(I.scalarmul(F.convert(Y**exponent))))
        assert product.is_zero_matrix,name
    commutator=TA.matmul(TB).sub(TB.matmul(TA)); assert not commutator.is_zero_matrix
    witness=next((i,j,str(v)) for i,row in enumerate(commutator.to_list()) for j,v in enumerate(row) if v!=F.zero)
    print('Checking all short-word matrix contractions',flush=True)
    states={'':L}
    for i,word in enumerate(data_words,1):
        if word:
            # Recompute any absent prefix; required Hankel words need not be prefix closed.
            prefix=word[:-1]
            if prefix not in states:
                state=L
                for letter in prefix:state=state.matmul(TA if letter=='A' else TB)
                states[prefix]=state
            states[word]=states[prefix].matmul(TA if word[-1]=='A' else TB)
        assert states[word].matmul(R).to_list()[0][0]==F.convert(raw(word)),word
        if i%50==0: print(f'TRANSFER {i}/{len(data_words)}',flush=True)
    formula={'variable':'Y','field':str(F),'basis':basis,'H':payload(H),'H_A':payload(HA),'H_B':payload(HB),
             'H_inverse':payload(inv),'T_A':payload(TA),'T_B':payload(TB),'L_transpose':payload(L),'R':payload(R),
             'historical_freeze_verified':False,'run_kind':'new_reproduction_audit_not_historical_discovery'}
    formula['sha256_without_hash']=hashlib.sha256(canonical(formula).encode()).hexdigest()
    write_json(out/'exact_transfer_matrices.json',formula)
    certificate.update(symbolic_transfer_checked=True, H_inverse_identity=True,H_TA_equals_HA=True,H_TB_equals_HB=True,
                       cubic_A_zero=True,cubic_B_zero=True,commutator_nonzero=True,commutator_witness=witness,
                       all_short_transfer_contractions_match=True,formula_sha256=formula['sha256_without_hash'])
    write_json(out/'certificate.json',certificate)
    freeze={'recorded_utc':datetime.now(timezone.utc).isoformat(),'formula_sha256':formula['sha256_without_hash'],
            'planned_words':plan,'scope':'new audit freeze only; all words and formulas are previously published in local manuscript/notebook'}
    write_json(out/'new_audit_freeze.json',freeze)
    tests=[]
    for i,item in enumerate(plan[:args.heldout],1):
        word=item['word'];state=L
        for letter in word:state=state.matmul(TA if letter=='A' else TB)
        predicted=state.matmul(R).to_list()[0][0]
        print(f'HELDOUT {i}/{args.heldout} length={len(word)} {item["design"]}: direct evaluation',flush=True)
        actual=ns['evaluate_word_raw_direct'](word)
        actual_expr=sp.Add(*[sp.Integer(v)*Y**int(k) for k,v in json.loads(actual['raw_coefficients_json']).items()])
        passed=predicted==F.convert(actual_expr)
        tests.append({**item,'raw_coefficient_identity':passed,'actual':actual})
        write_json(out/'new_long_word_checks.json',tests)
        assert passed,item
    certificate['new_long_word_checks']=len(tests)
    certificate['new_long_word_checks_all_pass']=all(x['raw_coefficient_identity'] for x in tests) if tests else None
    write_json(out/'certificate.json',certificate)
    print('COMPLETE',flush=True)


if __name__=='__main__':
    main()
