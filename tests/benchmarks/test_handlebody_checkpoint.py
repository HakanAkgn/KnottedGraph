"""Checkpoint identity must survive neither a seed change nor a code change."""
import ast
import copy
import hashlib
import json
from pathlib import Path

import pytest


def helpers():
    checkout = Path(__file__).resolve().parents[2]
    notebook = checkout / 'User_guide/benchmarks/04_thick_handlebody_validation.ipynb'
    if not notebook.exists():
        notebook = Path(__file__).with_name('04_thick_handlebody_validation.ipynb')
    content = json.loads(notebook.read_text())
    functions = []
    for cell in content['cells']:
        if cell['cell_type'] == 'code':
            functions.extend(node for node in ast.parse(''.join(cell['source'])).body
                             if isinstance(node, ast.FunctionDef) and node.name in
                             {'graph_payload_sha256', 'cached_result_is_current'})
    namespace = {'json': json, 'hashlib': hashlib}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(notebook), 'exec'), namespace)
    return namespace


def valid_cache():
    ns = helpers()
    record = {'seed_graph_sha256': 'seed1', 'recovery_source_sha256': 'code1',
              'recovered_graph': {'nodes': [], 'edges': []}, 'recovered_yamada': '1'}
    digest = ns['graph_payload_sha256'](record['recovered_graph'])
    record.update(recovered_graph_sha256=digest, recovered_yamada_graph_sha256=digest)
    return ns['cached_result_is_current'], copy.deepcopy(record), record


def test_exact_seed_recovered_and_source_identity_can_resume():
    check, row, record = valid_cache()
    assert check(row, record, 'seed1', 'code1')


@pytest.mark.parametrize('mutation', ['seed', 'source', 'graph', 'row_hash', 'polynomial_hash', 'missing_identity'])
def test_changed_or_unverifiable_checkpoint_is_recomputed(mutation):
    check, row, record = valid_cache()
    seed, source = 'seed1', 'code1'
    if mutation == 'seed':
        seed = 'seed2'
    elif mutation == 'source':
        source = 'code2'
    elif mutation == 'graph':
        record['recovered_graph']['nodes'].append({'id': 0, 'pos': [0, 0, 0]})
    elif mutation == 'row_hash':
        row['recovered_graph_sha256'] = 'different'
    elif mutation == 'polynomial_hash':
        row['recovered_yamada_graph_sha256'] = 'different'
    else:
        row.pop('seed_graph_sha256')
    assert not check(row, record, seed, source)
