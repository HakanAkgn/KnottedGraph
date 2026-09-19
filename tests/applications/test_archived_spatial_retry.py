from pathlib import Path
import subprocess
import sys


def test_archived_spatial_retry_preserves_geometry_and_validates_hashes():
    root=Path(__file__).resolve().parents[2]
    code='''
from pathlib import Path
import tempfile, json, gzip
from hashlib import sha256
from retry_archived_spatial_evaluations import load_graph, polynomial
from knotted_graph.projection.pd_code import sample_projections
with tempfile.TemporaryDirectory() as tmp:
    d=Path(tmp)
    data={'nodes':[{'id':0,'pos':[0,0,0]},{'id':1,'pos':[1,0,0]}],
          'edges':[{'u':0,'v':1,'key':0,'pts':[[0,0,0],[.5,0,0],[1,0,0]]}]}
    raw=gzip.compress(json.dumps(data).encode())
    (d/'graph.json.gz').write_bytes(raw)
    (d/'collapse.npz').write_bytes(b'test')
    row={'graph_sha256':sha256(raw).hexdigest(),'reconstruction':{'archive_sha256':sha256(b'test').hexdigest()},
         'graph':{'vertices':2,'edges':1,'components':1,'cycle_rank':0,'max_degree':1}}
    graph=load_graph(d,row)
    assert len(graph[0][1][0]['pts']) == 3
    views=sample_projections(graph,num_rotation_samples=4)
    assert polynomial(views[0])['terms'] == polynomial(views[1])['terms']
    (d/'graph.json.gz').write_bytes(raw+b'corrupt')
    try: load_graph(d,row)
    except ValueError: pass
    else: raise AssertionError('corrupted retained graph accepted')
'''
    subprocess.run([sys.executable,'-c',code],cwd=root/'dev',check=True)
