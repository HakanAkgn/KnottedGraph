#!/usr/bin/env python3
"""Replay a historical-nine output using independent rational trigonometry."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import hashlib
import inspect
import json
from pathlib import Path
from time import perf_counter

from knotted_graph.core.field_isotopy import tpms_problem
from knotted_graph.core.field_isotopy_rational import verify_rational


def work(payload):
    root, row = payload
    problem = tpms_problem(row['family'], row['lambda_bounds'], row['c_bounds'])
    path = Path(root)/row['certificate']
    raw = path.read_bytes()
    cert = json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)
    started = perf_counter()
    result = verify_rational(problem, cert)
    return {'family':row['family'],'lambda_bounds':row['lambda_bounds'],'c_bounds':row['c_bounds'],
            'certificate':row['certificate'],'certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'elapsed_seconds':perf_counter()-started, **result}


def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,choices=range(1,5),default=2);a=p.parse_args()
    if a.out.exists():p.error('use a new output file')
    path=a.data/'records.json'
    if path.exists(): rows=json.loads(path.read_text())
    else: rows=json.loads((a.data/'summary.json').read_text())
    rows=[r for r in rows if r['status']=='certified'];results=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futures=[ex.submit(work,(str(a.data),r)) for r in rows]
        for f in as_completed(futures):
            r=f.result();results.append(r);print(r['family'],r['lambda_bounds'],r['valid'],flush=True)
            a.out.write_text(json.dumps({'complete':len(results)==len(rows),'checked':len(results),'planned':len(rows),'all_valid':all(q['valid'] for q in results),'backend_source_sha256':hashlib.sha256(Path(inspect.getfile(verify_rational)).read_bytes()).hexdigest(),'results':results},indent=2)+'\n')


if __name__=='__main__':main()
