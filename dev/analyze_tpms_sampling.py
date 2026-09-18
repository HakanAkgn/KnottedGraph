"""Audit nested parameter-grid subsampling and voxel resolution; no convergence claim."""
import argparse
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path


def regions(grid, ignored=None):
    seen = set()
    count = 0
    sizes = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i,j) in seen or grid[i][j] == ignored:
                continue
            count += 1
            queue = deque([(i,j)])
            seen.add((i,j))
            size = 0
            while queue:
                a,b = queue.popleft()
                size += 1
                for c,d in ((a+1,b),(a-1,b),(a,b+1),(a,b-1)):
                    if 0 <= c < len(grid) and 0 <= d < len(grid[0]) and (c,d) not in seen and grid[c][d] == grid[i][j]:
                        seen.add((c,d)); queue.append((c,d))
            sizes.append(size)
    return {"connected_regions": count, "singleton_regions": sizes.count(1)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True, type=Path)
    p.add_argument("--stage-records", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    source = json.loads(args.records.read_text())
    families = defaultdict(list)
    for r in source: families[r['family']].append(r)
    results = {}
    for family, records in families.items():
        lambdas = sorted({r['lam'] for r in records})
        thresholds = sorted({r['threshold_c'] for r in records})
        lookup = {(r['threshold_c'],r['lam']):r for r in records}
        results[family] = {}
        for stride in (1,2,4):
            ls, cs = lambdas[::stride], thresholds[::stride]
            grid = [[lookup[c,lam] for lam in ls] for c in cs]
            selected = [r for row in grid for r in row]
            betti_grid = [[(r['interior_components'],r['handle_rank'],r['void_components']) for r in row] for row in grid]
            signatures = [[r['phase_signature'] if r['source']!='error' else None for r in row] for row in grid]
            results[family][str(len(ls))] = {
                "cells":len(selected), "lambda_spacing":ls[1]-ls[0],"c_spacing":cs[1]-cs[0],
                "betti_tuple_count":len({value for row in betti_grid for value in row}),
                "betti_regions_4neighbor":regions(betti_grid),
                "nonerror_signature_count":len({r['phase_signature'] for r in selected if r['source']!='error'}),
                "nonerror_signature_regions_4neighbor":regions(signatures),
                "cavity_positive_cells":sum(r['void_components']>0 for r in selected),
                "unresolved_cells":sum(r['source']=='error' for r in selected),
                "selected_lambda_indices":list(range(0,len(lambdas),stride)),
                "selected_c_indices":list(range(0,len(thresholds),stride)),
            }
    stages = json.loads(args.stage_records.read_text())['records']
    matching = defaultdict(dict)
    for r in stages:
        matching[r['family'],r['lambda'],r['c']][r['dimension']] = [r['mask'][b] for b in ('b0','b1','b2')]
    resolution = {
        "matched_parameter_cells":len(matching),
        "family_endpoint_duplicates_included":True,
        "pairwise_changed_cells":{f'{a}_{b}':sum(v[a]!=v[b] for v in matching.values()) for a,b in ((24,64),(64,96),(24,96))},
        "changed_between_any":sum(len({tuple(x) for x in v.values()})>1 for v in matching.values()),
        "all_matched_cells":[{'family':f,'lambda':lam,'c':c,'betti_by_dimension':v} for (f,lam,c),v in matching.items()],
    }
    output = {
        "interpretation":"Nested subsampling measures loss of resolved features at coarser parameter spacing; it is not independent parameter-grid convergence. Voxel Betti variation measures the masks themselves, independently of graph extraction.",
        "source_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.records,args.stage_records)},
        "parameter_subsampling":results,"voxel_resolution":resolution,
    }
    args.output.write_text(json.dumps(output,indent=2))
    print(json.dumps({"parameter_subsampling":results,"resolution_pairwise_changes":resolution['pairwise_changed_cells']},indent=2))


if __name__=='__main__': main()
