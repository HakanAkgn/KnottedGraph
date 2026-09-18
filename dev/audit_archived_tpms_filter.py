"""Extract archived contraction grids and reproduce every changed display cell."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from legacy_tpms_filter import stable_labels


def main():
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--html',type=Path)
    group.add_argument('--extracted-input',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.html:
        raw=args.html.read_bytes()
        match=re.search(r'const payload = (.*?);\nconst transitions',raw.decode())
        if not match:raise ValueError('Archived payload not found')
        payload=json.loads(match.group(1))
        source={'archived_html_filename':args.html.name,'archived_html_sha256':hashlib.sha256(raw).hexdigest(),
            'scope':'Archived abstract-contraction groupings and display changes only; not spatial equivalence.',
            'families':{r['key']:{'lambdas':r['lambdas'],'thresholds':r['thresholds'],
                'raw_contraction_grid':r['modes']['contraction']['z'],
                'archived_display_grid':r['modes']['stable_contraction']['z']} for r in payload['transitions']}}
    else:
        source=json.loads(args.extracted_input.read_text())
    for f,data in source['families'].items():
        raw=np.asarray(data['raw_contraction_grid'],dtype=int)
        displayed=np.asarray(data['archived_display_grid'],dtype=int)
        recreated,edits=stable_labels(raw,min_cells=4)
        if not np.array_equal(recreated,displayed):
            raise ValueError('Frozen filter does not reproduce archived grid: '+f)
        changed=raw!=displayed
        data.update({'threshold':4,'reproduced_display_exactly':True,
            'raw_label_count':len(set(raw.ravel())),'display_label_count':len(set(displayed.ravel())),
            'reassignment_operations':edits,'distinct_changed_cells':int(changed.sum()),
            'changed_cell_mask':changed.tolist(),
            'changed_cells':[{'row':int(i),'col':int(j),'lambda':data['lambdas'][j],
                'c':data['thresholds'][i],'raw_label':int(raw[i,j]),'display_label':int(displayed[i,j])} for i,j in zip(*np.nonzero(changed))]})
    source['total_distinct_changed_cells']=sum(x['distinct_changed_cells'] for x in source['families'].values())
    source['frozen_filter_sha256']=hashlib.sha256(Path(__file__).with_name('legacy_tpms_filter.py').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(source,indent=2))
    print(json.dumps({f:{k:d[k] for k in ('raw_label_count','display_label_count','distinct_changed_cells','reassignment_operations')} for f,d in source['families'].items()},indent=2))


if __name__=='__main__':main()
