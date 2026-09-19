#!/usr/bin/env python3
"""Publish only generated result documents to the dedicated research branch.

No forced update, protected-branch write, source-code overwrite or secret output.
The remote branch must still point to the exact checked-out workflow commit.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

REPOSITORY = 'HakanAkgn/KnottedGraph'
BRANCH = 'research/solid-spine-maps-41dc87c'
PREFIX = 'doc/reproducibility/completion_results/'


def api(method,path,data=None):
    body = None if data is None else json.dumps(data).encode()
    request = Request('https://api.github.com/repos/'+REPOSITORY+'/'+path,
                      data=body,method=method,headers={
                          'Authorization':'Bearer '+os.environ['GH_TOKEN'],
                          'Accept':'application/vnd.github+json','Content-Type':'application/json',
                          'X-GitHub-Api-Version':'2022-11-28','User-Agent':'KnottedGraph-result-publisher'})
    with urlopen(request,timeout=60) as response:
        return json.load(response)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args = parser.parse_args()
    if os.environ.get('GITHUB_REPOSITORY')!=REPOSITORY or os.environ.get('GITHUB_REF_NAME')!=BRANCH:
        raise RuntimeError('publication is restricted to the requested repository and research branch')
    head = os.environ['GITHUB_SHA']
    ref_path = 'git/refs/heads/'+quote(BRANCH,safe='/')
    current = api('GET',ref_path)
    if current['object']['sha']!=head:
        raise RuntimeError('branch changed since the build; refusing to overwrite newer work')
    qa = json.loads((args.directory/'report_QA.json').read_text())
    if qa['undefined_reference_or_citation_errors']:
        raise ValueError('report QA contains unresolved references')
    names = ['Results_Report.pdf','Results_Report.tex','manuscript_insert.tex','numbers.tex',
             'table_outcomes.tex','table_events.tex','table_links.tex','provenance.tex',
             'completion_summary.json','figure_source_data.json','README.md','RESULTS_SUMMARY.md',
             'REPORT_TEXT.txt','report_QA.json','package_manifest.json','pending_polynomial_details.json']
    names += [str(p.relative_to(args.directory)) for p in sorted((args.directory/'Figures').glob('*'))
              if p.suffix in ('.pdf','.png')]
    names += ['data/continuum_vs_digital_witnesses.json']
    paths = [args.directory/name for name in names]
    if not all(p.is_file() and not p.is_symlink() for p in paths):
        raise ValueError('missing publication file')
    if sum(p.stat().st_size for p in paths)>80*1024*1024:
        raise ValueError('result publication exceeds the explicit size bound')
    commit = api('GET','git/commits/'+head)
    entries = []
    for name,path in zip(names,paths):
        target = PREFIX+name
        if '..' in Path(target).parts or not target.startswith(PREFIX):
            raise ValueError('publication escaped the results directory')
        blob = api('POST','git/blobs',{'content':base64.b64encode(path.read_bytes()).decode(),'encoding':'base64'})
        entries.append({'path':target,'mode':'100644','type':'blob','sha':blob['sha']})
    tree = api('POST','git/trees',{'base_tree':commit['tree']['sha'],'tree':entries})
    created = api('POST','git/commits',{
        'message':'Publish compiled full-grid results report and restored source-linked figures\n\nGenerated from the recorded frozen inputs after complete tests, data validation,\nPDF build and automated layout checks. Preserve original source files, cavity\nimplementation, formula work and protected integration branch.',
        'tree':tree['sha'],'parents':[head]})
    if api('GET',ref_path)['object']['sha']!=head:
        raise RuntimeError('branch advanced while preparing results; no ref was changed')
    api('PATCH',ref_path,{'sha':created['sha'],'force':False})
    receipt = {'commit':created['sha'],'parent':head,'branch':BRANCH,'paths':[e['path'] for e in entries],
               'protected_integration_branch_written':False}
    args.receipt.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print('PUBLISHED_RESULTS '+json.dumps(receipt),flush=True)


if __name__=='__main__':
    main()
