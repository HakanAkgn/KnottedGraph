#!/usr/bin/env python3
"""Check compiled result files and package their actual numerical sources.

This performs automated layout/text/data checks and renders every report page
for visual review. It does not describe automated checks as visual approval.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

import fitz


def write(path, data):
    Path(path).write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False)+'\n')


def summarize(directory):
    summary = json.loads((directory/'completion_summary.json').read_text())
    witness = json.loads((directory/'data/continuum_vs_digital_witnesses.json').read_text())
    lines = ['# Completed computation and restored figures', '',
             'This package reports executed numerical results. It is not the unavailable current 42-page manuscript.', '',
             '## Grid outcomes', '', '| Outcome | Hamiltonian | TPMS |', '|---|---:|---:|']
    modes = summary['maps']
    for status in ('evaluated','fixed_diagram_evaluated','time_budget','existing_cavity_route_not_revised'):
        values = [modes[m]['final_statuses'].get(status,0) for m in ('hamiltonian','tpms')]
        lines.append('| '+status+' | '+' | '.join(map(str,values))+' |')
    lines += ['', '## Analytic TPMS rectangles', '',
              '| Family | Regular area (%) | Event-containing-cell area (%) | Unknown area (%) |',
              '|---|---:|---:|---:|']
    for family, data in summary['events']['by_family'].items():
        lines.append('| '+family+' | '+' | '.join(f"{100*data['area_fractions'][k]:.6f}" for k in ('regular','event','unknown'))+' |')
    lines += ['', 'Event outcomes: '+json.dumps(summary['events']['event_statuses'], sort_keys=True)+'.', '',
              'These are regularity and event-existence certificates, not a complete classification into distinct ambient-isotopy classes.', '',
              '## Closed-voxel manifold tests', '', '| Input | Tested | PL manifold | Nonmanifold |','|---|---:|---:|---:|']
    for mode, data in summary['voxel_links'].items():
        lines.append('| '+mode+' | '+' | '.join(str(data[k]) for k in ('cases','pl_manifold_cases','nonmanifold_cases'))+' |')
    lines += ['', 'The nonmanifold count refers to the represented closed voxel complexes, not to all analytic level sets.', '',
              '## Explicit source-to-digital comparison witnesses']
    for category in ('digital_homology_conflict_examples','spine_signature_variation_examples'):
        examples = witness[category]
        lines += ['', '### '+category, '', f'{len(examples)} examples retained (a capped example set, not an exhaustive total).']
        for example in examples[:3]:
            first, second = example['first'], example['second']
            lines += ['', f"Family `{example['family']}`, regular analytic rectangle `{example['regular_leaf_id']}`.",
                      f"First sample: lambda={first['lambda']!r}, c={first['level']!r}; digital Betti tuple {example['first_betti']}.",
                      f"Second sample: lambda={second['lambda']!r}, c={second['level']!r}; digital Betti tuple {example['second_betti']}.",
                      'Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.',
                      'Certificate: `'+example['regularity_certificate']+'`.']
    lines += ['', '## Remaining scope', '',
              'Unavailable polynomial values, higher-valence fixed-diagram values and excluded cavity-route cases remain separate.',
              'Unknown analytic rectangles are not relabeled as inequivalent. A selected-spine polynomial is not a complete source-solid invariant.',
              'The original manuscript source, author metadata and all-word formula work were not replaced.',
              'The existing cavity implementation and protected integration branch were not modified.', '',
              'The figures retain the original grids and representative geometry, with no small-region relabeling or abstract polynomial fallback.',
              'All figure arrays, exact coefficient legends, source hashes, local tests and example witnesses are provided in this package.']
    (directory/'RESULTS_SUMMARY.md').write_text('\n'.join(lines)+'\n')
    pending = []
    for mode in ('hamiltonian','tpms'):
        rows = json.loads((directory/f'data/{mode}_records.json').read_text())
        for row in rows:
            if row['status']=='time_budget':
                retry = row.get('retry_result',{})
                pending.append({'id':row['id'],'family':row['family'],'lambda':row['lambda'],'level':row['level'],
                                'graph':row.get('graph'),'retry_status':retry.get('status'),
                                'retry_error':retry.get('error'),'candidates':retry.get('candidate_projections',[])[:8],
                                'attempts':retry.get('attempts',[])})
    write(directory/'pending_polynomial_details.json',pending)
    print('FINAL_PENDING_POLYNOMIALS '+json.dumps(pending),flush=True)
    print('EXAMPLE_DIGITAL_CONFLICT '+json.dumps(witness['digital_homology_conflict_examples'][:1]),flush=True)
    print('FINAL_RESULT_COUNTS '+json.dumps(summary,sort_keys=True),flush=True)


def check(directory):
    directory = Path(directory)
    log = (directory/'Results_Report.log').read_text(errors='replace')
    fatal = [message for message in ('There were undefined references','There were undefined citations','! LaTeX Error') if message in log]
    if fatal:
        raise RuntimeError('unresolved report build errors: '+repr(fatal))
    output = directory/'QA'
    output.mkdir(exist_ok=True)
    summary = {'automated_checks_only':True,'visual_inspection_claimed':False,'pdf_files':{},
               'latex_overfull_box_count':log.count('Overfull \\hbox'),
               'latex_underfull_box_count':log.count('Underfull \\hbox'),
               'undefined_reference_or_citation_errors':fatal}
    text = []
    files = [directory/'Results_Report.pdf',*sorted((directory/'Figures').glob('*.pdf'))]
    for path in files:
        with fitz.open(path) as pdf:
            if not len(pdf):
                raise ValueError('empty PDF')
            if path.parent.name=='Figures' and len(pdf)!=1:
                raise ValueError('a standalone figure unexpectedly has multiple pages')
            pages = []
            for index,page in enumerate(pdf):
                extracted = page.get_text()
                if path.name=='Results_Report.pdf' and len(extracted.strip())<50:
                    raise ValueError('report contains an empty or near-empty page')
                spans = [span for block in page.get_text('dict')['blocks'] if block.get('type')==0
                         for line in block['lines'] for span in line['spans']]
                outside = [s['text'] for s in spans if fitz.Rect(s['bbox']).x0 < -1
                           or fitz.Rect(s['bbox']).x1 > page.rect.width+1
                           or fitz.Rect(s['bbox']).y0 < -1 or fitz.Rect(s['bbox']).y1 > page.rect.height+1]
                if outside:
                    raise ValueError('text lies outside PDF page: '+repr(outside[:3]))
                pixmap = page.get_pixmap(matrix=fitz.Matrix(1.7,1.7),alpha=False)
                image = output/f'{path.stem}-page-{index+1}.png'
                pixmap.save(image)
                pages.append({'page':index+1,'width_pt':page.rect.width,'height_pt':page.rect.height,
                              'text_characters':len(extracted),'out_of_page_text':outside,
                              'preview':str(image.relative_to(directory))})
                if path.name=='Results_Report.pdf':
                    text.extend([f'--- PAGE {index+1} ---',extracted])
            summary['pdf_files'][str(path.relative_to(directory))] = {
                'pages':pages,'bytes':path.stat().st_size,'sha256':sha256(path.read_bytes()).hexdigest()}
    (directory/'REPORT_TEXT.txt').write_text('\n'.join(text))
    write(directory/'report_QA.json',summary)
    summarize(directory)
    codes = directory/'code'
    codes.mkdir(exist_ok=True)
    for path in sorted(Path('dev').glob('*.py')):
        if any(word in path.name for word in ('completion','voxel_link','tpms','full_map','full_resolution','archived_spatial')):
            shutil.copyfile(path,codes/path.name)
    for path in (Path('dev/full_map_collapse.cpp'),Path('uv.lock'),Path('pyproject.toml')):
        shutil.copyfile(path,codes/path.name)
    (codes/'source_revision.txt').write_text(subprocess.check_output(['git','rev-parse','HEAD'],text=True))
    omitted = {'.aux','.log','.out','.fls','.fdb_latexmk','.synctex.gz'}
    for path in directory.iterdir():
        if path.is_file() and any(path.name.endswith(suffix) for suffix in omitted):
            path.unlink()
    manifest = {str(path.relative_to(directory)):sha256(path.read_bytes()).hexdigest()
                for path in sorted(directory.rglob('*')) if path.is_file() and path.name!='package_manifest.json'}
    write(directory/'package_manifest.json',manifest)
    archive = directory.parent/'KnottedGraph_results_and_figures.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as saved:
        for path in sorted(directory.rglob('*')):
            if path.is_file():
                saved.write(path,'KnottedGraph_results/'+str(path.relative_to(directory)))
    print('RESULT_PACKAGE_QA '+json.dumps(summary),flush=True)
    print('RESULT_PACKAGE_ZIP '+json.dumps({'name':archive.name,'bytes':archive.stat().st_size,
          'sha256':sha256(archive.read_bytes()).hexdigest()}),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    args = parser.parse_args()
    check(args.directory)
