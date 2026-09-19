#!/usr/bin/env python3
"""Assemble the complete text-only manuscript and check source invariants.

Only the document parts in this directory are read. No figure, existing paper,
scientific implementation or original bibliography is modified. The static
checks are not described as a successful full-project LaTeX compilation.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import zipfile

PARTS = ('00_main.tex', '10_representation.tex', '20_projection_yamada.tex',
         '30_validation.tex', '35_parameter_scans.tex', '40_layout_discovery.tex',
         '50_architecture_limits.tex')
FIGURES = {
    'GeneralVersion/Figures/InputFormatsOverview_Pipeline.pdf': r'width=\linewidth,keepaspectratio',
    'GeneralVersion/Figures/Provided_Functionalities.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/KnottedGraphVsTopoly.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/PorousMaterialPhaseMap.pdf': r'width=0.88\linewidth',
    'GeneralVersion/Figures/MathematicalGraphPatterns/LLMBasedFormulaDiscovery_Yamada.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/InputYamadaExamples.pdf': r'width=\linewidth,keepaspectratio',
    'GeneralVersion/Figures/InputSkeletonizationBeyondYamada.pdf': r'width=\linewidth,keepaspectratio',
    'GeneralVersion/Figures/SkeletonizationSteps.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/PDCodeGeneration.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/PDcode_to_yamada.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/YamadaResolutions.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/HandleBody_YamadaFramework_Stresstest.pdf': r'width=\linewidth,keepaspectratio',
    './GeneralVersion/Figures/Time_distributions.pdf': r'width=\linewidth,keepaspectratio',
    'GeneralVersion/Figures/PhaseDIagrams.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/Repulsive_curves.pdf': r'width=\linewidth',
    'GeneralVersion/Figures/MathematicalGraphPatterns/NonAbelianWorkedExample.pdf': r'width=0.62\linewidth',
    'GeneralVersion/Figures/Architecture.pdf': r'width=\linewidth',
}
TIKZ = r'''\begin{tikzpicture}[scale=0.9, line cap=round, line join=round]

% Reidemeister-II pair
\draw[line width=0.8pt]
    (0,0.45) -- (0.8,-0.45) -- (1.6,0.45) -- (2.4,0.45);
\draw[line width=0.8pt]
    (0,-0.45) -- (0.8,0.45) -- (1.6,-0.45) -- (2.4,-0.45);

% First crossing: upper strand
\draw[white,line width=3.2pt] (0.30,0.1125) -- (0.50,-0.1125);
\draw[line width=0.8pt]       (0.30,0.1125) -- (0.50,-0.1125);

% Second crossing: opposite strand is upper
\draw[white,line width=3.2pt] (1.10,0.1125) -- (1.30,-0.1125);
\draw[line width=0.8pt]       (1.10,0.1125) -- (1.30,-0.1125);

% Arrow
\draw[->,line width=0.8pt] (2.8,0) -- (3.8,0)
    node[midway,above] {$\mathrm{RII}$};

% After cancellation
\draw[line width=0.8pt] (4.2,0.45) -- (6.6,0.45);
\draw[line width=0.8pt] (4.2,-0.45) -- (6.6,-0.45);

\end{tikzpicture}'''
FORMULA_LABELS = (
    'supp:deriv:eq:factors', 'supp:deriv:eq:omega', 'supp:deriv:eq:transfers',
    'supp:deriv:eq:lambda_family', 'supp:deriv:eq:down_family', 'supp:deriv:eq:updown_family',
    'supp:deriv:eq:mixed_word_counts', 'supp:deriv:eq:mixed_commuting_transfer',
    'supp:deriv:eq:mixed_channels', 'supp:deriv:eq:abelian_final',
    'supp:deriv:eq:pure_braid_word', 'supp:deriv:eq:observed_reversal',
    'supp:deriv:eq:hankel_basis', 'supp:deriv:eq:H_definition', 'supp:deriv:eq:hankel_rank_15',
    'supp:deriv:eq:hankel_matrices', 'supp:deriv:eq:hankel_examples', 'supp:deriv:eq:TA_TB',
    'supp:deriv:eq:raw_word_transfer', 'supp:deriv:eq:AAABA_state_path', 'supp:deriv:eq:cubic',
    'supp:deriv:eq:projector_general', 'supp:deriv:eq:run_power', 'supp:deriv:eq:run_decomp',
    'supp:deriv:eq:raw_nonabelian_final', 'supp:deriv:eq:normalization',
    'supp:deriv:eq:AAABA_projector_prediction', 'supp:deriv:eq:TA3_reduction',
    'supp:deriv:eq:AAABA_reduced', 'supp:deriv:eq:AAABA_basis_contractions',
    'supp:deriv:eq:AAABA_prediction',
)


def uncomment(text):
    lines = []
    for line in text.splitlines():
        cut = len(line)
        for index, char in enumerate(line):
            if char != '%':
                continue
            slashes = 0
            prior = index-1
            while prior >= 0 and line[prior] == '\\':
                slashes += 1
                prior -= 1
            if slashes % 2 == 0:
                cut = index
                break
        lines.append(line[:cut])
    return '\n'.join(lines)


def check_braces(text):
    depth, index = 0, 0
    while index < len(text):
        char = text[index]
        if char == '\\':
            match = re.match(r'\\(?:[A-Za-z]+|.)', text[index:])
            if match:
                index += len(match.group())
                continue
        depth += (char == '{') - (char == '}')
        if depth < 0:
            raise ValueError('unmatched closing brace near '+text[max(0,index-60):index+30])
        index += 1
    if depth:
        raise ValueError('unmatched opening brace count '+str(depth))


def check_environments(text):
    stack, counts = [], Counter()
    for match in re.finditer(r'\\(begin|end)\{([^}]+)\}', text):
        operation, name = match.groups()
        if operation == 'begin':
            stack.append(name)
            counts[name] += 1
        elif not stack or stack.pop() != name:
            raise ValueError('mismatched environment near '+match.group())
    if stack:
        raise ValueError('unclosed environments '+repr(stack))
    return dict(counts)


def cross_references(text, labels):
    patterns = (
        (r'\\(?:Figref|Figpanelref|SuppFigref|SuppFigpanelref|SuppNoteref|SuppSubsecref|SuppTabref|Eqref|ref\*?|eqref)\s*\{([^{}]+)\}', 1),
        (r'\\(?:Eqsref|SuppSubsecsref)\s*\{([^{}]+)\}\s*\{([^{}]+)\}', 2),
        (r'\\hyperref\[([^\]]+)\]', 1),
    )
    targets = []
    for expression, count in patterns:
        for match in re.finditer(expression, text):
            targets += [match.group(i) for i in range(1,count+1) if not match.group(i).startswith('#')]
    missing = sorted(set(targets)-set(labels))
    if missing:
        raise ValueError('undefined static cross-reference targets '+repr(missing))
    return sorted(set(targets))


def test_bibliography_helper(root):
    path = root/'apply_available_bib_fixes.py'
    spec = importlib.util.spec_from_file_location('bib_text_fixes',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sample = '''@article{kg_galeski2022lifshitz, title={A {nested, title}}, pages={bad}, year={2022}}
@article{kg_shi2017lifshitz, title="Other, title", year={2017}}
@article{kg_li2019density, title={Density}, doi={old}, year={2019}}
@misc{untouched, title={Do not edit}, year={2020}}
'''
    changed, history = module.revise(sample)
    assert len(history)==3
    assert 'pages = {7418}' in changed and 'pages = {14988}' in changed
    assert 'doi = {10.1134/S0081543819030076}' in changed
    assert '@misc{untouched, title={Do not edit}, year={2020}}' in changed
    assert 'title={A {nested, title}}' in changed and 'title="Other, title"' in changed
    repeated, _ = module.revise(changed)
    assert repeated == changed
    try:
        module.revise('@article{other, title={Missing targets}}')
    except ValueError:
        pass
    else:
        raise AssertionError('missing bibliography targets silently accepted')


def build(root, out):
    root, out = Path(root).resolve(), Path(out).resolve()
    out.mkdir(parents=True,exist_ok=False)
    parts = [(name,(root/'parts'/name).read_text()) for name in PARTS]
    text = ''.join(content for name,content in parts)
    clean = uncomment(text)
    check_braces(clean)
    environments = check_environments(clean)
    labels = re.findall(r'\\label\{([^{}]+)\}',clean)
    duplicates = [key for key,count in Counter(labels).items() if count>1]
    if duplicates:
        raise ValueError('duplicate labels '+repr(duplicates))
    targets = cross_references(clean,labels)
    if not set(FORMULA_LABELS).issubset(labels):
        raise ValueError('a retained family-formula block is missing')
    graphics = re.findall(r'\\includegraphics\s*(?:\[([^]]*)\])?\s*\{([^{}]+)\}',clean)
    observed = {path:re.sub(r'\s+','',options) for options,path in graphics}
    if len(graphics)!=17 or observed!=FIGURES:
        raise ValueError('figure path/options changed: '+repr(observed))
    drawing = re.findall(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}',text,re.S)
    if drawing != [TIKZ]:
        raise ValueError('the retained inline drawing changed')
    if clean.index('GeneralVersion/Figures/PorousMaterialPhaseMap.pdf') < clean.index('\\BeginSupplementaryInformation'):
        raise ValueError('TPMS artwork remains in the main text')
    if clean.index('\\label{supp:limitations}') < clean.index('\\label{supp:software_architecture}'):
        raise ValueError('limitations did not follow software architecture')
    abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',clean,re.S).group(1)
    plain = re.sub(r'\\texttt\{([^{}]+)\}',r'\1',abstract).replace('$','')
    plain = re.sub(r'\s+',' ',plain).strip()
    if '\\' in plain or not plain.isascii() or len(plain)>1920:
        raise ValueError('plain abstract is not ASCII or exceeds 1920 characters: '+str(len(plain)))
    citation_keys = sorted({key.strip() for group in re.findall(r'\\cite(?:\[[^]]*\])*\{([^{}]+)\}',clean)
                            for key in group.split(',')})
    test_bibliography_helper(root)
    (out/'Third_draft.tex').write_text(text)
    (out/'Third_draft.txt').write_text(text)
    (out/'arxiv_abstract.txt').write_text(plain+'\n')
    for name in ('text_only_references.bib','README.md','AUDIT_RESPONSE.md','BIBLIOGRAPHY_NOTES.md','apply_available_bib_fixes.py'):
        shutil.copyfile(root/name,out/name)
    shutil.copyfile(Path(__file__),out/'assemble_text_revision.py')
    shutil.copytree(root/'parts',out/'parts')
    qa = {'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
          'full_manuscript_bytes':len(text.encode()),'full_manuscript_lines':len(text.splitlines()),
          'source_sha256':sha256(text.encode()).hexdigest(),
          'parts_sha256':{name:sha256(content.encode()).hexdigest() for name,content in parts},
          'all_17_original_figure_paths_and_options_preserved':True,
          'figure_asset_files_written_or_modified':False,
          'source_pdf_graphic_bytes_independently_rehashed':False,
          'inline_tikz_matches_transcribed_original_verbatim':True,
          'inline_tikz_sha256':sha256(TIKZ.encode()).hexdigest(),
          'retained_family_formula_labels':list(FORMULA_LABELS),
          'mathematical_expressions_preserved_by_editorial_review_not_a_universal_proof':True,
          'balanced_braces':True,'balanced_environments':environments,
          'label_count':len(labels),'duplicate_labels':duplicates,
          'static_cross_reference_targets':targets,'undefined_static_targets':[],
          'citation_keys':citation_keys,'original_bibliography_available_for_key_validation':False,
          'abstract_characters':len(plain),'abstract_ascii':True,
          'remaining_author_placeholders':re.findall(r'\\placeholder\{([^{}]+)\}',text),
          'tpms_figure_moved_to_supplement':True,'limitations_after_software_architecture':True,
          'bibliography_helper_synthetic_tests_passed':True,
          'full_original_project_compiled':False,'new_manuscript_pdf_created':False,
          'missing_original_build_dependencies':['GeneralVersion/macros.tex or root macros.tex','GeneralVersion/references.bib','17 unchanged standalone figure PDFs'],
          'no_arxiv_submission_performed':True}
    (out/'TEXT_REVISION_QA.json').write_text(json.dumps(qa,sort_keys=True,indent=2)+'\n')
    manifest = {str(p.relative_to(out)):sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
    (out/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    archive = out.parent/'KnottedGraph_text_only_revision.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as saved:
        for path in sorted(out.rglob('*')):
            if path.is_file():
                saved.write(path,'KnottedGraph_text_only_revision/'+str(path.relative_to(out)))
    print('TEXT_REVISION_QA '+json.dumps(qa,sort_keys=True),flush=True)
    print('TEXT_REVISION_PACKAGE '+json.dumps({'path':str(archive),'bytes':archive.stat().st_size,
          'sha256':sha256(archive.read_bytes()).hexdigest()}),flush=True)
    return qa


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    build(Path(__file__).resolve().parent,args.out)
