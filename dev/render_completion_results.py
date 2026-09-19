#!/usr/bin/env python3
"""Render source-linked categorical maps and retain complete representative geometry.

The original three-map/five-geometry TPMS design is restored in two explicitly
separate versions: selected-spine polynomial signatures and certified analytic
continuation components. Neither figure uses abstract contraction grouping or
small-region relabeling. Every exact exposed voxel face and graph point is
retained. Rasterization of the surface artist changes display encoding only.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.ndimage import label

HAMILTONIAN = ('hopf_to_trefoil','hopf_to_solomon','unknot_to_trefoil','unknot_to_solomon','trefoil_to_cinquefoil')
TPMS = ('schwarz_p_to_diamond','gyroid_to_schwarz_p','gyroid_to_diamond')
NAMES = {'hopf_to_trefoil':'Hopf link → Trefoil','hopf_to_solomon':'Hopf link → Solomon link',
         'unknot_to_trefoil':'Unknot → Trefoil','unknot_to_solomon':'Unknot → Solomon link',
         'trefoil_to_cinquefoil':'Trefoil → Cinquefoil','schwarz_p_to_diamond':'Schwarz-P → Diamond',
         'gyroid_to_schwarz_p':'Gyroid → Schwarz-P','gyroid_to_diamond':'Gyroid → Diamond'}
BASE_COLORS = ['#ffaa0e','#d62728','#1f77b4','#9467bd','#2ca02c']
VIEWS = [(26.,-64.),(20.,-30.),(28.,-48.),(18.,-76.),(30.,-24.)]
MISSING = '#e2e5e9'
EVENT = '#bfc5cc'
UNKNOWN = '#ffffff'


def dump(path, value):
    Path(path).write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def palette(count):
    if count <= 5:
        return BASE_COLORS[:max(1,count)]
    if count <= 20:
        return [mpl.colors.to_hex(plt.get_cmap('tab20')(i)) for i in range(count)]
    return [mpl.colors.to_hex(plt.get_cmap('turbo')(v)) for v in np.linspace(.04,.96,count)]


def sample_edges(values):
    values = np.asarray(values,dtype=float)
    if len(values)<2 or np.any(np.diff(values)<=0):
        raise ValueError('strictly ordered sample coordinates required')
    return np.r_[values[0],(values[:-1]+values[1:])/2,values[-1]]


def point_grid(rows, family):
    selected = [r for r in rows if r['family']==family]
    xs = sorted({r['lambda'] for r in selected})
    ys = sorted({r['level'] for r in selected})
    if len(selected)!=len(xs)*len(ys):
        raise ValueError('incomplete original plotting grid')
    lookup = {(r['level_index'],r['lambda_index']):r for r in selected}
    if len(lookup)!=len(selected):
        raise ValueError('duplicate plotting input')
    # Labels are local to each panel. Equality is exact coefficient equality,
    # not similarity of strings or abstract-graph isomorphism.
    available = [r for r in selected if r['status']=='evaluated']
    keys = sorted({tuple(map(tuple,r['polynomial_terms'])) for r in available})
    ids = {key:i+1 for i,key in enumerate(keys)}
    values = np.full((len(ys),len(xs)),-1,dtype=int)
    for (iy,ix),row in lookup.items():
        if row['status']=='evaluated':
            values[iy,ix] = ids[tuple(map(tuple,row['polynomial_terms']))]
    colors = palette(len(keys))
    legend = [{'id':i+1,'terms':[list(t) for t in key],'color':colors[i],
               'polynomial':next(r['yamada'] for r in available if tuple(map(tuple,r['polynomial_terms']))==key)}
              for i,key in enumerate(keys)]
    return {'family':family,'labels':values.tolist(),'x_edges':sample_edges(xs).tolist(),
            'y_edges':sample_edges(ys).tolist(),'sample_x':xs,'sample_y':ys,'legend':legend,
            'statuses':dict(Counter(r['status'] for r in selected)),
            'colored_cells':int((values>0).sum()),'uncolored_cells':int((values<0).sum())}


def refined_edges(values,depth):
    values = list(values)
    for _ in range(depth):
        result = []
        for left,right in zip(values[:-1],values[1:]):
            result.extend((left,left/2+right/2))
        values = result+[values[-1]]
    return values


def cover_grids(leaves, events):
    statuses = {r['id']:r['status'] for r in events}
    depth = max(r['depth'] for r in leaves)
    grids = {}
    for family in TPMS:
        size = 20*2**depth
        coverage = np.zeros((size,size),dtype=np.uint8)
        states = np.full((size,size),-1,dtype=int)
        leaf_indices = np.full((size,size),-1,dtype=int)
        for index,row in enumerate(leaves):
            if row['family']!=family:
                continue
            width = 2**(depth-row['depth'])
            xx,yy = row['tile_x']*width,row['tile_y']*width
            region = np.s_[yy:yy+width,xx:xx+width]
            coverage[region] += 1
            leaf_indices[region] = index
            states[region] = 0 if row['status']=='certified' else (2 if statuses[row['id']]=='unknown' else 1)
        if not np.all(coverage==1):
            raise ValueError('adaptive parameter cover is incomplete or overlapping')
        components,count = label(states==0,structure=np.ones((3,3),dtype=int))
        labels = np.where(states==0,components,np.where(states==1,-1,-2))
        grids[family] = {'family':family,'labels':labels.tolist(),'states':states.tolist(),
                         'leaf_indices':leaf_indices.tolist(),'component_count':int(count),
                         'x_edges':refined_edges(np.linspace(0,1,21),depth),
                         'y_edges':refined_edges(np.linspace(0,.3,21),depth),
                         'colors':palette(int(count))}
    return grids


def boundaries(ax,xs,ys,values):
    values = np.asarray(values)
    segments = []
    for yy,xx in np.argwhere((values[:,:-1]!=values[:,1:]) & (values[:,:-1]>0) & (values[:,1:]>0)):
        segments.append(((xs[xx+1],ys[yy]),(xs[xx+1],ys[yy+1])))
    for yy,xx in np.argwhere((values[:-1,:]!=values[1:,:]) & (values[:-1,:]>0) & (values[1:,:]>0)):
        segments.append(((xs[xx],ys[yy+1]),(xs[xx+1],ys[yy+1])))
    ax.add_collection(LineCollection(segments,colors='black',linewidths=.65))


def draw_map(fig,ax,grid,letter,ylabel,*,cover=False):
    values = np.asarray(grid['labels'],dtype=int)
    count = grid['component_count'] if cover else len(grid['legend'])
    colors = grid['colors'] if cover else [r['color'] for r in grid['legend']]
    if not colors:
        colors = palette(1)
    # Negative categories never receive a polynomial/component color.
    coded = np.where(values>0,values+1,np.where(values==-2,0,1))
    cmap = ListedColormap([UNKNOWN,EVENT if cover else MISSING]+colors)
    norm = BoundaryNorm(np.arange(-.5,len(colors)+2.5),cmap.N)
    image = ax.pcolormesh(grid['x_edges'],grid['y_edges'],coded,shading='flat',cmap=cmap,norm=norm,
                         antialiased=False,rasterized=False)
    boundaries(ax,grid['x_edges'],grid['y_edges'],values)
    ax.set(xlim=(grid['x_edges'][0],grid['x_edges'][-1]),ylim=(grid['y_edges'][0],grid['y_edges'][-1]),
           xlabel=r'$\lambda$',ylabel=ylabel)
    ax.set_title(NAMES[grid['family']],fontsize=13,pad=9)
    ax.text(-.12,1.03,letter,transform=ax.transAxes,fontsize=15,fontweight='bold')
    ax.set_xticks(np.linspace(0,1,6))
    ax.tick_params(labelsize=10,width=1.1,length=3.5)
    if count:
        scalar = mpl.cm.ScalarMappable(norm=BoundaryNorm(np.arange(.5,count+1.5),count),
                                        cmap=ListedColormap(colors))
        bar = fig.colorbar(scalar,ax=ax,fraction=.035,pad=.025,drawedges=True)
        bar.set_ticks([])
        bar.ax.set_title('C' if cover else r'$\Upsilon$',fontsize=13,pad=5)
        bar.dividers.set_linewidth(.4)
    return image


def cubical_faces(mask,axes):
    """Vectorized exact exposed-face enumeration; not marching-cubes decimation."""
    quads = []
    if any(len(a)!=n+1 for a,n in zip(axes,mask.shape)):
        raise ValueError('witness gridlines do not match the occupied-cell mask')
    for axis in range(3):
        other = [a for a in range(3) if a!=axis]
        for side in (0,1):
            neighbor = np.zeros_like(mask)
            source,target = [slice(None)]*3,[slice(None)]*3
            source[axis] = slice(None,-1) if side==0 else slice(1,None)
            target[axis] = slice(1,None) if side==0 else slice(None,-1)
            neighbor[tuple(target)] = mask[tuple(source)]
            indices = np.argwhere(mask & ~neighbor)
            corners = []
            for aa,bb in ((0,0),(1,0),(1,1),(0,1)):
                offset = np.zeros(3,dtype=int)
                offset[axis],offset[other[0]],offset[other[1]] = side,aa,bb
                corner = indices+offset
                corners.append(np.column_stack([axes[k][corner[:,k]] for k in range(3)]))
            quads.append(np.stack(corners,axis=1))
    return np.concatenate(quads,axis=0)


def geometry_data(directory, record):
    raw = (directory/'graph.json.gz').read_bytes()
    if sha256(raw).hexdigest()!=record['graph_sha256']:
        raise ValueError('representative graph hash changed')
    archive = directory/'collapse.npz'
    if sha256(archive.read_bytes()).hexdigest()!=record['reconstruction']['archive_sha256']:
        raise ValueError('representative witness hash changed')
    with np.load(archive,allow_pickle=False) as data:
        mask = data['mask']
        if sha256(mask.astype(np.uint8).tobytes()).hexdigest()!=record['source_mask_sha256']:
            raise ValueError('representative source mask changed')
        axes = [data[key] for key in ('x','y','z')]
        faces = cubical_faces(mask,axes)
    return faces,json.loads(gzip.decompress(raw))


def draw_geometry(ax,faces,graph,color,view):
    surface = Poly3DCollection(faces,facecolors=color,edgecolors='none',linewidths=0,alpha=.16)
    surface.set_rasterized(True)
    ax.add_collection3d(surface)
    for edge in graph['edges']:
        points = np.asarray(edge['pts'])
        ax.plot(*points.T,color='#111111',linewidth=1.25,solid_capstyle='round')
    nodes = np.asarray([n['pos'] for n in graph['nodes']])
    if len(nodes):
        ax.scatter(*nodes.T,color='#d62728',s=11,edgecolors='white',linewidths=.25,depthshade=False)
    limit = max(5.25,float(np.abs(faces).max())*1.015)
    ax.set(xlim=(-limit,limit),ylim=(-limit,limit),zlim=(-limit,limit))
    ax.set_box_aspect((1,1,1),zoom=1.05)
    ax.set_proj_type('ortho')
    ax.view_init(elev=view[0],azim=view[1])
    ax.set_axis_off()


def cover_at(grid,lam,level):
    xs,ys,values = np.asarray(grid['x_edges']),np.asarray(grid['y_edges']),np.asarray(grid['labels'])
    xx = np.flatnonzero((xs[:-1]<=lam)&(lam<=xs[1:]))
    yy = np.flatnonzero((ys[:-1]<=level)&(level<=ys[1:]))
    positives = {int(values[y,x]) for x in xx for y in yy if values[y,x]>0}
    if len(positives)>1:
        raise ValueError('closed regular rectangles sharing one point have inconsistent component labels')
    return next(iter(positives),-1)


def save(fig,path):
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight',pad_inches=.06,dpi=360)
    fig.savefig(path.with_suffix('.png'),bbox_inches='tight',pad_inches=.06,dpi=170)
    plt.close(fig)


def write_tables(out,summary):
    maps = summary['maps']
    text = ['\\begin{tabular}{lrr}\\toprule','Outcome & Hamiltonian & TPMS \\\\ \\midrule']
    statuses = [('Recorded grid points',None),('Normalized spatial / isolated-vertex','evaluated'),
                ('Higher-valence fixed diagram','fixed_diagram_evaluated'),('Uncompleted polynomial','time_budget'),
                ('Existing cavity route, not recomputed','existing_cavity_route_not_revised')]
    for name,key in statuses:
        numbers = [sum(maps[m]['final_statuses'].values()) if key is None else maps[m]['final_statuses'].get(key,0)
                   for m in ('hamiltonian','tpms')]
        text.append(name+' & '+' & '.join(f'{n:,}' for n in numbers)+' \\\\')
    text += ['\\bottomrule\\end{tabular}']
    (out/'table_outcomes.tex').write_text('\n'.join(text))
    text = ['\\begin{tabular}{lrrr}\\toprule','Family & Regular area & Event-cell area & Unknown area \\\\ \\midrule']
    for family in TPMS:
        area = summary['events']['by_family'][family]['area_fractions']
        name = NAMES[family].replace('→',r'$\to$')
        text.append(name+' & '+' & '.join(f'{100*area[k]:.4f}\\%' for k in ('regular','event','unknown'))+' \\\\')
    text += ['\\bottomrule\\end{tabular}']
    (out/'table_events.tex').write_text('\n'.join(text))
    text = ['\\begin{tabular}{lrrr}\\toprule','Source & Cases & PL manifold & Nonmanifold \\\\ \\midrule']
    for mode in ('hamiltonian','tpms'):
        counts = summary['voxel_links'][mode]
        text.append(mode.title()+' & '+' & '.join(f'{counts[k]:,}' for k in ('cases','pl_manifold_cases','nonmanifold_cases'))+' \\\\')
    text += ['\\bottomrule\\end{tabular}']
    (out/'table_links.tex').write_text('\n'.join(text))
    event_counts = summary['events']['event_statuses']
    remaining = event_counts.get('unknown',0)
    positive = sum(v for k,v in event_counts.items() if k!='unknown')
    definitions = {'EventCells':positive,'UnknownCells':remaining,
                   'NewEvents':summary['events']['new_event_certificates_replayed'],
                   'RegularLeaves':summary['events']['regular_leaf_count'],
                   'HomologyConflicts':summary['continuum_vs_digital']['retained_digital_homology_conflicts'],
                   'SignatureVariations':summary['continuum_vs_digital']['retained_spine_signature_variations']}
    (out/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+key+'}{'+f'{value:,}'+'}' for key,value in definitions.items()))


def build(out,data,artifacts):
    out,data,artifacts = Path(out),Path(data),Path(artifacts)
    out.mkdir(parents=True,exist_ok=False)
    figures = out/'Figures'
    figures.mkdir()
    summary = json.loads((data/'completion_summary.json').read_text())
    hrows = json.loads((data/'hamiltonian_records.json').read_text())
    trows = json.loads((data/'tpms_records.json').read_text())
    leaves = json.loads((data/'adaptive_leaves.json').read_text())
    events = json.loads((data/'event_records.json').read_text())
    hg = {key:point_grid(hrows,key) for key in HAMILTONIAN}
    tg = {key:point_grid(trows,key) for key in TPMS}
    cg = cover_grids(leaves,events)
    mpl.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'mathtext.fontset':'cm',
                         'axes.linewidth':1.1,'axes.labelsize':15,'pdf.fonttype':42,'ps.fonttype':42,
                         'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
    fig,axes = plt.subplots(5,1,figsize=(8.4,11.75),layout='constrained')
    for i,(key,ax) in enumerate(zip(HAMILTONIAN,axes)):
        draw_map(fig,ax,hg[key],f'({chr(97+i)})',r'$E$')
    save(fig,figures/'PhaseDIagrams')
    fig = plt.figure(figsize=(14.8,8.0),layout='constrained')
    gs = fig.add_gridspec(2,6,hspace=.11,wspace=.10)
    slots = (gs[0,:2],gs[0,2:4],gs[0,4:],gs[1,1:3],gs[1,3:5])
    for i,(key,slot) in enumerate(zip(HAMILTONIAN,slots)):
        draw_map(fig,fig.add_subplot(slot),hg[key],f'({chr(97+i)})',r'$E$')
    save(fig,figures/'PhaseDIagrams_Compact')
    representatives = []
    geometry = []
    for index,li in enumerate((0,1,5,11,17)):
        row = next(r for r in trows if r['family']=='gyroid_to_diamond' and r['lambda_index']==li and r['level_index']==7)
        source = artifacts/row['artifact_relative_path']
        faces,graph = geometry_data(source,row)
        geometry.append((faces,graph))
        destination = out/'representatives'/row['id']
        destination.mkdir(parents=True)
        for name in ('source.npz','collapse.npz','graph.json.gz'):
            shutil.copyfile(source/name,destination/name)
        p_label = int(tg['gyroid_to_diamond']['labels'][7][li])
        c_label = cover_at(cg['gyroid_to_diamond'],row['lambda'],row['level'])
        representatives.append({'id':row['id'],'lambda':row['lambda'],'c':row['level'],'status':row['status'],
                                'graph_sha256':row['graph_sha256'],'source_mask_sha256':row['source_mask_sha256'],
                                'witness_sha256':row['reconstruction']['archive_sha256'],'boundary_quads':len(faces),
                                'graph_nodes':len(graph['nodes']),'graph_edges':len(graph['edges']),
                                'graph_polyline_points':sum(len(e['pts']) for e in graph['edges']),
                                'polynomial_label':p_label,'continuation_label':c_label,'view':VIEWS[index],
                                'geometry_relative_path':str(destination.relative_to(out)),
                                'faces_or_polyline_points_removed':False})
    for cover,stem in ((False,'PorousMaterialPhaseMap'),(True,'PorousMaterialContinuationMap')):
        fig = plt.figure(figsize=(15.3,7.3),layout='constrained')
        gs = fig.add_gridspec(2,15,height_ratios=(1.,.95),hspace=.04,wspace=.04)
        for i,key in enumerate(TPMS):
            ax = fig.add_subplot(gs[0,i*5:(i+1)*5])
            grid = cg[key] if cover else tg[key]
            draw_map(fig,ax,grid,f'({chr(97+i)})',r'$c$',cover=cover)
            if key=='gyroid_to_diamond':
                ax.axhline(representatives[0]['c'],color='black',lw=.8,ls=(0,(3,2)))
                ax.scatter([r['lambda'] for r in representatives],[r['c'] for r in representatives],
                           color='black',edgecolors='white',linewidths=.5,s=20,zorder=10,clip_on=False)
        for i,(entry,(faces,graph)) in enumerate(zip(representatives,geometry)):
            ax = fig.add_subplot(gs[1,3*i:3*(i+1)],projection='3d')
            key = 'continuation_label' if cover else 'polynomial_label'
            index = entry[key]
            colors = cg['gyroid_to_diamond']['colors'] if cover else [r['color'] for r in tg['gyroid_to_diamond']['legend']]
            color = colors[index-1] if index>0 else '#929aa6'
            draw_geometry(ax,faces,graph,color,VIEWS[i])
            ax.set_title(r'$\lambda='+f"{entry['lambda']:.2f}"+'$',fontsize=14,pad=0)
        save(fig,figures/stem)
    fig,axes = plt.subplots(1,3,figsize=(15.3,3.65),layout='constrained')
    for i,(key,ax) in enumerate(zip(TPMS,axes)):
        draw_map(fig,ax,tg[key],f'({chr(97+i)})',r'$c$')
    save(fig,figures/'TPMS_SpatialSignatures')
    fig,axes = plt.subplots(1,3,figsize=(15.3,3.65),layout='constrained')
    evidence_colors = ['#1f77b4','#ffaa0e','#ffffff']
    for i,(key,ax) in enumerate(zip(TPMS,axes)):
        grid = cg[key]
        ax.pcolormesh(grid['x_edges'],grid['y_edges'],grid['states'],shading='flat',
                      cmap=ListedColormap(evidence_colors),norm=BoundaryNorm([-.5,.5,1.5,2.5],3))
        ax.set(xlim=(0,1),ylim=(0,.3),xlabel=r'$\lambda$',ylabel=r'$c$')
        ax.set_title(NAMES[key],fontsize=13,pad=9)
        ax.text(-.12,1.03,f'({chr(97+i)})',transform=ax.transAxes,fontweight='bold',fontsize=15)
    fig.legend(handles=[Patch(facecolor=c,edgecolor='black',label=t) for c,t in zip(evidence_colors,
               ('Certified regular rectangle','Contains a certified event','Unknown'))],
               loc='outside lower center',ncol=3,frameon=False,fontsize=11)
    save(fig,figures/'TPMS_EventEvidence')
    source = {'hamiltonian':hg,'tpms_signatures':tg,'tpms_continuation':cg,'representatives':representatives,
              'small_region_relabeling':False,'abstract_contraction_grouping':False,
              'labels_are_complete_source_isotopy_classes':False,'source_commit':summary['source_commit']}
    dump(out/'figure_source_data.json',source)
    dump(out/'completion_summary.json',summary)
    shutil.copytree(data,out/'data')
    write_tables(out,summary)
    template = Path(__file__).with_name('completion_report_template.tex')
    shutil.copyfile(template,out/'Results_Report.tex')
    for mode,rows in (('hamiltonian',hrows),('tpms',trows)):
        fields = ('id','family','lambda','level','lambda_hex','level_hex','baseline_status','status','yamada','source_mask_sha256','graph_sha256')
        with (out/f'{mode}_all_cells.csv').open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=fields)
            writer.writeheader()
            writer.writerows({k:r.get(k) for k in fields} for r in rows)
    (out/'provenance.tex').write_text('\\noindent Collector revision: \\texttt{'+summary['source_commit']+'}.\\par\n'+
        'Renderer revision: \\texttt{'+subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()+'}.\\par\n')
    methods = Path(__file__).with_name('completion_manuscript_insert.tex')
    shutil.copyfile(methods,out/'manuscript_insert.tex')
    shutil.copyfile(Path(__file__),out/'render_completion_results.py')
    readme = ['# Full-grid results package', '',
              'Results_Report.pdf explains the numerical results, figures, certificate scope and remaining gaps.',
              'Figures/ contains vector categorical maps and full-geometry representative views. Surface artists are rasterized without deleting faces.',
              'PhaseDIagrams.pdf retains the five vertically stacked source-notebook panels; PhaseDIagrams_Compact.pdf supplies a balanced manuscript alternative.',
              'PorousMaterialPhaseMap.pdf restores three spatial-signature maps plus all five original representative parameter points.',
              'PorousMaterialContinuationMap.pdf uses analytic regularity-cover components instead. These are separate scientific quantities.',
              'No small-region filter or abstract polynomial fallback is used.', '',
              'The current 42-page manuscript source is not included. manuscript_insert.tex is an insertion draft, not a reconstructed complete paper.',
              'The all-word formulas and existing cavity implementation are unchanged.', '',
              'Build the report from this directory with:',
              'latexmk -pdf -no-shell-escape -interaction=nonstopmode -halt-on-error Results_Report.tex', '',
              'All full record arrays, event certificates, voxel-link tests and example comparison witnesses are retained in data/.',
              'The five complete source/witness/graph archives are retained under representatives/.',
              'Input-run IDs and frozen environment are documented in the report and repository workflow.', '',
              'Scientific counts:',json.dumps(summary,sort_keys=True,indent=2)]
    (out/'README.md').write_text('\n'.join(readme)+'\n')
    dump(out/'package_manifest.json',{str(p.relative_to(out)):sha256(p.read_bytes()).hexdigest()
                                    for p in out.rglob('*') if p.is_file()})
    print('FIGURE_PACKAGE '+json.dumps({'figures':[p.name for p in figures.glob('*.pdf')],
          'representatives':representatives,'summary':summary},sort_keys=True),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--artifacts',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    arguments = parser.parse_args()
    build(arguments.out,arguments.data,arguments.artifacts)
