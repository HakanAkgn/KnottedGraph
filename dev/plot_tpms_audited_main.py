"""Generate the unfiltered TPMS main figure from audited cell records."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as effects
import numpy as np


def edges(values):
    a = np.asarray(values)
    return np.r_[a[0], (a[:-1]+a[1:])/2, a[-1]]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--records', type=Path, required=True)
    p.add_argument('--output-prefix', type=Path, required=True)
    args = p.parse_args()
    records = json.loads(args.records.read_text())
    order = [('schwarz_p_to_diamond','Schwarz-P → Diamond'),
             ('gyroid_to_schwarz_p','Gyroid → Schwarz-P'),
             ('gyroid_to_diamond','Gyroid → Diamond')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.7})
    fig, axes = plt.subplots(2,3,figsize=(11.4,7.2),sharex=True,sharey=True)
    fig.subplots_adjust(left=.065,right=.93,bottom=.25,top=.88,wspace=.16,hspace=.36)
    colors = ['#2166ac','#e5a93d','#b5b7bc','#fbf4ee']
    source_ids = {'yamada':0,'vertex':0,'diagram-yamada':1,'large-core':2,'error':3}
    status_cmap = ListedColormap(colors)
    status_norm = BoundaryNorm(np.arange(-.5,4.5),4)
    vmax=max(r['handle_rank'] for r in records)
    figure_data = {'source_records_sha256':hashlib.sha256(args.records.read_bytes()).hexdigest(),
        'raw_unfiltered':True,'family_order':[v[0] for v in order],
        'interpretation':'Finite-voxel homology and evaluation scope, not embedding-equivalence classes.',
        'families':{}}
    for col,(family,title) in enumerate(order):
        rows=[r for r in records if r['family']==family]
        ls=sorted({r['lam'] for r in rows});cs=sorted({r['threshold_c'] for r in rows})
        lookup={(round(r['lam'],12),round(r['threshold_c'],12)):r for r in rows}
        grid=[[lookup[round(l,12),round(c,12)] for l in ls] for c in cs]
        b0=np.array([[r['interior_components'] for r in row] for row in grid])
        b1=np.array([[r['handle_rank'] for r in row] for row in grid])
        b2=np.array([[r['void_components'] for r in row] for row in grid])
        status=np.array([[source_ids[r['source']] for r in row] for row in grid])
        ex,ey=edges(ls),edges(cs)
        top=axes[0,col];bottom=axes[1,col]
        image=top.pcolormesh(ex,ey,b1,cmap='viridis',vmin=0,vmax=vmax,shading='flat',rasterized=False)
        bottom.pcolormesh(ex,ey,status,cmap=status_cmap,norm=status_norm,shading='flat',rasterized=False)
        for i,j in zip(*np.nonzero(b2)):
            for ax in (top,bottom):
                ax.add_patch(Rectangle((ex[j],ey[i]),ex[j+1]-ex[j],ey[i+1]-ey[i],
                    facecolor='none',edgecolor='#83322e',linewidth=.45,hatch='////',zorder=4))
        for i,j in zip(*np.nonzero(b0>1)):
            top.text((ex[j]+ex[j+1])/2,(ey[i]+ey[i+1])/2,str(b0[i,j]),
                     ha='center',va='center',fontsize=6.8,color='white',fontweight='bold',
                     path_effects=[effects.withStroke(linewidth=1.5,foreground='#202020')])
        for row,ax in enumerate((top,bottom)):
            ax.set_xlim(0,1);ax.set_ylim(0,.3)
            ax.set_xticks([0,.25,.5,.75,1]);ax.set_xticklabels(['0','.25','.50','.75','1'])
            ax.set_yticks([0,.1,.2,.3]);ax.tick_params(length=3,width=.7)
            ax.text(-.04,1.03,chr(ord('a')+row*3+col),transform=ax.transAxes,
                    ha='left',va='bottom',fontweight='bold',fontsize=12)
            if col==0:ax.set_ylabel(r'Level offset $c$')
            if row==1:ax.set_xlabel(r'Interpolation $\lambda$')
        top.set_title(title,fontsize=11.3,pad=19)
        bottom.set_title('Evaluation status',fontsize=10,pad=8)
        counts=Counter(r['source'] for r in rows)
        figure_data['families'][family]={'lambdas':ls,'thresholds':cs,'b0_grid':b0.tolist(),
            'b1_grid':b1.tolist(),'b2_grid':b2.tolist(),'status_grid':status.tolist(),
            'source_counts':dict(counts)}
    top_box=axes[0,2].get_position()
    cax=fig.add_axes([.945,top_box.y0,.014,top_box.height])
    cb=fig.colorbar(image,cax=cax,ticks=[0,5,10,15,20,25,30,35])
    cb.set_label(r'Voxel first Betti number $\beta_1$',fontsize=10)
    fig.text(.065,.985,'Finite-resolution TPMS homology and evaluation status',fontsize=14,fontweight='bold')
    fig.text(.065,.145,r'Top row: $\beta_0=1$ unless marked by a numeral; hatching marks $\beta_2>0$.',fontsize=10)
    handles=[Patch(facecolor=colors[0],label='Normalized subcubic Yamada polynomial'),
             Patch(facecolor=colors[1],label='Fixed-diagram polynomial (valence > 3)'),
             Patch(facecolor=colors[2],label='Abstract graph summary (above evaluation limit)'),
             Patch(facecolor=colors[3],edgecolor='#83322e',hatch='////',label='Cavities: graph-spine reconstruction unsupported')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.055,.018),ncol=2,
               frameon=False,fontsize=9.6,columnspacing=2.8,handlelength=2.3,labelspacing=.75)
    for suffix in ('.pdf','.png'):
        path=Path(str(args.output_prefix)+suffix)
        fig.savefig(path,dpi=220,bbox_inches='tight',facecolor='white')
        print(path)
    Path(str(args.output_prefix)+'_source_data.json').write_text(json.dumps(figure_data,indent=2))
    plt.close(fig)


if __name__=='__main__':main()
