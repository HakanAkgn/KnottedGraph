#!/usr/bin/env python3
"""Render unfiltered finite-box Betti numbers and extraction status from pilot rows."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm,ListedColormap
from matplotlib.patches import Patch
def main():
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 plan=json.loads((a.data/'plan.json').read_text());rows=json.loads((a.data/'records.json').read_text());families=list(plan['families']);lambdas=np.array(plan['lambda_array']);allgammas=np.array(plan['candidate_gammas'])
 names=['Hopf link\n→ trefoil','Hopf link\n→ Solomon link','Unknot\n→ trefoil','Unknot\n→ Solomon link','Trefoil\n→ cinquefoil']
 statuses=['not-computed','betti-consistent-extraction','unsupported-boundary','unsupported-cavity','empty','unavailable-extraction']
 labels=['Not evaluated','Graph passes Betti checks','Touches box boundary','Cavity: graph spine unsupported','Empty mask','Extraction unavailable']
 colors=['#d9d9d9','#278066','#eeaa44','#7654a3','#aec7e8','#bc4545']
 grids={};maxb=[0,0]
 for key in families:
  shape=(plan['families'][key],len(lambdas));b0=np.full(shape,np.nan);b1=b0.copy();b2=b0.copy();status=np.zeros(shape,int)
  for r in rows:
   if r['transition']!=key:continue
   ij=(r['gamma_index'],r['lambda_index']);status[ij]=statuses.index(r['status']);v=r.get('volume')
   if v:
    b0[ij]=v['components'];b1[ij]=v['handle_rank'];b2[ij]=v['enclosed_voids']
  grids[key]={'beta0':b0,'beta1':b1,'beta2':b2,'status':status}
  for i,b in enumerate([b0,b1]):
   if np.isfinite(b).any():maxb[i]=max(maxb[i],int(np.nanmax(b)))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'pdf.fonttype':42,'ps.fonttype':42})
 fig,axs=plt.subplots(3,5,figsize=(17.5,8.5),layout='constrained')
 for col,key in enumerate(families):
  gammas=allgammas[:plan['families'][key]]
  for row,field in enumerate(['beta0','beta1','status']):
   ax=axs[row,col];arr=grids[key][field]
   if row<2:
    cmap=plt.get_cmap('Blues' if row==0 else 'viridis',maxb[row]+1).copy();cmap.set_bad('#d9d9d9');norm=BoundaryNorm(np.arange(maxb[row]+2)-.5,cmap.N);mesh=ax.pcolormesh(lambdas,gammas,np.ma.masked_invalid(arr),cmap=cmap,norm=norm,shading='nearest',rasterized=True)
   else:
    cmap=ListedColormap(colors);norm=BoundaryNorm(np.arange(len(colors)+1)-.5,len(colors));mesh=ax.pcolormesh(lambdas,gammas,arr,cmap=cmap,norm=norm,shading='nearest',rasterized=True)
   ax.set_xlim(0,1);ax.set_ylim(float(gammas[0]),float(gammas[-1]));ax.set_xticks([0,.5,1]);ax.tick_params(labelsize=13)
   if row==0:ax.set_title(names[col],fontsize=18,pad=9)
   if row==2:ax.set_xlabel(r'$\lambda$',fontsize=19)
   if col==0:ax.set_ylabel([r'$\beta_0$: energy $E$',r'$\beta_1$: energy $E$','Status: energy $E$'][row],fontsize=17)
   if col==4 and row<2:
    cbar=fig.colorbar(mesh,ax=axs[row,:].tolist(),fraction=.015,pad=.01,ticks=np.arange(maxb[row]+1));cbar.ax.tick_params(labelsize=13)
 present=set(r['status'] for r in rows)
 if len(rows)<plan['planned_cells']:present.add('not-computed')
 legend=[Patch(facecolor=colors[i],label=labels[i]) for i,s in enumerate(statuses) if s in present]
 fig.legend(handles=legend,loc='outside lower center',ncols=3,fontsize=13,frameon=False)
 fig.suptitle(f'Finite-window Hamiltonian pilot: {len(rows):,} samples at '+r'$64^3$ voxels',fontsize=21)
 pdf=a.out/'Hamiltonian_volume_homology_pilot.pdf';png=a.out/'Hamiltonian_volume_homology_pilot.png';fig.savefig(pdf);fig.savefig(png,dpi=160);plt.close(fig)
 arrays={key+'_'+field:v for key,g in grids.items() for field,v in g.items()};np.savez_compressed(a.out/'figure_arrays.npz',lambdas=lambdas,candidate_gammas=allgammas,**arrays)
 metadata={'completed_cells':len(rows),'planned_cells':plan['planned_cells'],'source_records_sha256':hashlib.sha256((a.data/'records.json').read_bytes()).hexdigest(),'status_index':dict(enumerate(statuses)),'status_counts':dict(Counter(r['status'] for r in rows)),'filtered':False,'yamada_classification':False,'continuum_convergence_established':False,'beta0_max':maxb[0],'beta1_max':maxb[1],'figure_pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()}
 (a.out/'figure_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n');print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()

