from pathlib import Path
import subprocess
import sys


def test_complete_plot_arrays_and_no_removed_cubical_faces():
    root = Path(__file__).resolve().parents[2]
    code = '''
import numpy as np
from render_completion_results import point_grid,cubical_faces,cover_grids
rows=[]
for iy in range(2):
    for ix in range(2):
        rows.append({'family':'test','level_index':iy,'lambda_index':ix,'lambda':float(ix),'level':float(iy),
                     'status':'evaluated','yamada':'A+1','polynomial_terms':[[0,1],[1,1]]})
rows[-1]['status']='fixed_diagram_evaluated'
g=point_grid(rows,'test')
assert g['labels']==[[1,1],[1,-1]]
assert g['colored_cells']==3 and g['uncolored_cells']==1
assert len(cubical_faces(np.ones((1,1,1),dtype=bool),(np.arange(2),)*3))==6
assert len(cubical_faces(np.ones((2,1,1),dtype=bool),(np.arange(3),np.arange(2),np.arange(2))))==10
print('exact source masks and plotting-category distinctions passed')
'''
    subprocess.run([sys.executable,'-c',code],cwd=root/'dev',check=True)
