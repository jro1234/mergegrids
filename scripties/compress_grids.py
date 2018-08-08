
import sys
import os
import shlex
import subprocess

sim_name = sys.argv[1]

os.chdir('../runs/{0}'.format(sim_name))
os.mkdir('ogrids')
missing = list()

dirs = os.listdir(os.getcwd())
for traj_dir in dirs:
    print traj_dir
    if not os.path.isdir(traj_dir) or traj_dir.split('/')[-1] == 'ogrids':
        continue
    if os.path.exists('{0}/resid3d.bin.grd'.format(traj_dir)):
        os.rename('{0}/resid3d.bin.grd'.format(traj_dir),'ogrids/resid3d.bin.grd.{0}'.format(traj_dir.split('/')[-1]))
    else:
        missing.append(traj_dir.split('/')[-1])

with open('ogrids/missing.txt','w') as missing_out:
    missing_out.write('List of grids missing from simulation\n\n')
    missing_out.write('Simulation Name:\n             {0}\n'.format(sim_name))
    for missing_file in missing:
        missing_out.write('\n')
        missing_out.write(missing_file)
