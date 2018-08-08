
import sys
import os
import subprocess
import shlex


sim_name = sys.argv[1]

os.chdir('./{0}/ogrids'.format(sim_name))

f_grid_list = os.listdir(os.getcwd())

for f_grid in f_grid_list:
    with open('{0}.convert.out'.format(sim_name),'w') as convert_out:
        subprocess.call(shlex.split('../../../mnt/sda_flex/bin/convert_grid resid3d.bin.grd.{0} resid3d.grd.{0} -convert'.format(f_grid.split('.')[-1])), stdout=convert_out)
