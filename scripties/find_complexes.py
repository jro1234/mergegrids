
import sys
import os


sim_name = sys.argv[1]

os.chdir('./analyze_runs/{0}'.format(sim_name))

dirs = os.listdir(os.getcwd())
for dir in dirs:
    if not os.path.isdir(dir):
        continue
    size = os.path.getsize('{0}/complexes'.format(dir))
    print dir,'  ',size
