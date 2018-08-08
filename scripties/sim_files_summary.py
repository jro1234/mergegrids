
import sys
import os



def count_file_lines(f_name, start_line=0):
    if abs(float(start_line) - int(start_line)) < 1e-8:
        line_count = -int(start_line)
    else:
        print "Must provide integer value for line start, {0} is not".format(start_line)
        return
    try:
        with open(f_name, 'r') as f_in:
            for line in f_in:
                line_count += 1
    except IOError as e:
        print "Could not open file {0}\nError Message".format(f_name)
        print e
        return
    return line_count


sim_name = sys.argv[1]
os.chdir('../runs/{0}'.format(sim_name))
dirs = os.listdir(os.getcwd())

for dir in dirs:

    if not os.path.isdir(dir):
        continue

#    f_complexes = '{0}/complexes'.format(dir)
#    complexes_size = os.path.getsize(f_complexes)
#    complexes_count = count_file_lines(f_complexes, 4)

    f_bootstrap = '{0}/assoc_bootstrap.log'.format(dir)
    bootstrap_size = os.path.getsize(f_bootstrap)
    bootstrap_count = count_file_lines(f_bootstrap, 8)

#    f_grid = '{0}/resid3d.bin.grd'.format(dir)
#    grid_size = os.path.getsize(f_grid)

#    f_trajectories = '{0}/trajectories'.format(dir)
#    trajectories_size = os.path.getsize(f_trajectories)

    print dir
#    print 'complexes\tfile size:\t   ',complexes_size,'\t  num records:\t   ',complexes_count
    print 'bootstrap\tfile size:\t   ',bootstrap_size,'\t  num records:\t   ',bootstrap_count
#    print 'grid     \tfile size:\t   ',grid_size
#    print 'trjctries\tfile size:\t   \n',trajectories_size
