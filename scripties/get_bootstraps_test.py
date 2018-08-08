
import sys
import os
from numpy import array as np_array
from numpy.linalg import norm as np_norm


sim_name = sys.argv[1]
n_run = sys.argv[2]

os.chdir('./analyze_runs/{0}'.format(sim_name))

dirs = os.listdir(os.getcwd())
count = 0
with open('assoc_bootstrap.log','w') as data_out:
    for traj_dir in dirs:
        if not os.path.isdir(traj_dir):
            continue
        os.chdir(traj_dir)
        print '\n'+traj_dir
        with open('assoc_bootstrap.log','r') as data_in:
            if count == 0:
                count += 1
                for i in range(7):
                    data_out.write(next(data_in))
                data_out.write(next(data_in)[:14]+n_run+'\n')
            else:
#                try:
#                    [next(data_in) for i in range(8)]
#                    for line in data_in:
#                        data_out.write(line)
#                        count += 1
#                except:
                center_loc = (0,0,0)
                final_distances_strings = dict()
                with open('trajectories','r') as traj_in:
                    [next(traj_in) for i in range(2)]
                    line3 = next(traj_in)
                    line4 = next(traj_in)
                    site_loc = np_array([float(line4.split()[i]) for i in range(1,4)]) - np_array([float(line3.split()[i]) for i in range(1,4)])
                    for entry_line in traj_in:
                        entry_list = entry_line.split()
                        coords = entry_list[2:5]
                        traj_num = int(entry_list[0])
                        current_loc = np_array([float(coords[i]) for i in range(3)])
                        if traj_num not in final_distances_strings.keys():
                            final_distances_strings.update({traj_num:None})
                        final_distances_strings[traj_num] = current_loc
                    for key in final_distances_strings.keys():
                        site_dist = np_norm(final_distances_strings[key] - site_loc)
                        entry = -1
                        if site_dist < 52.49:
                            entry = int(site_dist)
                        print entry
                        data_out.write(entry)

        print "number of lines:  {0}".format(count)
        os.chdir('../')
