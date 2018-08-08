

import multiprocessing
import sys
import os
import time


sim_name = sys.argv[1]
n_tmp_files = 15
n_files_proc = 12


def merge_grids(n_grids,n_tmp_files,n_lines_file,grid_list=None,core_idx=0,n_procs=1,grids_header='nothing\nlives\nin\nthis\ncave\n'):
    #################################################################################
    #
    #  
    #  
    #
    #################################################################################
    tmp_list = []
    master = False
    for j in range(n_tmp_files):
        tmp_list.append([])
        for i in range(n_procs):
            tmp_list[j].append('resid3d.grd.{0}.tmp.{1}'.format(i+core_idx,j))
    if not grid_list:
        grid_list = tmp_list
        out_list = ['resid3d.grd.merged']
        master = True
        core_idx = -1
        print "--Starting master process to merge to intermediate files"
#        print grid_list
    else:
        print "--Starting process {0} to average values from files:\n   '{1}'".format(core_idx,"'\n   '".join(grid_list))
        grid_list = [grid_list]
        out_list = [item for sublist in tmp_list for item in sublist]
#    print "\n\n--->test, exiting early\n\n"
#    return
    for ii in range(1,1+len(grid_list)):
        f_block = grid_list[ii-1]
        ready = False
        any_grids = True
        f_block_state = None
#        annoyed = 0
        while not ready:
            f_block_state_fresh = []
            for f_grid in f_block:
                if os.path.isfile(f_grid):
                    f_size = os.stat(f_grid).st_size
                    f_block_state_fresh.append(f_size)
            print f_block_state_fresh
            if not f_block_state:
                print "Process {0} Validating input files ready for merging".format(core_idx)
#                annoyed += 1
#                if annoyed > 50:
#                    return
                f_block_state = f_block_state_fresh
                if master:
                    time.sleep(10)
            elif master and set(f_block_state) == set(f_block_state_fresh):
                print "Input files ready for final merger, block {0} of {1}".format(ii,len(grid_list))
                ready = True
            elif not master and len(set(f_block_state_fresh)) == 1:
                print "Block {0} of {1} ready for merger".format(ii,len(grid_list))
                ready = True
            else:
                print "Waiting for all input files to finalize, block {0} of {1}".format(ii,len(grid_list))
                f_block_state = f_block_state_fresh
                time.sleep(10)
        try:
            f_grids = []
            n_closed = []
            for i in range(len(f_block)):
                f_grids.append(open(f_block[i], 'r'))
                for j in range(5):
                    next(f_grids[i])
            print "Grids open for reading values"
            for out_file in out_list:
#            for out_file in [out_list[0]]:
                if os.path.isfile(out_file):
                    mode = 'a'
                else:
                    mode = 'w'
                with open(out_file, mode) as f_out:
                    print "Opened file {0} for writing".format(out_file)
#                    print grids_header
                    if mode == 'w':
                        f_out.write(grids_header)
                    n_lines_written = 0
                    while (n_lines_written < n_lines_file):
#                    while (n_lines_written < 5230):
                        grid_sum = [0.,0.,0.,0.,0.,0.]
                        len_grid_line = 6
                        new_line = None
                        for f_grid in f_grids:
#                        for f_grid in f_grids[:2]:
                            try:
                                line = next(f_grid)
                                if line.strip() == "0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00":
#                                    print "line of zeros"
                                    continue
                                else:
#                                   The following block should work, Error variables have a local / global namespace issue for some reason in this block
#                                    try:
 #                                       grid_sum += [float(line.split()[i])/n_grids for i in range(len_grid_line)]
  #                                  except IndexError:
   #                                     try:
    #                                        if [int(line.split()[i]) for i in range(1,3)] == [im,jm]:
     #                                           new_line = line
       #                                 except ValueError,IndexError:
        #                                    #  This will detect first instance of grid values line with 0 < entries < 6 (at end of blocks)
         #                                   #  Afterwards falls in first try: statement as len(grid_sum) has now changed from 6
          #                                  len_grid_line = len(line.split())
           #                                 grid_sum = [float(line.split()[i])/n_grids for i in range(len_grid_line)]   '''
                                    if len(line.split()) == len_grid_line:
#                                        print "non-zero line"
                                        grid_sum = [grid_sum[i] + float(line.split()[i])/n_grids for i in range(len_grid_line)]
                                    elif len(line.split()) == 3 and [int(line.split()[i]) for i in range(1,3)] == [im,jm]:
                                        new_line = line
#                                        print "how many times", f_grid, n_lines_written
                                    else:
                                        #  This will detect first instance of grid values line with 0 < entries < 6 (at end of blocks)
                                        #  Afterwards falls in first try: statement as len(grid_sum) has now changed from 6
#                                        print "end of block type line"
                                        len_grid_line = len(line.split())
                                        if len_grid_line > 1:
                                            grid_sum = [float(line.split()[i])/float(n_grids) for i in range(len_grid_line)]
                                        elif len_grid_line == 1:
                                            grid_sum=[float(line.strip())/float(n_grids)]
                            except StopIteration:
                                print "Past last line of UHBD files".format(f_grid)
                                raise
                        if not new_line:
                            if grid_sum == [0.,0.,0.,0.,0.,0.]:
                                new_line = line
                            else:
                                new_line = ''
#                                print grid_sum, n_lines_written
                                for i in range(len_grid_line):
                                    new_line += 'E'.join(['  {:.5e}'.format(grid_sum[i]).split('e')[j] for j in range(2)])
                                new_line += '\n'
                        n_lines_written += 1
#                        print "end file cycle, {0} lines written".format(n_lines_written)
                        f_out.write(new_line)
        except IOError as e:
            print "--Something went wrong: {0}\nwith opening file: {1}".format(e.strerror,f_grid)
            return
        except StopIteration as e:
            print "--\n--\nStopIteration\n--Closing all files opened by process {0}, reached last line.\n--List of file closures:\n".format(core_idx)
            for index,f_grid in enumerate(f_grids):
                #print f_grid
                f_grid.close()
                print f_grid.name,'\n',f_grid
                n_closed.append(index)
            if len(n_closed) == len(f_grids):
                any_grids = False
        finally:
            print "--\n--\nGeneral Closure Check\n--Closing all files still opened by process {0}.\n--List of file closures:\n".format(core_idx)
            if any_grids:
                print "Odd that these grids were not already closed:\n"
                for index in set(range(len(f_grids))).difference(set(n_closed)):
                    print f_grid
                    f_grid.close()
                    print f_grid
            else:
                print "Everything seems cleared up with block {0}".format(ii)
        if master:
            for f_grid in f_grids:
                os.remove(f_grid.name)
                print "Removed file {0}".format(f_grid.name)

###
#  Go to grids directory "sim_name" within analyze_grids folder,
#  go to "ogrids" folder and put ascii grids into "f_grid_list"
###
os.chdir('../{0}/ogrids'.format(sim_name))
f_list = os.listdir(os.getcwd())
f_grid_list = list()
for i in range(len(f_list)):
    if f_list[i].startswith('resid3d.grd') and '.tmp.' not in f_list[i]:
        f_grid_list.append(f_list[i])
###
#  get "n_grids", "n_procs"
#  divy up blocks of "f_grid_list" with "f_grid_breakup"
#  "n_lines_block" from header calculation of file lines
###
n_grids = len(f_grid_list)
n_cores = multiprocessing.cpu_count()
n_procs = n_cores - 1
inc_grids_list = n_grids / n_procs
f_grid_breakup = [i*inc_grids_list for i in range(n_procs)]
f_grid_breakup.append(n_grids)
#print f_grid_breakup
new_breakup = [f_grid_breakup[0]]
for i in range(1,len(f_grid_breakup)):
    len_segment = f_grid_breakup[i]-f_grid_breakup[i-1]
    for j in range(len_segment/n_files_proc):
        new_breakup.append(n_files_proc*(j+1)+f_grid_breakup[i-1])
    new_breakup.append(f_grid_breakup[i])
f_grid_breakup = new_breakup
n_procs = len(f_grid_breakup)-1
#print f_grid_breakup, n_procs
###
#  --Based off f77 code by David Sept from
#  https://searchcode.com/codesearch/raw/22774144/
#
#  "grids_header" read once to get important 'metadata'
#  
#  "im" dimension of grid in x-direction
#  "jm" dimension of grid in y-direction
#  "km" dimension of grid in z-direction - dimension of value block organization in file
#  "h"  spacing grid points
#  "ox" x-coord of grid origin in PDB
#  "oy" y-coord of grid origin in PDB
#  "oz" z-coord of grid origin in PDB
#
#  origin of grid is ?lower-left-front?
###
grids_header = str()
with open(f_grid_list[0],'r') as read_header:
    for i in range(5):
        header_line = next(read_header)
        grids_header += header_line
        if i == 2:
            im, jm, km = [int(header_line.split()[j]) for j in range(3)]
#            h = float(header_line[21:33])
 #           ox = float(header_line[33:45])
  #          oy = float(header_line[45:57])
   #         oz = float(header_line[57:69])
###
#  "n_lines_block" from header calculation of file lines, excludes header
###
n_blocks = km
n_lines_block = (im*jm)/6 + 1
if (im*jm)%6 > 0:
    n_lines_block += 1
n_lines = n_lines_block*n_blocks
######
#   completely arbitrary, works for "n_tmp_files" 10 but should be fixed up a bit
######
n_lines_file = (n_lines / n_tmp_files) + (n_lines % n_tmp_files) / 2 + 2


print "\n--Merging {0} grids from simulation {1}\n  by plain value-averaging at each point".format(n_grids,sim_name)
print "-\n--This will use {0} cores, {1} to merge original files and 1 master to integrate".format(n_cores,n_procs)
print "-\n--Grids are of x,y,z dimensions {0},{1},{2}".format(im,jm,km)
print "-\n--Will use {0} intermediate file segments for master process to integrate".format(n_tmp_files)
print "-\n--This gives about {0} lines printed per each intermediate file\n-\n--\n---\n--\n-\n-".format(n_lines_file)


procs = [multiprocessing.Process(target=merge_grids, args=(n_grids,n_tmp_files,n_lines_file,f_grid_list[f_grid_breakup[i]:f_grid_breakup[i+1]],i)) for i in range(n_procs)]

# Start initial averaging processes
for i in range(n_procs):
    procs[i].start()

# Start the master integrator process
#   n_tmp_files,n_lines_file,grid_list=None,core_idx=0,n_procs=1,grids_header='nothing\nlives\nin\nthis\ncave\n'
#print n_procs
master = multiprocessing.Process(target=merge_grids, args=(1,n_tmp_files,n_lines_file,None,0,n_procs,grids_header))
master.start()



