#!/usr/bin/python

from pylab import *
from numpy import *
from os import listdir
from os import path
import matplotlib.pyplot as plt
import re

data_dir = "../traj_data" 
max_num_file = 10
# font size
fs = 30

# find all trajectory files  
files = [f for f in listdir(data_dir) if  re.match( r'traj\_[0-9]+\.txt', f)]

# define the key used to sort the name of trajectory files
convert = lambda text : int(text) if text.isdigit() else text 
alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
# sort file names
files.sort(key=alphanum_key)

print "\n%d trajectory files are found in directory: %s:" % (len(files), data_dir) 
num_file = min( max_num_file, len(files) )
print "the following (%d) files will be used : " % (num_file) 
for f in files[:num_file] :
  print f

# read the dimension of system from the first file
traj_file = open("%s/%s" % (data_dir, files[0]), 'r')
dim, = [ int(x) for x in traj_file.readline().split() ] 
traj_file.close()

# read the time vector from each trajectory file
t_vec = []
for f in files[:num_file] :
  traj_file = open("%s/%s" % (data_dir, f), 'r')
  itmp, = [ int(x) for x in traj_file.readline().split() ] 
  t_vec_tmp = np.loadtxt(traj_file, usecols=0, unpack=True)
  traj_file.close()
  t_vec.append(t_vec_tmp)

# compute the largest time of all trajectories
tot_t_vec = [(x[-1]) for x in t_vec]
max_t = max(tot_t_vec)
print "\nMaximal time = %.6f\n" % max_t

# print each component
for idx in range(dim) :
    print "print the %dth component... " % (idx+1)

    fig = figure()
    ax = gca() 

#   read each trajectory file
    f_idx = 0
    for f in files[:num_file] :
        traj_file = open("%s/%s" % (data_dir, f), 'r')
        itmp, = [ int(x) for x in traj_file.readline().split() ] 
#    assert(itmp == dim) 
        x = np.loadtxt(traj_file, usecols=(idx+1), unpack=True)
        traj_file.close()
        plot( t_vec[f_idx], x, linewidth=2, label="traj %d" % f_idx)
        f_idx = f_idx + 1

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.xlim( 0, max_t )

    xlabel('$t$', fontsize=fs)
    ylabel('$x_%d$' % idx, fontsize=fs)
    lg = legend(bbox_to_anchor=(-0.0, -0.05, 1.0, 1.0), prop={'size':18})
    lg.draw_frame(False)
    plt.gcf().subplots_adjust(bottom=0.2) 
   
    savefig("../fig/x_%d_traj.eps" % idx)
