#!/usr/bin/python

"""
  This code plots the copy-number of each component in the trajectory data as a function of
  time.

"""

from pylab import *
from numpy import *
from os import listdir
from os import path
import matplotlib.pyplot as plt
import re

data_dir = "../traj_data" 
max_num_file = 5
lc = ['b', 'k', 'g', 'r', 'y']
# font size
fs = 30

# the final time of each trajectory
max_t = 10
print "\nMaximal time = %.6f\n" % max_t

# find all trajectory files  
file_name_idx = [68, 16, 2, 7, 15] 
files = [ "%s/traj_%d.txt" % (data_dir, i) for i in file_name_idx]

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
  t_vec_tmp = np.loadtxt(traj_file, usecols=(0,), unpack=True)
  t_vec_tmp = np.append( t_vec_tmp, [max_t] )
  traj_file.close()
  t_vec.append(t_vec_tmp)

print "\n"

# print each component
for idx in range(dim) :
    print "print the %dth component... " % (idx+1)

    fig = figure()
    ax = gca() 

#   read each trajectory file
    f_idx = 0
    ub = 0 
    for f in files[:num_file] :
        traj_file = open("%s/%s" % (data_dir, f), 'r')
        itmp, = [ int(x) for x in traj_file.readline().split() ] 
        x = np.loadtxt(traj_file, usecols=(idx+1,), unpack=True)
        x = np.append(x, x[-1])
        ub = max(ub, max(x))
        traj_file.close()
        plot( t_vec[f_idx], x, linewidth=3, color=lc[f_idx], label="traj %d" % f_idx)
        f_idx = f_idx + 1

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    plt.xlim( 0, max_t )
#    plt.ylim( 0, ub+2)

    if idx==0 :
#      plt.ylim( 1e-2, ub+2)
      ax.set_yscale('log')
    else :
     plt.ylim( 0, ub+2)

    xlabel('$t$', fontsize=fs, labelpad=-5)
    title('$x^{(%d)}$' % (idx+1), fontsize=fs)
#    lg = legend(bbox_to_anchor=(-0.0, -0.05, 1.0, 1.0), prop={'size':18})
#    lg.draw_frame(False)
#    plt.gcf().subplots_adjust(bottom=0.2) 

    fig.tight_layout()
   
    savefig("../fig/ex2-x%d-traj.eps" % (idx+1) )
