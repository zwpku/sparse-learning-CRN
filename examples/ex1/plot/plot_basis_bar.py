#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D   

data = np.loadtxt("../output/ratio_of_basis_in_ai.txt") 
channel_num = len(data[:,0])
basis_num = len(data[0,:])
# Data generation
c_vec = np.linspace(0, channel_num, channel_num)
basis_vec = np.linspace(0, basis_num, basis_num)

np.clip(data, 0.01, None, out=data)
print(data)
# Plotting
fig = plt.figure()
ax = fig.add_subplot(111)
width=0.45
pad = basis_num * width + 0.1

for idx in range(channel_num):
    c_vec[idx] = idx * pad

ax.bar(c_vec, data[:,0], width, color='r', align='edge', label=r'$1$')
ax.bar(c_vec+width, data[:,1], width, color='b', align='edge', label=r'$x^{(1)}$')
ax.bar(c_vec+2*width, data[:,2], width, color='g',align='edge', label=r'$x^{(2)}$')
ax.bar(c_vec+3*width, data[:,3], width, color='k', align='edge', label=r'$(x^{(1)})^2$')
ax.bar(c_vec+4*width, data[:,4], width, color='y', align='edge', label=r'$x^{(1)}x^{(2)}$')
ax.bar(c_vec+5*width, data[:,5], width, color='m', align='edge', label=r'$(x^{(2)})^2$')

plt.axvline(x=(c_vec[0]+6*width+c_vec[1])*0.5, color='gray',   linestyle='--')
plt.axvline(x=(c_vec[1]+6*width+c_vec[2])*0.5, color='gray',   linestyle='--')
plt.axvline(x=(c_vec[2]+6*width+c_vec[3])*0.5, color='gray',   linestyle='--')

lg = legend(bbox_to_anchor=(0.30, 0.00, 0.48, 1.02), prop={'size':18}, ncol=2)
lg.draw_frame(False)

#ax.set_xlabel('Channel')
ax.set_xticks(c_vec+width*3)
ax.set_xticklabels((r'Channel $1$', r'Channel $2$',r'Channel $3$',r'Channel $4$'))
plt.xlim(0, pad * channel_num)
plt.ylim(0.0, 1.43)
#ax.bar3d(Xi, Yi, Zi, dx, dy, dz, color = 'w')

savefig("../fig/bar2d-for-basis.eps" )


