#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from pylab import *
from numpy import *
from os import listdir
from os import path
import re

data_dir = "../output" 
max_num_file = 5
lc = ['b', 'k', 'g', 'r', 'y']
# font size
fs = 30

fig = figure()
ax = gca() 


average_weights = np.loadtxt("%s/ratio_of_basis_in_ai.txt" % data_dir)
channel_num = size(average_weights[:,0])
basis_num = size(average_weights[0,:])
print( "channel_num=%d, basis_num=%d\n" % (channel_num, basis_num) )

x=range(channel_num)
for i in range(basis_num):
    print(x)
    print(average_weights[:,i])
    if i==0:
        plt.bar(x, height=average_weights[:,i], width=0.45, alpha=0.8)
    if i > 0:
        plt.bar(x, height=average_weights[:,i], width=0.45, alpha=0.8, bottom=average_weights[:,i-1])

ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

#plt.xlim( 0, max_t )
#plt.ylim( -1, ub+3)

xlabel('$t$', fontsize=fs, labelpad=-5)
#title('$x^{(%d)}$' % (idx+1), fontsize=fs)
#lg = legend(bbox_to_anchor=(-0.0, -0.01, pos[idx], 1.0), prop={'size':18})
#lg.draw_frame(False)
#plt.gcf().subplots_adjust(bottom=0.2) 

fig.tight_layout()

savefig("../fig/weights-of-basis.eps" )
