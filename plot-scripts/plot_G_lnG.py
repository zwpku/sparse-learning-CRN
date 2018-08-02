#!/usr/bin/python

"""

  This code plots the functions G(x), G'(x), (ln G), and (ln G)' to files.

"""

from pylab import *
from numpy import *

lc = ['b', 'r', 'k', 'c', 'm', 'y']

# activation function: max(x, 0)
def activation(x) :
  if x >= 0:
      return x 
  else:
      return 0.0 

# G(x) = delta * ln(1+exp(x/delta))
def Geps(x) :
  if x > g_cut * delta :
      return x 
  if x < -g_cut * delta :
      return 0.0 
  return delta * log( 1 + exp( x / delta ) ) 

# G'(x) : the derivative of the function G(x)
def d_Geps(x) :
  if x > g_cut * delta :
      return 1.0 
  if x < -g_cut * delta :
      return 0.0 
  return 1.0 / (1 + exp(-x / delta)) 

# (ln G)'(x)
def d_logGeps(x) :
  if x > g_cut * delta :
      return 1.0 / x 
  if x < -g_cut * delta :
      return 1.0 / delta 
  return 1.0 / (delta * (1 + exp(-x / delta)) * log(1 + exp(x / delta) ) ) 

# range of x, [-xb, xb] 
xb = 3

# number of discrete interval
n = 200

g_cut = 30 
    
xvec = np.linspace( -xb, xb, n, endpoint=True )

# max(x,0)
y_vec_0 = [ activation(x) for x in xvec ]   

# G(x) and G'(x), delta = 0.4
delta = 0.4 
y_vec = [ Geps(x) for x in xvec ]   
dy_vec = [ d_Geps(x) for x in xvec ]   

# G(x) and G'(x), delta = 1.0
delta = 1.0 
y_vec_1 = [ Geps(x) for x in xvec ]   
dy_vec_1 = [ d_Geps(x) for x in xvec ]   

ax = gca() 
plot( xvec, y_vec_0, linestyle='-', linewidth=3, color=lc[0], label='$\epsilon=0$' ) 
plot( xvec, y_vec, linestyle=':', linewidth=3, color=lc[1], label='$\epsilon=0.4$' ) 
plot( xvec, y_vec_1, linestyle='-', linewidth=3, color='k', label='$\epsilon=1$' ) 

ax.set_aspect('auto')
xticks(np.arange(-xb, xb+1, step=1), fontsize=25)
yticks(fontsize=24)
legend(bbox_to_anchor=(0.01,0.89), loc="upper left", frameon=False, fontsize=23)

plt.xlim(-xb, xb)
plt.gcf().subplots_adjust(bottom=0.18) 
fs = 30
xlabel('x', fontsize=fs, labelpad=-5)
ax.set_title('$G_\epsilon(x)$', fontsize=fs)
ax.title.set_position([.5, 1.01])

#plt.tight_layout()
fig_file_name = '../fig/geps.eps' 
savefig(fig_file_name)

plt.clf()
ax= gca()
plot( xvec, dy_vec, linestyle=':', linewidth=3, color=lc[1], label='$\epsilon=0.4$' ) 
plot( xvec, dy_vec_1, linestyle='-', linewidth=3, color='k', label='$\epsilon=1$' ) 

ax.set_aspect('auto')
xticks(np.arange(-xb, xb+1, step=1), fontsize=25)
yticks(fontsize=24)
legend(bbox_to_anchor=(0.01,0.89), loc="upper left", frameon=False, fontsize=23)

plt.xlim(-xb, xb)
plt.gcf().subplots_adjust(bottom=0.18) 
fs = 30
xlabel('x', fontsize=fs, labelpad=-5)
title('$G_\epsilon\'(x)$', fontsize=fs)
ax.title.set_position([.5, 1.02])

#plt.tight_layout()
fig_file_name = '../fig/d_geps.eps' 
savefig(fig_file_name)

# (ln G)' for different values delta
xvec_pos = np.linspace( 0.01, xb, n, endpoint=True )
delta = 0.4 
d_logy_vec = [ d_logGeps(x) for x in xvec ]
delta = 0.6 
d_logy_vec_1 = [ d_logGeps(x) for x in xvec ]
delta = 1.0 
d_logy_vec_2 = [ d_logGeps(x) for x in xvec ]

d_logy_vec_3 = [ (1.0 / x)  for x in xvec_pos ]

plt.clf()
ax= gca()
plot( xvec, d_logy_vec, linestyle=':', linewidth=3, color=lc[1], label='$\epsilon=0.4$' ) 
plot( xvec, d_logy_vec_1, linestyle='-', linewidth=3, color=lc[2], label='$\epsilon=0.6$' ) 
plot( xvec, d_logy_vec_2, linestyle='-', linewidth=3, color=lc[3], label='$\epsilon=1$' ) 
plot( xvec_pos, d_logy_vec_3, linestyle='-', linewidth=3, color=lc[4], label='$1/x$' ) 

ax.set_aspect('auto')
xticks(np.arange(-xb, xb+1, step=1), fontsize=25)
yticks(fontsize=24)
legend(bbox_to_anchor=(0.01,0.89), loc="upper left", frameon=False, fontsize=23)

plt.xlim(-xb, xb)
plt.ylim(-0.2, 6)
plt.gcf().subplots_adjust(bottom=0.18) 
fs = 30
xlabel('x', fontsize=fs, labelpad=-5)
title('$G_\epsilon\'/G_\epsilon$', fontsize=fs)
ax.title.set_position([.5, 1.02])

#plt.tight_layout()
fig_file_name = '../fig/d_logg_eps.eps' 
savefig(fig_file_name)
