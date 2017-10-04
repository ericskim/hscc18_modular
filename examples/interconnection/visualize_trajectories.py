import sys

import numpy as np
import matplotlib.pyplot as plt

filename = "traj.txt"
pngname = "coop_reach"
dim = 6

def get_traj_array(filename, s_or_i = 'state'):
    traj = []
    inputs = []
    with open(filename) as f:
        l = 0 
        for line in f:
            x = [i for i in line.rstrip(' \n').split(' ')]
            if 'State:' in x:
                x = [float(i) for i in x[1:]]
                traj.append(x)
            if 'Input:' in x:
                x = [float(i) for i in x[1:]]
                inputs.append(x)
            # cur_pre_iter = np.rint(x[0])
    traj = np.array(traj)
    inputs = np.array(inputs)
    if s_or_i == 'state':
        return traj
    elif s_or_i == 'input':
        return inputs
T= 19
traj = get_traj_array("traj_diff_inits.txt", 'state')[0:T, :]
inputs = get_traj_array("traj_diff_inits.txt", 'input')[0:T, :]

div_traj = get_traj_array("traj_diff_inits_diverge.txt")[0:T,:]


fig = plt.figure()
# ax = fig.add_subplot(111)
ax2 = fig.add_subplot(111)

t= np.arange(T)

lw = 1
# ax.plot(t, traj[:,0],
#         t, traj[:,1], 
#         t, traj[:,2], 
#         t, traj[:,3], 
#         t, traj[:,4], 
#         t, traj[:,5], linewidth = 2*lw)
# ax.plot(t, div_traj[:,0], '--',
#         t, div_traj[:,1], '--', 
#         t, div_traj[:,2], '--', 
#         t, div_traj[:,3], '--', 
#         t, div_traj[:,4], '--', 
#         t, div_traj[:,5], '--', linewidth = lw, alpha = .4)

# td =.06
d = .07
# ax2.plot(t-2.5*td, inputs[:,0]-2.5*d,
#          t-1.5*td, inputs[:,1]-1.5*d, 
#          t-.5*td, inputs[:,2]-.5*d, 
#          t+.5*td, inputs[:,3]+.5*d, 
#          t+1.5*td, inputs[:,4]+1.5*d, 
#          t+2.5*td, inputs[:,5]+2.5*d, linestyle = '-', linewidth = 2*lw, drawstyle='steps-post', marker='o') 
ax2.set_xlim(-td-.1,T-1)
ax2.set_ylim(-2.25,2.25)

err = [np.zeros(T), np.ones(T)]
ax2.errorbar(t-2.5*td, inputs[:,0]-2.5*d,xerr=err, fmt='o',linewidth = 2*lw)
ax2.errorbar(t-1.5*td, inputs[:,1]-1.5*d,xerr=err, fmt='o',linewidth = 2*lw)
ax2.errorbar(t-0.5*td, inputs[:,2]-0.5*d,xerr=err, fmt='o',linewidth = 2*lw)
ax2.errorbar(t+0.5*td, inputs[:,3]+0.5*d,xerr=err, fmt='o',linewidth = 2*lw)
ax2.errorbar(t+1.5*td, inputs[:,4]+1.5*d,xerr=err, fmt='o',linewidth = 2*lw)
ax2.errorbar(t+2.5*td, inputs[:,5]+2.5*d,xerr=err, fmt='o',linewidth = 2*lw)


# ax2.plot(t, np.zeros(T), '--',
#         t, np.zeros(T), '--',
#         t, np.zeros(T), '--',
#         t, np.zeros(T), '--',
#         t, np.zeros(T), '--',
#         t, np.zeros(T), linewidth = lw, alpha = .4, drawstyle='steps')


# ax.set_xlim(0,T-1)
# ax.set_ylim(0,31)
plt.show()

