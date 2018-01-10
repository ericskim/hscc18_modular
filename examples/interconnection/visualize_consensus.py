import sys

import numpy as np
import matplotlib.pyplot as plt

T= 19

def get_traj_array(filename, s_or_i = 'state'):
    """
    Parses output of consensus.cc script.

    Toggle s_or_i to 'state' or 'input' to get inputs or states
    """
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


t= np.arange(T)

# Visualize state trajectories 
statefig = plt.figure()
traj = get_traj_array("traj_active.txt", 'state')[0:T, :]
div_traj = get_traj_array("traj_passive.txt", 'state')[0:T,:]
ax = statefig.add_subplot(111)
lw = 1
ax.plot(t, traj[:,0],
        t, traj[:,1], 
        t, traj[:,2], 
        t, traj[:,3], 
        t, traj[:,4], 
        t, traj[:,5], linewidth = 2*lw)
ax.plot(t, div_traj[:,0], '--',
        t, div_traj[:,1], '--', 
        t, div_traj[:,2], '--', 
        t, div_traj[:,3], '--', 
        t, div_traj[:,4], '--', 
        t, div_traj[:,5], '--', linewidth = lw, alpha = .4)
ax.set_title("State Trajectories")

# Visualize inputs
inputs = get_traj_array("traj_active.txt", 'input')[0:T, :]
inputfig = plt.figure()
ax2 = inputfig.add_subplot(111)
td =.06 # temporal offset for visualization
d = .07 # vertical offset for visualization
ax2.plot(t-2.5*td, inputs[:,0]-2.5*d,
         t-1.5*td, inputs[:,1]-1.5*d, 
         t-.5*td, inputs[:,2]-.5*d, 
         t+.5*td, inputs[:,3]+.5*d, 
         t+1.5*td, inputs[:,4]+1.5*d, 
         t+2.5*td, inputs[:,5]+2.5*d, linestyle = '-', linewidth = 2*lw, drawstyle='steps-post', marker='o') 
ax2.set_xlim(-td-.1,T-1)
ax2.set_ylim(-2.25,2.25)
ax2.set_title("Input Trajectories")

plt.show()

