import os
import sys
import shutil
import matplotlib.pyplot as plt; plt.ion()

plt.ion()
sys.path.append("C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Code_star")

from params import *
from functions.make_star import mkStarD, getDensitypos
from functions.evolve_star import getAcceleration, leapfrog
from functions.plots import ini_plots, update_plots


if __name__ == "__main__": # main program, not executed if you import this file

    # Create output subfolder to store snapshots from the simulation, make sure it is empty
    path = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Plots\\Plots_star\\"
    
    if os.path.exists(path+"nbody_output") == True:
        files = os.listdir(path+"nbody_output")
        [os.remove(path+"nbody_output\\"+files[i]) for i in range(len(files))]
    
    os.makedirs(path+'nbody_output', exist_ok=True)
    
    # Initialisation of the model
    p, v, mp = mkStarD(n_star_1, R_1, M_1, n_dim, G) #positions, velocity, particle mass

    
    pos = getDensitypos(n_grid, n_dim, R_1)         #positions at which to evaluate the density
    acc, rho = getAcceleration(p, v, mp, n_dim, n_star_1, n_grid, G, R_1, pos) #acceleration


    # Evolution of the model
    k = 0           #steps
    f, ax1, ax2, d0, d1 = ini_plots(p, rho, R_1, path)   #initialise the plots

    for tt in range(1, n_step):
        # Implement your Leapfrog algorithm to take a step in v and p
        p, v, acc, rho = leapfrog(v, dt, acc, p, mp, n_dim, n_star_1, n_grid, G, R_1, pos)

        # Update the plots with the new particle positions
        # to save some time, only update the plot every 4 iterations
        if(tt%4 == 0):
            update_plots(d0, d1, p, rho, tt, R_1, f, ax1, ax2, path)
            k += 1