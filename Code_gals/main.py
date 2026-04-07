import os
import sys
import matplotlib.pyplot as plt; plt.ion()

plt.ion()
sys.path.append("C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Code_gals")

from params import *
from functions.make_gal import mkGalaxy2D, InitialConditions
from functions.evolve_gal import getAcceleration, leapfrog
from functions.plots import ini_plots, update_plots



if __name__ == "__main__": # main program, not executed if you import this file

    # Create output subfolder to store snapshots from the simulation
    path = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Plots\\Plots_gal\\"
    os.makedirs(path+'nbody_output', exist_ok=True)
    
    # Initialisation of the model
    p, v, mp = InitialConditions(n_galaxy_1, n_galaxy_2, M_1, M_2, R_1, R_2, x02, y02, vx02, vy02, n_dim, G, mkGalaxy2D) #positions, velocity, particle mass
    acc = getAcceleration(p, mp, n_dim, n_particle, G, Rg) #acceleration

    # Evolution of the model
    k = 0           #steps
    f, ax, d0, d1 = ini_plots(p, n_galaxy_1, Rg, path)   #initialise the plots

    for tt in range(1, n_step):
        # Implement your Leapfrog algorithm to take a step in v and p
        p, v, acc = leapfrog(v, dt, acc, p, mp, n_dim, n_particle, G, Rg)

        # Update the plots with the new particle positions
        # to save some time, only update the plot every 4 iterations
        if(tt%10 == 0):
            update_plots(d0, d1, p, k, n_galaxy_1, Rg, f, ax, path)
            k += 1