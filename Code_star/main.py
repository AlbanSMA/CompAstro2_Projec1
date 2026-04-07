import os
import matplotlib.pyplot as plt; plt.ion()

from Code_star.params import *
from Code_star.functions.make_gal import mkGalaxy2D
from Code_star.functions.evolve_gal import getAcceleration, leapfrog
from Code_star.functions.plots import ini_plots, update_plots

plt.ion()


if __name__ == "__main__": # main program, not executed if you import this file

    # Create output subfolder to store snapshots from the simulation
    path = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Plots\\Plots_star\\"
    os.makedirs(path+'nbody_output', exist_ok=True)
    
    # Initialisation of the model
    p, v, mp = mkGalaxy2D(n_star_1, M_1, R_1, n_dim) #positions, velocity, particle mass
    acc = getAcceleration(p, mp, n_dim, n_star_1, G, Rsun) #acceleration

    # Evolution of the model
    k = 0           #steps
    f, ax, d0 = ini_plots(p, Rsun, path)   #initialise the plots

    for tt in range(1, n_step):
        # Implement your Leapfrog algorithm to take a step in v and p
        p, v, acc = leapfrog(v, dt, acc, p, mp, n_dim, n_star_1, G, Rsun)

        # Update the plots with the new particle positions
        # to save some time, only update the plot every 4 iterations
        if(tt%10 == 0):
            update_plots(d0, p, k, Rsun, f, ax, path)
            k += 1