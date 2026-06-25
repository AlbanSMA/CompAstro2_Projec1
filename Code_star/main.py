import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt; plt.ion()

plt.ion()
sys.path.append("C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Code_star")

from params import *
from functions.make_star import mkStarD, getDensitypos
from functions.evolve_star import getAcceleration, leapfrog
from functions.plots import ini_plots, update_plots


if __name__ == "__main__": # main program, not executed if you import this file
    timestart = time.time()

    # Create output subfolder to store snapshots from the simulation, make sure it is empty
    path = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project1\\Plots\\Plots_star\\"
    folder = "Msun0.3Rsun"
    if os.path.exists(path+folder) == True:
        files = os.listdir(path+folder)
        [os.remove(path+folder+"\\"+files[i]) for i in range(len(files))]
    
    os.makedirs(path+folder, exist_ok=True)

    allrho = np.zeros((n_step, int(n_grid)))
    totacc = np.zeros((n_step))
    totgrav = np.zeros((n_step))
    totDP = np.zeros((n_step))
    totvisc = np.zeros((n_step))
    
    # Initialisation of the model
    p, v, mp = mkStarD(n_star_1, R_1, M_1, n_dim) #positions, velocity, particle mass

    
    pos = getDensitypos(n_grid, n_dim, R_1)         #positions at which to evaluate the density
    acc, rho, grav, grad_P, visc = getAcceleration(p, v, mp, n_dim, n_star_1, n_grid, G, R_1, pos) #acceleration

    allrho[0] = rho[int(n_grid/2)]


    # Evolution of the model
    k = 0           #steps
    f, ax1, ax2, d0, d1 = ini_plots(p, rho, R_1, path+folder, n_grid)   #initialise the plots

    for tt in range(1, n_step):
        try:
            # Implement your Leapfrog algorithm to take a step in v and p
            p, v, acc, rho, grav, grad_P, visc = leapfrog(v, dt, acc, p, mp, n_dim, n_star_1, n_grid, G, R_1, pos)

            # Update the plots with the new particle positions
            # to save some time, only update the plot every 4 iterations
            if (tt%50 == 0):
                totacc[k] = np.mean(acc)
                totgrav[k] = np.mean(grav)
                totDP[k] = np.mean(grad_P)
                totvisc[k] = np.mean(visc)
                allrho[k] = rho[int(n_grid/2)]
                update_plots(d0, d1, p, rho, tt, R_1, f, ax1, ax2, path+folder)
                k += 1
        except KeyboardInterrupt:
            break
    
    # Make some additional plots
    newfig, newax = plt.subplots(1,1, figsize = (20, 5))

    img = newax.imshow(allrho[:k].T, cmap="magma")
    newax.set_xlabel("Time")
    newax.set_ylabel(r"Density")
    newax.set_title("Evolution of the density as a function of time")
    newax.grid(alpha=0.5)
    newfig.colorbar(img, ax=newax)
    newfig.savefig(path+folder+"\\DensityVtime.pdf", dpi=300, bbox_inches = "tight")

    tt = int(tt/50)
    newfig2, (newax2, newax3) = plt.subplots(2,1)
    t = np.linspace(0, tt, tt)

    newax2.plot(t[totacc[:tt] != 0], totacc[:tt][totacc[:tt] != 0]/totacc[:tt][totacc[:tt] != 0], color = "black", label="Total acceleration")
    newax2.plot(t[totacc[:tt] != 0], totgrav[:tt][totacc[:tt] != 0]/totacc[:tt][totacc[:tt] != 0], color="red", label="Gravity")
    newax2.plot(t[totacc[:tt] != 0], totDP[:tt][totacc[:tt] != 0]/totacc[:tt][totacc[:tt] != 0], color="green", label="Pressure gradient")
    newax2.plot(t[totacc[:tt] != 0], totvisc[:tt][totacc[:tt] != 0]/totacc[:tt][totacc[:tt] != 0], color="blue", label="Viscosity")
    newax2.legend()
    newax2.grid(alpha=0.5)
    newax2.set_xlabel("Time (steps)")
    newax2.set_ylabel("Acceleration")
    newax2.set_title ("Mean acceleration of the particles over time")

    newax3.plot(t, totacc[:tt], color = "black", linewidth=3, label="Total acceleration")
    newax3.plot(t, totgrav[:tt], color="red", label="Gravity")
    newax3.plot(t, totDP[:tt], color="green", label="Pressure gradient")
    newax3.plot(t, totvisc[:tt], color="blue", label="Viscosity")
    newax3.legend()
    newax3.grid(alpha=0.5)
    newax3.set_xlabel("Time (steps)")
    newax3.set_ylabel("Acceleration")

    newfig2.savefig(path+folder+"\\AccVtime.pdf", dpi=300, bbox_inches="tight")

    timestop = time.time()

    timing = timestop-timestart
    print(timing/60)

