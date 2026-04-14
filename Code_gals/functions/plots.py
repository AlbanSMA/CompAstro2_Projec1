import matplotlib.pyplot as plt
import numpy as np

def ini_plots(p, n_galaxy_1, Rg, path):
    f, ax = plt.subplots(figsize=(6,6))

    # Plot galaxy 1, all particles from 0:n_galaxu_1
    d0, = ax.plot(p[0,0:n_galaxy_1]/Rg, p[1,0:n_galaxy_1]/Rg, 'o', color='orangered', ms=1.7, alpha=0.35, linewidth=0, mew=0)

    # Plot galaxy 2 using a different color: particles from n_galaxy_1:(n_galaxu_1+n_galaxy2)
    d1, = ax.plot(p[0,n_galaxy_1::]/Rg, p[1,n_galaxy_1::]/Rg, 'o', color='dodgerblue', ms=1.7, alpha=0.35, linewidth=0, mew=0)

    
    # plot labels
    ax.set_ylabel("y/Rg")
    ax.set_xlabel("x/Rg")
    ax.set_title("t_step={0}".format(0))

    # set lim
    ax.set_ylim(-20, 20)
    ax.set_xlim(-20, 20)

    f.savefig(path+"nbody_output_offset/img{0:05d}.png".format(0), dpi=300, format='png')
    return f, ax, d0, d1


def update_plots(d0, d1, p, k, n_galaxy_1, Rg, f, ax, path):
    # update the data points in each of the plots,
    # keeping all labels and axes exactly the same
    d0.set_data(p[0,0:n_galaxy_1]/Rg, p[1,0:n_galaxy_1]/Rg)
    d1.set_data(p[0,n_galaxy_1::]/Rg, p[1,n_galaxy_1::]/Rg)

    # Update title with the time step
    ax.set_title(f"t_step={k}")
            
    # Force re-drawing the figure    
    f.canvas.draw()
    f.canvas.flush_events()

    # save image from the simulation to your drive as png
    f.savefig(path+"nbody_output_offset/img{0:05d}.png".format(k), dpi=300, format='png')
    return d0, d1