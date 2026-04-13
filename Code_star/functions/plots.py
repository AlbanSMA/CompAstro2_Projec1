import matplotlib.pyplot as plt
import numpy as np

def ini_plots(p, rho, Rg, path):
    f, (ax1, ax2) = plt.subplots(1,2, figsize=(6,6))

    # Plot star
    d0, = ax1.plot(p[0,:]/Rg, p[1,:]/Rg, 'o', color='orangered', ms=1.7, alpha=0.35, linewidth=0, mew=0)
    rho0 = ax2.imshow(rho)

    # plot labels
    ax1.set_ylabel("y/Rg")
    ax1.set_xlabel("x/Rg")
    ax1.set_title("t_step={0}".format(0))

    ax2.set_ylabel("y")
    ax2.set_xlabel("x")
    ax2.set_title("t_step={0}".format(0))

    # set lim
    ax1.set_ylim()
    ax1.set_xlim()

    ax2.set_ylim()
    ax2.set_xlim()

    f.savefig(path+"nbody_output/img{0:05d}.png".format(0), dpi=300, format='png')
    return f, ax1, ax2, d0, rho0


def update_plots(d0, rho0, p, rho, k, Rg, f, ax1, ax2, path):
    # update the data points in each of the plots,
    # keeping all labels and axes exactly the same
    d0.set_data(p[0,:]/Rg, p[1,:]/Rg)
    rho0.set_data(rho)

    # Update title with the time step
    ax1.set_title(f"t_step={k}")
    ax2.set_title(f"t_step={k}")
            
    # Force re-drawing the figure    
    f.canvas.draw()
    f.canvas.flush_events()

    # save image from the simulation to your drive as png
    f.savefig(path+"nbody_output/img{0:05d}.png".format(k), dpi=300, format='png')