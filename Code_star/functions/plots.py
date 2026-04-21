import matplotlib.pyplot as plt

def ini_plots(p, rho, Rstar, path):
    f, (ax1, ax2) = plt.subplots(1,2, figsize=(14,6))

    # Plot star
    d0, = ax1.plot(p[0,:]/Rstar, p[1,:]/Rstar, 'o', color='orangered', ms=1.7, alpha=0.35, linewidth=0, mew=0)
    d1 = ax2.imshow(rho)

    # plot labels
    ax1.set_ylabel("y/Rstar")
    ax1.set_xlabel("x/Rstar")
    ax1.set_title("t_step={0}".format(0))

    ax2.set_ylabel("y")
    ax2.set_xlabel("x")
    ax2.set_title("t_step={0}".format(0))

    # set lim
    ax1.set_ylim()
    ax1.set_xlim()

    ax2.set_ylim()
    ax2.set_xlim()

    ax1.grid(alpha=0.5)

    f.savefig(path+"nbody_output/img{0:05d}.png".format(0), dpi=300, format='png')
    return f, ax1, ax2, d0, d1


def update_plots(d0, d1, p, rho, k, Rstar, f, ax1, ax2, path):
    # update the data points in each of the plots,
    # keeping all labels and axes exactly the same
    d0.set_data(p[0,:]/Rstar, p[1,:]/Rstar)
    d1.set_data(rho)

    # Update title with the time step
    ax1.set_title(f"t_step={k}")
    ax2.set_title(f"t_step={k}")

    ax1.set_ylim()
    ax1.set_xlim()
            
    # Force re-drawing the figure    
    f.canvas.draw()
    f.canvas.flush_events()

    # save image from the simulation to your drive as png
    f.savefig(path+"nbody_output/img{0:05d}.png".format(k), dpi=300, format='png')