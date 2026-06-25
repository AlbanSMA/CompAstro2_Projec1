import numpy as np
import scipy as sp
#import numba

def getAcceleration(p, v, mp, n_dim, n_particle, n_grid, G, Rstar, pos):
    """
    Given an array of particle positions p(n_dim, n_particle),
    calculate the net acceleration at the location of each particle
    from their gravitational interaction with the rest
    """

    acc = np.zeros((n_dim, n_particle))
    grad_P = np.zeros((n_dim, n_particle))
    grav = np.zeros((n_dim, n_particle))
    rho = np.zeros((n_particle, n_grid, n_grid))

    #get some constants
    h = 0.02*Rstar
    C = 5/(14 * np.pi * h**2)
    kappa = 0.5e-2
    visc = 2e-3

    #Calculate the acceleration  
    for j in range(n_particle):
        mask = (np.arange(n_particle) != j)

        # Gravity
        r_ji_x = p[0,mask] - p[0,j]
        r_ji_y = p[1,mask] - p[1,j]

        r_ji = np.clip(np.sqrt((r_ji_x)**2 + (r_ji_y)**2), 0.1*Rstar, None)


            #Get the gravitational acceleration
        grav[0,j] = G * sum((mp[mask] / (r_ji)**3)*r_ji_x)
        grav[1,j] = G * sum((mp[mask] / (r_ji)**3)*r_ji_y)

        # P gradient
            #smoothing kernel
        dW = getdW(n_particle, n_dim, j, p, h, C)

            #sum
        grad_P[0,j] = 2*kappa*np.sum(mp[:]*dW[0,:], axis=0)
        grad_P[1,j] = 2*kappa*np.sum(mp[:]*dW[1,:], axis=0)

            
            # Density around j
        W = getKernel(n_grid, n_dim, j, p, h, C, pos)     
        rho[j,:,:] = mp[j]*W

    
    #Sum rho over all particles:
    rhosum = np.sum(rho, axis=0)

    #Get the acceleration
    acc[0,:] = grav[0,:] + grad_P[0,:] - visc*v[0,:]
    acc[1,:] = grav[1,:] + grad_P[1,:] - visc*v[1,:]

    
    print(f"acc:{np.mean(abs(acc)):.2e}", f"grav:{np.mean(abs(grav)):.2e}", 
          f"grad_P:{np.mean(abs(grad_P)):.2e}", f"visc:{np.mean(abs(visc*v)):.2e}")
    return acc, rhosum, grav, grad_P, visc*v


def leapfrog(v, dt, acc, p, mp, n_dim, n_particle, n_grid, G, Rstar, pos):
    """For a given set of particle, position, velocity, acceleration and mass, calculate the new 
    position after time dt from the velocity and acceleration"""

    #Update the parameters
    new_v = v + dt/2 * acc
    new_p = p + dt*new_v
    new_acc, rho, grav, grad_P, visc = getAcceleration(new_p, new_v, mp, n_dim, n_particle, n_grid, G, Rstar, pos)
    new_new_v = new_v + dt/2*new_acc

    #Reset for the next loop
    p = new_p
    v = new_new_v
    acc = new_acc
    print(f"v:{np.mean(abs(v)):.2e}")
    return p, v, acc, rho, grav, grad_P, visc




###############################################################


def getKernel(n_grid, n_dim, j, p, h, C, pos):
    """Get the smoothing kernel for particle j as a function of position, 
    radius to 0 and some constants h and C"""

    #Initialise
    W = np.zeros((n_grid, n_grid))
    dr = np.zeros((n_dim, n_grid, n_grid))
    r = np.zeros((n_grid, n_grid))

    # Get distances from all positions to the particle
    dr[0,:,:] = pos[0,:] - p[0,j]
    dr[1,:,:] = pos[1,:] - p[1,j]

    r[:,:] = np.sqrt((dr[:,:,:]**2).sum(axis=0))

    # Smooth
    q = r/h

    #Mask
    mask1 = (q<1)
    mask2 = (q>=1) & (q<2)

    #Get the kernel
    W[mask1] = C*((2-q[mask1])**3 - 4*(1-q[mask1])**3)
    W[mask2] = C*((2-q[mask2])**3)
    return W



###############################################################

def getdW(n_particle, n_dim, j, p, h, C):
    """Get the gradient of the smoothing kernel accross space for 
    all particles (n_particle) in every dimension (n_dim) as a function 
    of their position in space p and the constants h and C."""

    #Define the arrays first
    dW_temp = np.zeros((n_dim, n_particle))
    dr = np.zeros((n_dim, n_particle))
    r = np.zeros((n_particle))
    q = np.zeros((n_particle))

    #Get the distance to all other particles
    dr[0,:] = np.array([p[0,:]-p[0,j]])
    dr[1,:] = np.array([p[1,:]-p[1,j]])

    r[:] = np.sqrt((dr**2).sum(axis=0))

    q = r/h

    #Get dW
        #masks
    mask1 = (q>0) & (q<1)
    mask2 = (q>=1) & (q<2) 

        #get dW_temp
    dW_temp[:,mask1] = ((2.0-q[mask1])**2 - 4*(1-q[mask1])**2)/r[mask1]
    dW_temp[:,mask2] = ((2.0-q[mask2])**2)/r[mask2]

        #then multiply by dr and (-3*C/h)
    dW = (-3*C/h)*dW_temp*dr.reshape(n_dim,n_particle)
    return dW