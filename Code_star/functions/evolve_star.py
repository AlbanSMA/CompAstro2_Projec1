import numpy as np
#import numba

def getAcceleration(p, v, mp, n_dim, n_particle, G, Rg):
    """
    Given an array of particle positions p(n_dim, n_particle),
    calculate the net acceleration at the location of each particle
    from their gravitational interaction with the rest
    """

    n_dim, n_particle = p.shape
    acc = np.zeros((n_dim, n_particle))
    grad_P = np.zeros((n_dim, n_particle))
    grav = np.zeros((n_dim, n_particle))

    #get some constants
    h = 0.02*Rg
    C = 5/(14*np.pi*h**2)
    kappa = 1e-2
    visc = 2e-3

    #Calculate the acceleration  
        # P gradient
            #smoothing kernel
    dW = getdW(n_particle, n_dim, p, h, C)

            #sum
    grad_P[0,:] = -2*kappa*sum(mp*dW[0,:])
    grad_P[1,:] = -2*kappa*sum(mp*dW[1,:])

        # Gravity
    p_sub_0 = np.repeat(p[0,:], len(p[0,:]), axis=0).reshape(len(p[0,:]), 
                            len(p[0,:]))
    p_sub_1 = np.repeat(p[1,:], len(p[1,:]), axis=0).reshape(len(p[1,:]), 
                            len(p[1,:]))
    r_ji_x = p[0,:].transpose()-p_sub_0
    r_ji_y = p[1,:].transpose()-p_sub_1

            #Clip anything too close to 0
    r_ji = np.clip(np.sqrt((r_ji_x)**2 + (r_ji_y)**2), 0.1*Rg, None)

            #Get the acceleration
    grav[0,:] = - G * sum((mp[:] / (r_ji)**3)*r_ji_x)
    grav[1,:] = - G * sum((mp[:] / (r_ji)**3)*r_ji_y)
    
    #Get the acceleration
    acc[0,:] = -grad_P[0,:] - visc*v[0,:] + grav[0,:]
    acc[1,:] = -grad_P[1,:] - visc*v[1,:] + grav[1,:]

    print(acc, grad_P, visc*v, grav)
    #Also get the density
    W = getKernel(n_particle, n_dim, p, h, C)
    rho = sum(mp*W)

    return acc, rho


def leapfrog(v, dt, acc, p, mp, n_dim, n_particle, G, Rg):
    """For a given set of particle, position, velocity, acceleration and mass, calculate the new 
    position after time dt from the velocity and acceleration"""

    #Update the parameters
    new_v = v + dt/2 * acc
    new_p = p + dt*new_v
    new_acc, rho = getAcceleration(new_p, new_v, mp, n_dim, n_particle, G, Rg)
    new_new_v = new_v + dt/2*new_acc

    #Reset for the next loop

    p = new_p
    v = new_new_v
    acc = new_acc
    return p, v, acc, rho



def getKernel(n_particle, n_dim, p, h, C):
    """Get the smoothing kernel for particle j as a function of position, 
    radius to 0 and some constants h and C"""

    W = np.zeros((n_particle, n_dim, n_particle))
    p_sub_0 = np.repeat(p[0,:], len(p[0,:]), axis=0).reshape(len(p[0,:]), 
                            len(p[0,:]))
    p_sub_1 = np.repeat(p[1,:], len(p[1,:]), axis=0).reshape(len(p[1,:]), 
                            len(p[1,:]))
    r = np.concatenate((np.sqrt((p[0,:]-p_sub_0)**2), np.square((p[1,:]-p_sub_1)**2))).reshape(n_particle, n_dim, n_particle)
    q = abs(r)/h

    mask1 = (q>=0) & (q<1)
    mask2 = (q>=1) & (q<=2)

    W[mask1] = C*(2-q[mask1])**3 - 4*(1-q[mask1])**3
    W[mask2] = C*(2-q[mask2])**3
    #Above q=2, everything remains 0

    return W

def getdW(n_particle, n_dim, p, h, C):
    """Get the gradient of the smoothing kernel accross space for 
    all particles (n_particle) in every dimension (n_dim) as a function 
    of their position in space p and the constants h and C."""
    dW = np.zeros((n_particle, n_dim, n_particle))
    dr = np.zeros((n_particle*n_dim*n_particle))

    #Get variables
    p_sub_0 = np.repeat(p[0,:], len(p[0,:]), axis=0).reshape(len(p[0,:]), 
                            len(p[0,:]))
    p_sub_1 = np.repeat(p[1,:], len(p[1,:]), axis=0).reshape(len(p[1,:]), 
                            len(p[1,:]))
    dr = np.concatenate((p[0,:]-p_sub_0, p[1,:]-p_sub_1)).reshape(n_particle, n_dim, n_particle)
    r = np.sqrt((dr**2).sum())
    q = r/h

    #Get dW
    mask1 = (q>=0) & (q<1)
    mask2 = (q>=1) & (q<=2) 
    dW[mask1] = 3*(C/h)*((2-q)**2 - 4*(1-q)**2) * dr/r
    dW[mask2] = 3*(C/h)*((2-q)**2) * dr/r

    #dW should be of shape : dW per interactions with other particles * dW per dimension * dW per particles
    return dW

    