import numpy as np
#import numba

def getAcceleration(p, v, mp, n_dim, n_particle, n_grid, G, Rstar, pos):
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
    h = 0.08*Rstar
    C = -5/(14 * np.pi * h**2)
    kappa = 1e-2
    visc = 2e-3

    #Calculate the acceleration  
        # P gradient
            #smoothing kernel
    dW = getdW(n_particle, n_dim, p, h, C, Rstar)

            #sum
    grad_P[0,:] = -2*kappa*sum(mp*dW[0,:])
    grad_P[1,:] = -2*kappa*sum(mp*dW[1,:])

        # Gravity
    for j in range(n_particle):
        r_ji_x = p[0,j] - p[0,:]
        r_ji_y = p[1,j] - p[1,:]

            #Clip anything too close to 0
        r_ji = np.clip(np.sqrt((r_ji_x)**2 + (r_ji_y)**2), 0.01*Rstar, None)

            #Get the gravitational acceleration
        grav[0,j] = - G * sum((mp[:] / (r_ji)**3)*r_ji_x)
        grav[1,j] = - G * sum((mp[:] / (r_ji)**3)*r_ji_y)
    
    #Get the acceleration
    acc[0,:] = -grad_P[0,:] - visc*v[0,:] + grav[0,:]
    acc[1,:] = -grad_P[1,:] - visc*v[1,:] + grav[1,:]

    #Get the density
    rho = np.zeros((n_grid, n_grid))
    W = getKernel(n_grid, n_particle, n_dim, p, h, C, pos, Rstar)
    
    rho = mp[0]*W
    return acc, rho


def leapfrog(v, dt, acc, p, mp, n_dim, n_particle, n_grid, G, Rstar, pos):
    """For a given set of particle, position, velocity, acceleration and mass, calculate the new 
    position after time dt from the velocity and acceleration"""

    #Update the parameters
    new_v = v + dt/2 * acc
    new_p = p + dt*new_v
    new_acc, rho = getAcceleration(new_p, new_v, mp, n_dim, n_particle, n_grid, G, Rstar, pos)
    new_new_v = new_v + dt/2*new_acc

    #Also update the density map
    #rho = getDensity(n_particle, n_dim, p, Rg)

    #Reset for the next loop
    p = new_p
    v = new_new_v
    acc = new_acc
    print(p/Rstar)
    return p, v, acc, rho




###############################################################


def getKernel(n_grid, n_particle, n_dim, p, h, C, pos, Rstar):
    """Get the smoothing kernel for particle j as a function of position, 
    radius to 0 and some constants h and C"""

    #Initialise
    W = np.zeros((n_particle, n_grid, n_grid), dtype = np.float32)
    dr = np.zeros((n_particle, n_dim, n_grid, n_grid), dtype = np.float64)
    r = np.zeros((n_particle, n_grid, n_grid), dtype = np.float32)

    # Get distances from all positions to all particles
    for j in range(n_particle):
        # We make a n_grid by n_grid matrix which is identical in all points in one dimension 
        # for one particle
        # Then we do the same in the other dimension
        # When we combine the two, and later on sum over the particles, these identical arrays combine
        # with other identical arrays in different way, producing an n_grid by n_grid matrice which should
        # balance the weight and produce the right density in all points
        dr[j,0,:,:] = np.broadcast_to(pos[0,:] - p[0,j], (n_grid, n_grid))
        dr[j,1,:,:] = np.broadcast_to(pos[1,:] - p[1,j], (n_grid, n_grid))

        r[j] = np.abs(dr[j,:,:,:]).sum(axis=0)
        r[j] = np.clip(r[j], 0.01*Rstar, None)

    q = r/h

    mask1 = (q<1)
    mask2 = (q>=1) & (q<2)

    W[mask1] = C*((2-q[mask1])**3 - 4*(1-q[mask1])**3)
    W[mask2] = C*((2-q[mask2])**3)
    W = np.sum(W, axis=0)
    print(W)
    return W



###############################################################

def getdW(n_particle, n_dim, p, h, C, Rstar):
    """Get the gradient of the smoothing kernel accross space for 
    all particles (n_particle) in every dimension (n_dim) as a function 
    of their position in space p and the constants h and C."""

    #Define the arrays first
    dW_temp = np.zeros((n_particle, n_dim, n_particle))

    dr = np.zeros((n_particle, n_dim, n_particle))
    r = np.zeros((n_particle, n_particle))
    q = np.zeros((n_particle, n_particle))

    #Get variables
    for j in range(n_particle):
        dr[j,0,:] = np.array([p[0,j]-p[0,:]])
        dr[j,1,:] = np.array([p[1,j]-p[1,:]])
        dr = np.clip(dr, 0.01*Rstar, None)

        r[j,:] = np.clip(np.sqrt((dr[j,:,:]**2).sum(axis=0)), 0.01*Rstar, None)
    q = r/h

    #Get dW
        #masks
    mask1 = (q<1)
    mask2 = (q>=1) & (q<=2) 
        #make a double one for the bigger arrays
    doublemask1 = np.broadcast_to(mask1[:,None,:], (n_particle, n_dim, n_particle))
    doublemask2 = np.broadcast_to(mask2[:,None,:], (n_particle, n_dim, n_particle))

        #make doubles for the q arrays and r arrays
    q = np.broadcast_to(q[:,None,:], (n_particle, n_dim, n_particle))
    r = np.broadcast_to(r[:,None,:], (n_particle, n_dim, n_particle))

        #get dW_temp, without summing over n_particle yet
    dwtemp_mask1_1 = ((2.0-q[doublemask1])**2)/r[doublemask1]
    dwtemp_mask1_2 = - 4*(1-q[doublemask1])**2/r[doublemask1]
    dW_temp[doublemask1] = dwtemp_mask1_1 + dwtemp_mask1_2

    dWtemp_mask2_1 = ((2.0-q[doublemask2])**2)
    dW_temp[doublemask2] = dWtemp_mask2_1/r[doublemask2]

        #then sum over them to get the right arrays
    dW = np.sum((3*C/h)*dW_temp*dr, axis=0)
    return dW

    