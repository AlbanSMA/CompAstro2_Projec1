import numpy as np

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
    for j in range(n_particle):
        #some useful functions
        r_j = np.sqrt(p[0,j]**2 + p[1,j]**2)

        #The acceleration is : a_j = - grad_P - visc*v + g
        # P gradient
            # Smoothing kernel
        r = np.sqrt((p[0,j]-p[0,:])**2 + (p[1,j]-p[1,:]))
        q = abs(r - r_j)/h
        W1 = [C*(2-q[np.where(q>=0 and q<1)])**3 - 4*(1-q[np.where(q>=0 and q<1)])**3]
        W2 = [C*(2-q[np.where(q>=1 and q<2)])**3]
        W = np.array(W1+W2)
            
            #sum
        grad_P[0,j] = -2*kappa*sum(mp*W)
        grad_P[1,j] = -2*kappa*sum(mp*W)

        # Gravity
        r_ji_x = p[0,j]-p[0,:]
        r_ji_y = p[1,j]-p[1,:]
        r_ji = np.clip(np.sqrt((r_ji_x)**2 + (r_ji_y)**2), 0.1*Rg, None)

        #Get the acceleration
        grav[0,j] = - G * sum((mp[:] / (r_ji)**3)*r_ji_x)
        grav[1,j] = - G * sum((mp[:] / (r_ji)**3)*r_ji_y)

    acc[0,:] = -grad_P[0,:] - visc*v[0,:] + grav[0,:]
    acc[1,:] = -grad_P[1,:] - visc*v[1,:] + grav[1,:]

    return acc


def leapfrog(v, dt, acc, p, mp, n_dim, n_particle, G, Rg):
    """For a given set of particle, position, velocity, acceleration and mass, calculate the new 
    position after time dt from the velocity and acceleration"""

    #Update the parameters
    new_v = v + dt/2 * acc
    new_p = p + dt*new_v
    new_acc = getAcceleration(new_p, new_v, mp, n_dim, n_particle, G, Rg)
    new_new_v = new_v + dt/2*new_acc

    #Reset for the next loop
    p = new_p
    v = new_new_v
    acc = new_acc
    print(acc)
    return p, v, acc