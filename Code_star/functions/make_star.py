import numpy as np

def mkStarD(n_particle, R, M, n_dim, G):
    """
    Here you should implement a routine that places particles in the XY-plane
    and assigns an initial velocity to each of them. It should return three
    arrays: the positions array ((n_dim, n_particle)), the velocity array ((n_dim, n_particle))
    and the particle mass array (n_particle), where in this case n_dim = 2.

    Input parameters:
        n_particle: number of particles to distribute
        R: Radius of the star (it is the sigma parameter of a random Gaussian distribution) [m]
        M: Total mass of the star [kg]
    """

    #initialise the arrays
    p  = np.zeros((n_dim, n_particle))
    v  = np.zeros((n_dim, n_particle))
    mp = np.zeros(n_particle)

    #Position
        #create a gaussian spread of particle around the center at random angles
    r = np.random.normal(0, 0.34*R, n_particle)
    theta = np.random.uniform(0, 2*np.pi, n_particle)

        #change the coordinates
    p[0,:] = r*np.cos(theta)
    p[1,:] = r*np.sin(theta)

    #Mass particles
    mp[:] = M/n_particle

    #Velocity
    #calculate the mass inside radius r
    #M_r = np.array([i*(a+1) for a, i in enumerate(mp)])

    #v_temp = np.sqrt(G*M_r/abs(r))
    #v[0,:] = v_temp*np.cos(theta)
    #v[1,:] = v_temp*np.sin(theta)
    v[:] = 0
    return p, v, mp