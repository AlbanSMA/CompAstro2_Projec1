import numpy as np

def mkGalaxy2D(n_particle, R, M, n_dim, G, x0=0.0, y0=0.0, vx0=0.0, vy0=0.0):
    """
    Here you should implement a routine that places particles in the XY-plane
    and assigns an initial velocity to each of them. It should return three
    arrays: the positions array ((n_dim, n_particle)), the velocity array ((n_dim, n_particle))
    and the particle mass array (n_particle), where in this case n_dim = 2.

    Input parameters:
        n_particle: number of particles to distribute
                 R: Radius of the galaxy (it is the sigma parameter of a random Gaussian distribution) [m]
                 M: Total mass of the galaxy [kg]
                x0: x-coordinate of the center of the galaxy [m]
                y0: y-coordinate of the center of the galaxy [m]
               vx0: linear velocity component in the x-plane [m/s]
               vy0: linear velocity component in the y-plane [m/s]
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
    M_r = np.array([i*(a+1) for a, i in enumerate(mp)])

    # get the velocity from that mass and radius r and make it two vectors v_x and v_y in the right directions
    v_temp = np.sqrt(G*M_r/abs(r))
    v[0,:] = - v_temp*np.sin(theta)*0
    v[1,:] = v_temp*np.cos(theta)*0


    #add the global offsets in position and velocity
    p[0,:] += x0
    p[1,:] += y0
    v[0,:] += vx0
    v[1,:] += vy0

    return p, v, mp


def InitialConditions(n_galaxy1, n_galaxy2, M_1, M_2, R_1, R_2, x02, y02, vx02, vy02, n_dim, G, mkGalaxy2d):
    """
    Here you should create the initial conditions. You should call mkGalaxy2D
    twice to create two different galaxies at different coordinates
    """

    # place galaxies
    p1, v1, mp1 = mkGalaxy2d(n_galaxy1, R_1, M_1, n_dim, G)
    p2, v2, mp2 = mkGalaxy2d(n_galaxy2, R_2, M_2, n_dim, G, x02, y02, vx02, vy02)


    # now concatenate those arrays into one large array
    p  = np.zeros((n_dim, n_galaxy1+n_galaxy2))
    v  = np.zeros((n_dim, n_galaxy1+n_galaxy2))
    mp = np.zeros(n_galaxy1+n_galaxy2)

    p[:, 0:n_galaxy1] = p1
    p[:, n_galaxy1::] = p2
    v[:, 0:n_galaxy1] = v1
    v[:, n_galaxy1::] = v2
    mp[0:n_galaxy1] = mp1
    mp[n_galaxy1::] = mp2

    return p, v, mp