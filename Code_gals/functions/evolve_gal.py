import numpy as np

def getAcceleration(p, mp, n_dim, n_particle, G, Rg):
    """
    Given an array of particle positions p(n_dim, n_particle),
    calculate the net acceleration at the location of each particle
    from their gravitational interaction with the rest
    """

    n_dim, n_particle = p.shape
    acc = np.zeros((n_dim, n_particle))
    
    # Acceleration equation
    # F_j = m_j * sum(i != j) G * m_i/r^3_ji = m_j * a_j
    # So a_j = - sum(i != j) G * m_i/r^3_ji

    p_sub_0 = np.repeat(p[0,:], len(p[0,:]), axis=0).reshape(len(p[0,:]), 
                            len(p[0,:]))
    p_sub_1 = np.repeat(p[1,:], len(p[1,:]), axis=0).reshape(len(p[1,:]), 
                            len(p[1,:]))
    r_ji_x = p[0,:].transpose()-p_sub_0
    r_ji_y = p[1,:].transpose()-p_sub_1

    #Clip anything too close to 0
    r_ji = np.clip(np.sqrt((r_ji_x)**2 + (r_ji_y)**2), 0.1*Rg, None)

    #Get the acceleration
    acc[0,:] = - G * sum((mp[:] / (r_ji)**3)*r_ji_x)
    acc[1,:] = - G * sum((mp[:] / (r_ji)**3)*r_ji_y)
    
    return acc


def leapfrog(v, dt, acc, p, mp, n_dim, n_particle, G, Rg):
    """For a given set of particle, position, velocity, acceleration and mass, calculate the new 
    position after time dt from the velocity and acceleration"""

    #Update the parameters
    new_v = v + dt/2 * acc
    new_p = p + dt*new_v
    new_acc = getAcceleration(new_p, mp, n_dim, n_particle, G, Rg)
    new_new_v = new_v + dt/2*new_acc

    #Reset for the next loop
    p = new_p
    v = new_new_v
    acc = new_acc
    return p, v, acc