# Physical quantities
G = 6.67E-11                # N m^2 / kg^2
Msun = 1.99E30              # Kg
Rsun = 6.9E11                # m
tau = ((0.5*Rsun)**3 / (G*Msun))**0.5 # time unit

# General model quantities
n_dim = 2
n_step = 4000
n_grid = 100
dt = tau * 0.0001

# Star 1
n_star_1 = 2000
M_1 = Msun
R_1 = 0.3*Rsun

