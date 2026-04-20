# Physical quantities
G = 6.67E-11                # N m^2 / kg^2
Msun = 1.99E30              # Kg
Rsun = 6.9E11                 # m
tau = (Rsun**3 / (G*Msun))**0.5 # time unit

# General model quantities
n_dim = 2
n_step = 1800
n_grid = 50
dt = tau * 0.0001

# Galaxy 1
n_star_1 = 1500
M_1 = Msun
R_1 = 0.3*Rsun

