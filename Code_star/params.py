# Physical quantities
G = 6.67E-11                # N m^2 / kg^2
Msun = 1.99E30              # Kg
Rsun = 6.957E5              # Km
tau = ((Rsun)**3 / (G*Msun))**0.5 # time unit

# General model quantities
n_dim = 2
n_step = 8000
n_grid = 50
dt = tau * 0.0015

# Star 1
n_star_1 = 700
M_1 = Msun
R_1 = 0.3*Rsun

