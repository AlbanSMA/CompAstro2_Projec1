# Physical quantities
G = 6.67E-11                # N m^2 / kg^2
Msun = 1.99E30              # Kg
Mg = 1.54E12 * Msun         # kg
Rg = 9.5E17                 # Km
tau = (Rg**3 / (G*Mg))**0.5 # time unit

# General model quantities
n_dim = 2
n_step = 1800
dt = tau * 0.01

# Galaxy 1
n_galaxy_1 = 500
M_1 = Mg
R_1 = Rg


# Galaxy 2
n_galaxy_2 = 200
M_2 = 1.54E10*Msun   #Kg
R_2 = 3.8e17        #Km
x02 = 8*Rg        #Km
y02 = 2.5*Rg        #Km
vx02 = -Rg/(2*tau)        #Km/s
vy02 = 0       #Km/s

#Total
n_particle = n_galaxy_1+n_galaxy_2
