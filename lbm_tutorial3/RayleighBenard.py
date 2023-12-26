'''Showcase for thermalLBM.py Simulation.
Simulating a purely temperature-driven flow: Rayleigh-Benard convection for
given Rayleigh and Prandtl numbers.
'''

# importing necessary classes and packages
from thermalLBM import ThermalLBM
from boundary import Boundary
import numpy as np


# Defining parameters for simulation
Nx = 100
Ny = 50
gy = -9.81 # acceleration of gravity in the y-direction
beta = 0.000101937 # expansion coefficient
Pr = 0.71
Ra = 5000

Ttop=-1
Tbot=1
# Calculating kinetic viscosity and thermal diffusivity from Rayleigh and Prandtl numbers
nu = ((Pr*beta*(Tbot-Ttop)*abs(gy)*(Ny**3))/Ra)**(.5)
alpha = nu/Pr

# Setting up boundary object for flow
flow_boundary = Boundary()
boundary_array = np.zeros((Ny, Nx))
boundary_array[[0, -1], :] = 1

flow_boundary.init_bounceback(boundary_array)#

# Setting up boundary object for temperature
temp_boundary = Boundary()

boundary_array = np.zeros((Ny, Nx))
boundary_array[0, :] = 1#
temp_boundary.init_constant(boundary_array, Ttop)#

boundary_array = np.zeros((Ny, Nx))
boundary_array[-1, :] = 1#
temp_boundary.init_constant(boundary_array, Tbot)#

# Initialize thermal simulation for Rayleigh-Benard convection
RayleighBenard = ThermalLBM(Nx, Ny, nu, alpha, gy*beta, flow_boundary, temp_boundary)

# Run simulation with graphical output
RayleighBenard.run_vis(1)
