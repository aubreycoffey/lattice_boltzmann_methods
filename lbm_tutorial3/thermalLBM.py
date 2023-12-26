import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from lbm import BasicLBM


class ThermalLBM(BasicLBM):
    '''Class for thermal Lattice Boltzmann simulations.'''

    def __init__(self, Nx, Ny, nu, alpha, thermal_coeff, flow_boundary, temp_boundary):
        BasicLBM.__init__(self, Nx, Ny)

        # Calculate relaxation parameters
        self.omega_nu = 1 / (3*nu + .5)
        self.omega_alpha = 1 / (3*alpha + .5)

        # Output resulting simulation parameters
        print('Simulation initialized for:')
        print('Nx = %i \nNy = %i' % (Nx, Ny))
        print('nu = %f \nalpha = %f\nomega_nu = %f \nomega_alpha = %f' %
              (nu, alpha, self.omega_nu, self.omega_alpha))
        print('thermal_coeff = %f' % (thermal_coeff))

        self.thermal_coeff = thermal_coeff
        self.flow_boundary = flow_boundary
        self.temp_boundary = temp_boundary

        # Set up macroscopic and microscopic variables
        self.feq = np.empty((9, Ny, Nx))
        self.rho = np.ones((Ny, Nx))
        self.geq = np.empty((9, Ny, Nx))
        self.theta = np.zeros((Ny, Nx))
        self.ux = np.zeros((Ny, Nx))
        self.uy = np.zeros((Ny, Nx))

        # Slight perturbation of density
        self.rho += np.random.random((self.Ny, self.Nx))*.001  

        # Set initial populations to equilibrium
        self.update_eq(self.feq, self.rho, self.ux, self.uy)
        self.f = self.feq.copy()

        self.update_eq(self.geq, self.theta, self.ux, self.uy)
        self.g = self.geq.copy()

    def buoyancy(self):
        '''Including buoyant body force using Boussinesq approximation.'''
        f, theta, rho = self.f, self.theta, self.rho
        coeff = self.thermal_coeff

        f[2, :, :] += -1/3 * coeff * theta * rho
        f[4, :, :] += 1/3 * coeff * theta * rho
        f[[5, 6], :, :] += -1/12 * coeff * theta * rho
        f[[7, 8], :, :] += 1/12 * coeff * theta * rho
        
        #self.f=f
        
    def run(self, n):
        '''Running the simulation for n steps.'''
        #t1 = time.clock()  # Keep track of time for performance evaluation
        t1 = time.perf_counter()
        for _ in range(n):

            # Compute macroscopic quantities
            self.rho, self.theta = self.density(self.f), self.density(self.g)# np.sum(self.f,axis=0), np.sum(self.g, axis=0)
            self.ux, self.uy = self.velocity(self.f,self.rho)
            
            # Collision step
            self.g = self.collide(self.g,self.geq,self.omega_alpha,self.theta,self.ux,self.uy,boundary=self.temp_boundary)
            self.f = self.collide(self.f,self.feq,self.omega_nu,self.rho,self.ux,self.uy,boundary=self.flow_boundary)

            # Buoyany term addition
            self.buoyancy()
            
            # Streaming step
            #self.g=self.stream(self.g)
            #self.f=self.stream(self.f)
            #self.g = self.g +self.omega_alpha*(self.g - self.geq)
            self.stream(self.f)
            self.stream(self.g)



        #t2 = time.clock()
        t2 = time.perf_counter()
        # Printing time and resulting MLUPS (million lattice updates per second)
        print('Calculation time for %i loops: %fs \n %f MLUPS' %
              (n, t2-t1, n*self.Ny*self.Nx/(t2-t1)/1000000))

    def run_vis(self, n):
        '''Visualizing velocity and temperature results every n steps.'''

        fig = plt.figure()
        subplot1 = fig.add_subplot(1, 2, 1)
        subplot2 = fig.add_subplot(1, 2, 2)

        X, Y = np.meshgrid(np.linspace(0, self.Nx, self.Nx), np.linspace(self.Ny, 0, self.Ny))

        def update(i):
            self.run(n)
            subplot1.cla()
            subplot2.cla()
            subplot1.set_title('Velocity')
            subplot2.set_title('Temperature')

            u = np.sqrt(self.ux**2+self.uy**2)

            subplot1.imshow(u, cmap=plt.get_cmap('viridis'), extent=[0, self.Nx, 0, self.Ny])
            # subplot1.streamplot(X, Y, self.ux, self.uy, color='k')
            subplot2.imshow(self.theta, cmap=plt.get_cmap('jet'), vmin=-1,
                            vmax=1, extent=[0, self.Nx, 0, self.Ny], interpolation='bilinear')
            # subplot2.contour(X, Y, self.theta, colors='k')

            subplot1.axis([0, self.Nx, 0, self.Ny])
            subplot2.axis([0, self.Nx, 0, self.Ny])

        a = animation.FuncAnimation(fig, update, frames=1, repeat=True, interval=10)
        plt.show()
