import numpy as np





class BasicLBM:
    '''Basic class with elementary methods for Lattice Boltzmann simulations.'''

    def __init__(self, Nx, Ny):
        '''Initialization of grid size and temporary population array for bounceback boundary handling.
        '''
        self.Nx = Nx
        self.Ny = Ny
        self.tmp = np.zeros((9, Ny, Nx))

    def density(self, f):
        '''Calculation of the density for a given population f.'''
        return np.sum(f, axis=0)

    def velocity(self, f, density):
        '''Calculation of the velocities for a given population f.'''
        ux = ((f[1, :, :] + f[5, :, :] + f[8, :, :])-(f[6, :, :] + f[3, :, :] + f[7, :, :])) /density
        ux[np.isnan(ux)]=0
        uy = ((f[6, :, :] + f[5, :, :] + f[2, :, :])-(f[8, :, :] + f[4, :, :] + f[7, :, :])) /density
        uy[np.isnan(uy)]=0
        return ux, uy

    def update_eq(self, feq, density, ux, uy):
        '''Updating equilibrium distributions given macrosopic density and velocities.'''
        u = ux**2 + uy**2  # Calculating square of velocity magnitude

        w1 = 1/9  # Weights for lattice directions
        w5 = 1/36

        # Using discrete equilibrium distribution function for D2Q9 lattice
        feq[0, :, :] = 4/9 * density * (1 - 3/2*u)

        feq[1, :, :] = w1 * density * (1 + 3*ux + 9/2*ux**2 - 3/2*u)
        feq[2, :, :] = w1 * density * (1 + 3*uy + 9/2*uy**2 - 3/2*u)
        feq[3, :, :] = w1 * density * (1 - 3*ux + 9/2*ux**2 - 3/2*u)
        feq[4, :, :] = w1 * density * (1 - 3*uy + 9/2*uy**2 - 3/2*u)

        #feq[5, :, :] = w5 * density * (1 + 3*(ux+uy) + 9*ux*uy + 3*u)
        feq[5, :, :] = w5 * density * (1 + 3*(ux+uy) + 9/2*(ux+uy) - 3/2*u)
        feq[6, :, :] = w5 * density * (1 + 3*(-ux+uy) + 9/2*(-ux+uy)**2 - 3/2*u)
        feq[7, :, :] = w5 * density * (1 + 3*(-ux-uy) + 9/2*(-ux-uy)**2 - 3/2*u)
        feq[8, :, :] = w5 * density * (1 + 3*(ux-uy) + 9/2*(ux-uy)**2 - 3/2*u)

    def collide(self, f, feq, omega, density, ux, uy, boundary=False,):
        '''Collision using BGK operator.
        Relaxation from populations f to equilibrium feq according to relaxation parameter omega.
        Macrosopic variables are needed for boundary handling.
        Boundary object (cf. boundary.py) can be passed on.
        '''

        if boundary is False:
            '''Collision without any boundaries.'''
            #self.tmp=np.copy(f)
            self.update_eq(feq, density, ux, uy) # Update equilibrium
            #fOut(i,:,:) = fIn(i,:,:)-(omega*(fIn(i,:,:)-fEq(i,:,:)));
            f = f-omega*(f-feq)  # Relaxation of populations
            return f

        else:
            '''Collision conforming to given boundary object.'''
            boundary.enforce_boundaries(f, density, ux, uy)  # Enforce boundary (cf. boundary.py)

            self.update_eq(feq, density, ux, uy) # Update equilibrium
            f = f*(1-omega) + feq*omega  # Relaxation of populations
            self.tmp = np.copy(f)            # Store post collision values for anti bounce-back

            # Replace boundary populations according to bounce-back rule
            if boundary.bb_array is not False:
                i = boundary.bb_array
                f[1, i], f[3, i] = f[3, i], f[1, i]
                f[2, i], f[4, i] = f[4, i], f[2, i]
                f[5, i], f[7, i] = f[7, i], f[5, i]
                f[6, i], f[8, i] = f[8, i], f[6, i]

            # Replace constant boundary populations according to anti bounce-back
            if boundary.constant_array is not False:
                i = boundary.constant_array
                f[1, i] = -self.tmp[3, i] + 2/9*boundary.constant_values[i]
                f[2, i] = -self.tmp[4, i] + 2/9*boundary.constant_values[i]
                f[3, i] = -self.tmp[1, i] + 2/9*boundary.constant_values[i]
                f[4, i] =-self.tmp[2, i] + 2/9*boundary.constant_values[i]
                f[5, i] =-self.tmp[7, i] + 2/36*boundary.constant_values[i]
                f[6, i] = -self.tmp[8, i] + 2/36*boundary.constant_values[i]
                f[7, i] =-self.tmp[5, i] + 2/36*boundary.constant_values[i]
                f[8, i] = -self.tmp[6, i] + 2/36*boundary.constant_values[i]
            return f

    def stream(self, f):
        '''Streaming of populations f on a D2Q9 lattice.'''
        f[1, :, :] = np.roll(f[1, :, :], 1, axis=1)
        f[2, :, :] = np.roll(f[2, :, :], -1, axis=0)
        f[3, :, :] = np.roll(f[3, :, :], -1, axis=1)
        f[4, :, :] = np.roll(f[4, :, :], 1, axis=0)

        f[5, :, :] = np.roll(np.roll(f[5, :, :], 1, axis=1), -1, axis=0)
        f[6, :, :] = np.roll(np.roll(f[6, :, :], -1, axis=1), -1, axis=0)
        f[7, :, :] = np.roll(np.roll(f[7, :, :], -1, axis=1), 1, axis=0)
        f[8, :, :] = np.roll(np.roll(f[8, :, :], 1, axis=1), 1, axis=0)
