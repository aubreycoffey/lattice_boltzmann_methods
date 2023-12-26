class Boundary():
    '''Boundary class for thermalLBM and LBM simulations.'''

    def __init__(self):
        self.boundaries = []  # List of all applied boundaries (except bounce-back)
        self.bb_array = False
        self.constant_array = False

    def init_bounceback(self, bb_array):
        '''Initialize position of bounce-back cells.'''

        self.bb_array += (bb_array == 1)

    def init_constant(self, logic_array, constant):
        '''Initialize position and value of constant condition for populations like
        temperature or concentration.
        '''
        try:   # Check for existence of object variable constant_values
            self.constant_values[logic_array == 1] = constant  # Value array
        except AttributeError:
            self.constant_values = 0*logic_array
            self.constant_values[logic_array == 1] = constant

        self.constant_array += (logic_array == 1)  # Logic array indicating position

        if self.constant not in self.boundaries:
            self.boundaries.append(self.constant)

    def constant(self, f, rho, ux, uy):
        '''Enforcing constant value boundary.
        Intended for ADE-populations with constant values on walls.
        Therefore fluid velocity is set to zero.
        NOT FOR CONSTANT DENSITY/PRESSURE CONDITIONS!!!
        '''
        rho[self.constant_array == 1] = self.constant_values[self.constant_array == 1]
        ux[self.constant_array == 1] = 0
        uy[self.constant_array == 1] = 0

    def init_zhvl(self, value):
        '''Initialize Zou-He velocity boundary on the left side.
        Value can be matrix of shape Ny*1 (velocity profile).
        '''
        self.vl = value
        self.boundaries.append(self.zhvl)

    def init_zhpr(self, value):
        '''Initialize Zou-He density/pressure boundary on the right side.'''
        self.pr = value
        self.boundaries.append(self.zhpr)

    def zhvl(self, f, rho, ux, uy):
        '''Enforce Zou-He velocity boundary on the left side.'''
        ux[1:-1, 0] = self.vl
        uy[1:-1, 0] = 0
        rho[1:-1, 0] = (f[0, 1:-1, 0] + f[2, 1:-1, 0] + f[4, 1:-1, 0] + 2 *
                        (f[3, 1:-1, 0] + f[6, 1:-1, 0] + f[7, 1:-1, 0])) / (1-ux[1:-1, 0])

        f[1, 1:-1, 0] = f[3, 1:-1, 0] + 2/3*rho[1:-1, 0]*ux[1:-1, 0]
        f[5, 1:-1, 0] = f[7, 1:-1, 0] + .5 * \
            (f[4, 1:-1, 0] - f[2, 1:-1, 0]) + 1/6*rho[1:-1, 0]*ux[1:-1, 0]
        f[8, 1:-1, 0] = f[6, 1:-1, 0] - .5 * \
            (f[4, 1:-1, 0] - f[2, 1:-1, 0]) + 1/6*rho[1:-1, 0]*ux[1:-1, 0]

    def zhpr(self, f, rho, ux, uy):
        '''Enforce Zou-He density/pressure boundary on the left side.'''
        rho[1:-1, -1] = self.pr
        ux[1:-1, -1] = (f[0, 1:-1, -1] + f[2, 1:-1, -1] + f[4, 1:-1, -1] + 2 *
                        (f[1, 1:-1, -1] + f[5, 1:-1, -1] + f[8, 1:-1, -1]))/rho[1:-1, -1] - 1
        uy[1:-1, -1] = 0

        f[3, 1:-1, -1] = f[1, 1:-1, -1] - 2/3*rho[1:-1, -1]*ux[1:-1, -1]
        f[6, 1:-1, -1] = f[8, 1:-1, -1] - .5 * \
            (f[4, 1:-1, -1] - f[2, 1:-1, -1]) - 1/6*rho[1:-1, -1]*ux[1:-1, -1]
        f[7, 1:-1, -1] = f[5, 1:-1, -1] + .5 * \
            (f[4, 1:-1, -1] - f[2, 1:-1, -1]) - 1/6*rho[1:-1, -1]*ux[1:-1, -1]

    def enforce_boundaries(self, f, rho, ux, uy):
        '''Enforce all initialized boundary conditions in the boundary list.'''
        for boundary_method in self.boundaries:
            boundary_method(f, rho, ux, uy)
