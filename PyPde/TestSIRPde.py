import numpy as np
import numba as nb
from pde import UnitGrid, SphericalGrid, ScalarField, FieldCollection, PDEBase, PlotTracker

class SIRPDE(PDEBase):
    """ SIR model with spatial diffusivity for mobility """ 

    # Apparently we need to set this to be false for numba
    check_implementation = False

    def __init__(self, params=[0.3, 0.9, 0.1], bc='natural'):
        self.params = params # Parameters are infectivity rate, recovery rate, and diffusivity
        self.bc = bc # Boundary condition

    def get_initial_state(self, s, i):
        """ Generate the initial state"""
        norm = (s + i).data.max() # Normalise the data via the maximum
        if (norm > 1):
            s /= norm
            i /= norm
        s.label = 'susceptible'
        i.label = 'infected'

        # we assume a closed population so r = 1 - s - i
        r = ScalarField(s.grid, data = 1 - s - i, label='recovered')
        return FieldCollection([s, i, r])

    def evolution_rate(self, state, t = 0):
        """ Define the PDE equations """
        s, i, r = state # Define the dependent variables
        rhs = state.copy() 
        beta, gamma, diff = self.params # Define the parameters
        rhs[0] = diff * s.laplace(self.bc) - beta * s * i
        rhs[1] = diff * i.laplace(self.bc) + beta * s * i  - gamma * i
        rhs[2] = diff * r.laplace(self.bc) + gamma * i
        return rhs

    def _make_pde_rhs_numba(self, state):
        """ Overwritten class to increase compiiling speed via numba """
        # Define state and parameters
        # s_sub, i_sub, r_sub = state
        beta, gamma, diff = self.params

        # Define laplacian on grid
        state_laplace = state.grid.get_operator('laplace', bc=self.bc)

        @nb.jit 
        def pde_rhs(state_data, t):
            """ Compiled helper function to evaluate the RHS """

            # Get the states
            s = state_data[0]
            i = state_data[1]
            r = state_data[2]

            # Define the PDE RHS
            rate = np.empty_like(state_data)
            rate[0] = diff * state_laplace(s) - beta * s * i
            rate[1] = diff * state_laplace(i) + beta * s * i - gamma * i
            rate[2] = diff * state_laplace(r) + gamma * i
            return rate

        return pde_rhs

# Define the PDE with the parameters
sir_eq = SIRPDE(params=[2.0, 0.3, 0.01])

# Initialise the state
grid = UnitGrid([64, 64])
s = ScalarField(grid, 1)
i = ScalarField(grid, 0)
i.data[28:36, 28:36] = 1
initial_state = sir_eq.get_initial_state(s, i)

# Now solve the PDE
tracker = PlotTracker(interval=1, scale='unity')
sol = sir_eq.solve(initial_state, t_range=20, dt = 1e-2, tracker=['progress', tracker], method='scipy')

