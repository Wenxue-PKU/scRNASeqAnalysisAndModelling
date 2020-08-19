import numpy as np
import numba as nb
from pde import PDEBase, CartesianGrid, ScalarField, FieldCollection, PlotTracker, FileStorage

# Define the number of ligands and receptors
num_ligands = 2
num_receptors = 3
num_complexes = 4

# Define the PDE class for a single ligand and multiple receptors 
class LigandReceptorsPDE(PDEBase):
    """ PDE class consisting of a model that considers n ligand and m receptors species.
        We assume that each ligand can bind to multiple receptors and each receptor
        can bind to multiple ligands. Assuming single ligand-receptor complexes, this gives
        n * m possible complexes and n + m + n * m total species to consider. We will however
        refine the model as we go along.
        
        @param diffs, n x 1 vector of diffusivities of the ligand molecules
        @param prod_rates, n x 1 Field collection for production rates (based on expression) for ligand molecules
        @param bind_rates, n x m matrix of binding affinities between ligand-receptor pairs
        @param diss_rates, n x m matrix of disassociation rates of ligand-receptor complexes
        @param degra_rates, (n + m + l) x 1 vector of degradation rates of ligands, receptors, then complexes,
        @param bc, specifies the boundary condition. For now, we go with Neumann BCs
    """
    def __init__(self, diffs, prod_rates, bind_rates, diss_rates, degrad_rates, bc):
        """ Initialise the diffusivity of the ligand, the binding affinities of the receptors, the degradation
            rates of the molecules, and the natural (no-flux boundary conditions. """
        self.diffs = diffs # Ligand diffusivities
        self.prod_rates = prod_rates # Production rate sof ligands
        self.bind_rates = bind_rates # Binding affinities 
        self.diss_rates = diss_rates # Disassociation rates of complexes
        self.degrad_rates = degrad_rates # Degradation rates
        self.num_ligands = num_ligands # Number of ligand species
        self.num_receptors = num_receptors # Number of receptor species
        self.num_complexes = num_complexes # Number of bound complexes
        self.bc = bc # Boundary condition

    def evolution_rate(self, state, t=0):
        """ Function to define the PDEs for the ligand-receptor model. We assume only the ligand molecule
            can diffuse in space, which is then consumed by the receptors. """
        diffs = self.diffs
        bind_rates = self.bind_rates
        diss_rates = self.diss_rates
        prod_rates = self.prod_rates
        degrad_rates = self.degrad_rates

        # Initialise rhs
        rhs = state.copy()
        # First figure out where the non-zero binding (and dissociation rates) are.
        # We're kind of assuming that a complex will have  non-zero binding and
        # dissociation (may need to change this later).
        nonzero_indices = np.nonzero(bind_rates)
        bound_ligand_indices = nonzero_indices[0]
        bound_receptor_indices = nonzero_indices[1]
        nonzero_bind_rates = bind_rates[bound_ligand_indices, bound_receptor_indices]
        nonzero_diss_rates = diss_rates[bound_ligand_indices, bound_receptor_indices]

        for i in range(num_ligands + num_receptors + num_complexes):
            # Define the state 
            u = state[i]

            # Construct PDEs
            if (i < num_ligands): # For ligands, we add diffusion, a production term, and degradation
                u_t = diffs[i] * u.laplace(self.bc) + prod_rates[i] - degrad_rates[i]*u
                
                ligand_index = i # Ligands are first in the list, so match 1-1 with their index

                # Account for any binding to receptors that may occur
                if (ligand_index in bound_ligand_indices):
                    bound_receptors = np.where(bound_ligand_indices == ligand_index)[0] # Get the receptors bound to this ligand
                    for j in bound_receptors:
                        receptor_index = int(bound_receptor_indices[int(j)])
                        # Account for binding with all possible receptor partners
                        receptor_u = state[receptor_index + num_ligands] # Define the receptor state
                        complex_u = state[num_ligands + num_receptors + ligand_index * (num_receptors - 2) + receptor_index]  # Define the complex state based on natural ordering
                        u_t += diss_rates[ligand_index, receptor_index] * complex_u - bind_rates[ligand_index, receptor_index] * u * receptor_u
                    
                rhs[i] = u_t # Add equation to list

            elif ((i >= num_ligands)&(i < num_ligands + num_receptors)): # For receptors we consider binding, dissociation and degradation
                u_t = -1.0 * degrad_rates[i] * u

                receptor_index = i - num_ligands # Shift the index across to determine the receptor index (wrt the bindin gmatrix)

                if (receptor_index in bound_receptor_indices):
                    bound_ligands = np.where(bound_receptor_indices == receptor_index)[0] # Get the receptors bound to this ligand

                    for j in bound_ligands:
                        ligand_index = int(bound_ligand_indices[int(j)])
                        # Account for binding and dissociation
                        ligand_u = state[ligand_index]
                        complex_u = state[num_ligands + num_receptors + ligand_index * (num_receptors - 2) + receptor_index]
                        u_t += diss_rates[ligand_index, receptor_index] * complex_u - bind_rates[ligand_index, receptor_index] * ligand_u * u

                rhs[i] = u_t # Add equation to list

            else: # Complexes are formed by ligands binding to receptors and dissociate/degrade

                # Shift index over to get the complex index
                complex_index = i - (num_ligands + num_receptors)

                # Work out where the binding rates are equal to the current binding rate considered
                bind_rate = nonzero_bind_rates[complex_index]
                diss_rate = nonzero_diss_rates[complex_index]
                degrad_rate = degrad_rates[i]

                # The actual index pair is the pair such that p*m + q = i
                ligand_index = int(bound_ligand_indices[complex_index]) # # Ligand index
                receptor_index = int(bound_receptor_indices[complex_index]) + num_ligands# Receptor index
                u_ligand = state[ligand_index] # Associated ligand
                u_receptor = state[receptor_index] # Associated receptor
                u_t = bind_rate * u_ligand * u_receptor - diss_rate * u - degrad_rate * u # Binding, dissociation, then degradation
            
                rhs[i] = u_t

        return rhs


# Define the initial grid, which corresponds to 
grid = CartesianGrid(([0, 1], [0, 1]), (40, 40))

# Define the parameters
diffusivities = np.array([0.001, 0.001]) 
prod_scalars = np.array([0.0001, 0.0001])
# production_rates = [prod_scalars[0]*ScalarField.from_expression(grid, '0.5 * (tanh(50.0*(y - 40*0.2)) + tanh(50.0*(0.3*40 - y)))'), \
#     prod_scalars[1]*ScalarField.from_expression(grid, '0.5 * (tanh(50*(y - 0.7*40)) + tanh(50*(0.8*40 - y)))')]
production_rates = [prod_scalars[0]*ScalarField.from_expression(grid, '0.5 * (tanh(50.0*(y - 40*0.2)) + tanh(50.0*(0.3*40 - y)))'), \
    prod_scalars[1]*ScalarField.from_expression(grid, '0.5 * (tanh(50*(y - 0.2*40)) + tanh(50*(0.3*40 - y)))')]
binding_rates = np.array([[0.1, 0.0, 0.1], [0.1, 0.1, 0.0]])
dissociation_rates = np.array([[0.001, 0.0, 0.001], [0.001, 0.001, 0.0]])
degradation_rates = np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

# Define the initial data
initial_data = []
# Get the locations of the non-zero binding rates
nonzero_bind_indices = np.nonzero(binding_rates)
bound_ligand_indices = nonzero_bind_indices[0]
bound_receptor_indices = nonzero_bind_indices[1]

# Let's initialise the data
for i in range(num_ligands + num_receptors + num_complexes):

    if (i < num_ligands):

        init_label = 'Ligand ' + str(i + 1)
        init_data = ScalarField(grid, 0, label=init_label)
        init_data.data[:, (8):(13)] = 1.0

        initial_data.append(init_data)

    elif ( (i >= num_ligands)&(i < num_ligands + num_receptors) ):

        init_label = 'Receptor ' + str(i - num_ligands + 1)
        init_data = ScalarField.random_uniform(grid, 0, 1, label=init_label)
        initial_data.append(init_data)

    else:

        ligand_index = bound_ligand_indices[i - (num_ligands + num_receptors)] # # Ligand index
        receptor_index = bound_receptor_indices[i - (num_ligands + num_receptors)] # Receptor index
        init_label='Complex ' + str(ligand_index + 1) + ', ' + str(receptor_index + 1)
        init_data = ScalarField(grid, 0, label=init_label)
        initial_data.append(init_data)

initial_data = FieldCollection(initial_data)

# This code is used to store the initial data as a h5 file
# initialDataFileName = './Results/LigandsSameRegion/ligandreceptormodel_n_' + str(num_ligands) + '_m_' + str(num_receptors) \
#             + '_p_' + str(num_complexes) + '_sameregion_initialcondition.h5'
# initial_data.to_file(initialDataFileName)
# Load this initial data
# initial_data = FieldCollection.from_file(initialDataFileName)

# # # # Define the scaling factors for the various parameters we alter
# scaling_factors = [(10.0)**i for i in range(-2, 6)] # We will multiply ligand 2's factors by 10^(n), where n is specified by this array.

t_end = 10 # End timestep
t_track = 0.25 # How often we track simulations
dt = 0.025

# Define the PDE model
eq = LigandReceptorsPDE(diffs=diffusivities, \
                        prod_rates=production_rates, \
                        bind_rates=binding_rates, \
                        diss_rates=dissociation_rates, \
                        degrad_rates=degradation_rates, bc='natural')

fileName = './Results/LigandsSameRegion/ligandreceptormodel_n_' + str(num_ligands) + '_m_' + str(num_receptors) \
            + '_p_' + str(num_complexes) + '.h5'
storage_write = FileStorage(fileName)
plot_tracker = PlotTracker(interval=t_track, scale='unity') # To plot the solutions at certain time intervals
trackers = ['progress', plot_tracker, storage_write.tracker(t_track)] # To save the solutions
sol = eq.solve(initial_data, t_range=t_end, dt = dt, tracker=trackers, method='scipy') # Run the solution
