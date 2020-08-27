import numpy as np
import numba as nb
from pde import PDEBase, CartesianGrid, ScalarField, FieldCollection, PlotTracker, FileStorage

### This is code for a reaction-diffusion model to investigate cell-cell communication in the form of
### ligand-receptor binding, where the ligand diffuses and is secreted at one end.
### In particular, we are interested in how well a package like SoptSC differentiates between
### the ligand/receptor expressions in differentiating cell-cell communication.
### Ultimately, we want to extend SoptSC to differentiate between communiciation via cell-cell contact
### and communication via diffusion of a ligand.

# Define the number of ligands and receptors
num_ligands = 1
num_receptors = 1
num_complexes = 1

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
        # We're kind of assuming that a complex will have non-zero binding and
        # dissociation (may need to change this later).
        
        # If there's only one complex, so one ligand, one receptor, our job is easy
        if (np.size(bind_rates) == 1):

            # Define the state variables
            u_ligand = state[0]
            u_receptor = state[1]
            u_complex = state[2]

            # First equation models diffusion of the ligand, plus production of the ligand, as well as its binding, increase due to dissociation, and natural degradation
            rhs[0] = diffs[0] * u_ligand.laplace(self.bc) + prod_rates[0] - bind_rates[0] * u_ligand * u_receptor + diss_rates[0] * u_complex - degrad_rates[0] * u_ligand

            # Second equation models binding of the receptor, increase due to dissociation and natural degradation
            rhs[1] = -bind_rates[0] * u_ligand * u_receptor + diss_rates[0] * u_complex - degrad_rates[1] * u_receptor

            # Third equation models increase in complex due to binding, decrease due to dissociation and degradation
            rhs[2] = bind_rates[0] * u_ligand * u_receptor - diss_rates[0] * u_complex - degrad_rates[2] * u_complex

        else: # Else we need to do a bit more work
        
        # First figure out where the non-zero binding (and dissociation rates) are.
        # We're kind of assuming that a complex will have non-zero binding and
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


# Define the initial grid, which here will just be a 1D line.
grid = CartesianGrid([[0, 1]], 100)

# Some characteristic length scales
tissue_length = 1000.0 # Length of tissue
cell_length = 10.0 # Length of scale
wavelength = cell_length/tissue_length
# Define the parameters
diffusivities = np.array([0.001]) 
prod_scalars = np.array([0.01])
wavelength = cell_length/tissue_length # Needed for the initial condition
production_rates = [prod_scalars[0]*ScalarField.from_expression(grid, '0.5 * (1 + tanh(50.0*(' + str(2.0 * wavelength) + ' - x)))')]
binding_rates = np.array([0.1])
dissociation_rates = np.array([0.001])
degradation_rates = np.array([0.001, 0.001, 0.001])

# Define the initial data
initial_data = []

# Let's initialise the data
if (num_ligands + num_receptors + num_complexes) > 3: # If we have more than just one ligand, one receptor, and one bound complex, we can run this

    # Get the locations of the non-zero binding rates
    nonzero_bind_indices = np.nonzero(binding_rates)
    bound_ligand_indices = nonzero_bind_indices[0]
    bound_receptor_indices = nonzero_bind_indices[1]

    for i in range(num_ligands + num_receptors + num_complexes):

        if (i < num_ligands):

            init_label = 'Ligand ' + str(i + 1)
            init_data = ScalarField(grid, 0, label=init_label)
            init_data.data[0:9] = 1.0

            initial_data.append(init_data)

        elif ( (i >= num_ligands)&(i < num_ligands + num_receptors) ):

            init_label = 'Receptor ' + str(i - num_ligands + 1)
            init_data = ScalarField.from_expression(grid, '1 - sin(2 * pi * x /' + str(2 * wavelength) + ')', label=init_label)
            initial_data.append(init_data)

        else:

            ligand_index = bound_ligand_indices[i - (num_ligands + num_receptors)] # # Ligand index
            receptor_index = bound_receptor_indices[i - (num_ligands + num_receptors)] # Receptor index
            init_label='Complex ' + str(ligand_index + 1) + ', ' + str(receptor_index + 1)
            init_data = ScalarField(grid, 0, label=init_label)
            initial_data.append(init_data)
else:

    # Initialise the ligand
    init_label = 'Ligand 1'
    init_data = ScalarField(grid, 0, label=init_label)
    init_data.data[0:9] = 1.0

    initial_data.append(init_data)

    # Initialise the receptor
    init_label = 'Receptor 1'
    print('1.0 - sin(2.0 * pi * x /' + str(2.0 * wavelength) + ')')
    init_data = ScalarField.from_expression(grid, '0.5 * (1.0 - sin(2.0 * pi * x /' + str(2.0 * wavelength) + '))', label=init_label)
    print(init_data.data)
    initial_data.append(init_data)

    # Initialise the complex
    init_label = 'Complex 1, 1'
    init_data = ScalarField(grid, 0, label=init_label)
    initial_data.append(init_data)

initial_data = FieldCollection(initial_data)

init_data.grid
# This code is used to store the initial data as a h5 file
# initialDataFileName = './Results/LigandsSameRegion/ligandreceptormodel_n_' + str(num_ligands) + '_m_' + str(num_receptors) \
#             + '_p_' + str(num_complexes) + '_sameregion_initialcondition.h5'
# initial_data.to_file(initialDataFileName)
# Load this initial data
# initial_data = FieldCollection.from_file(initialDataFileName)

# # # # Define the scaling factors for the various parameters we alter
# scaling_factors = [(10.0)**i for i in range(-2, 6)] # We will multiply ligand 2's factors by 10^(n), where n is specified by this array.

t_end = 5 # End timestep
t_track = 0.1 # How often we track simulations
dt = 0.025

# Define the PDE model
eq = LigandReceptorsPDE(diffs=diffusivities, \
                        prod_rates=production_rates, \
                        bind_rates=binding_rates, \
                        diss_rates=dissociation_rates, \
                        degrad_rates=degradation_rates, bc='natural')

# fileName = './Results/LigandsSameRegion/ligandreceptormodel_n_' + str(num_ligands) + '_m_' + str(num_receptors) \
#             + '_p_' + str(num_complexes) + '.h5'
# storage_write = FileStorage(fileName)
plot_tracker = PlotTracker(interval=t_track) # To plot the solutions at certain time intervals
# trackers = ['progress', plot_tracker, storage_write.tracker(t_track)] # To save the solutions
sol = eq.solve(initial_data, t_range=t_end, dt = dt, tracker=['progress', plot_tracker], method='scipy') # Run the solution