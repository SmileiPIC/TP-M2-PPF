# -------------------
# MY PYTHON VARIABLES
# -------------------

import math

grid_length = 1.68   # length of the simulation grid

# plasma parameters
velocity=0.1
amplitude=0.0001
density=1.

# numerical parameters
nppc=200                # nb of particles per cell
dx = grid_length/64.    # cell length (spatial resolution)
dt = 0.95*dx            # timestep (0.95% of the CFL)
Tsim = 100.             # duration of the simulation
diagEvery=int(0.1/dt)   # diag output frequency

# perturbation on the electron densities
def ne_(x):
    return density/2. * ( 1.+amplitude*np.cos(2.*math.pi*x/grid_length) )


### ----------------------------------------------------
### HERE STARTS SMILEI VARIABLE DEFINITION: in Blocks()
### ----------------------------------------------------

# Main()
Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    simulation_time = Tsim,
    timestep = dt,
    grid_length = [grid_length],
    cell_length = [dx],
    number_of_patches = [ 4 ],
    print_every = int(Tsim/dt/10.),
    EM_boundary_conditions = [['periodic']],
    #solve_poisson=False
)

# External field (applied at t=0 only)
# ExternalField(
#     field = "Ex",
#     profile = ex_
#)

# Background ion species
Species(
    name = 'ion',
    number_density = density,
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1.0,
    charge = 1.0,
    time_frozen = 2.*Tsim,
    boundary_conditions = [['periodic']]
)

# Forward moving electron species 1
Species(
    name = "electron1",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.0,
    charge = -1.0,
    number_density = ne_,
    boundary_conditions = [['periodic']],
    mean_velocity = [velocity, 0, 0]
)

# Backward moving electron species 2
Species(
    name = "electron2",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = nppc,
    mass = 1.0,
    charge = -1.0,
    number_density = ne_,
    boundary_conditions = [['periodic']],
    mean_velocity = [-velocity, 0, 0]
)

### Diagnostics blocks(): frequency given by diagEvery=int(0.5/dt)
DiagScalar (
    precision = 3,
    every=diagEvery,
    vars=['Utot','Ukin','Uelm','Uelm_Ex','Uelm_Ey','Uelm_Ez','Uelm_Bx_m','Uelm_By_m','Uelm_Bz_m']
)

DiagFields(
    every = diagEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz','Jx','Rho']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diagEvery,
    species = ['electron1', 'electron2'],
    axes = [
        ["x", 0., grid_length, int(grid_length/dx/2.)],
        ["vx", -4*velocity, 4*velocity, 100]
    ]
)
