import numpy as np

### MY VARIABLES

# Plasma (0) & Beam (b) Electrons
vb = 0.12         # beam velocity
n0 = 1.            # plasma density
T0 = 4.e-4         # plasma thermal velocity
nb = 0.05          # beam density
Tb = 4.e-4         # beam temperature
v0 = -(nb/n0)*vb   # plasma drift velocity

# Numerical parameters
dx    = np.sqrt(T0/n0)      # spatial resolution
dt    = 0.95*dx             # temporal resolution
Lx    = 512*dx              # size of the simulation domain
Tsim  = 64                  # duration of the simulation
nppc0 = 500                 # number of particles_per_cell (plasma)
nppcb = 500                 # number of particles_per_cell (beam)

# Test condition on the instability
# print 'nb/n0 * (v0/vthb)**2 = ',nb/n0*vb**2/Tb,' > ',vb**3/np.sqrt(T0)**3 * np.exp(-vb**2/2./T0)


### SMILEI BLOCKS

### Main() block
Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    grid_length = [Lx],
    cell_length = [dx],
    simulation_time = Tsim,
    timestep = dt,
    number_of_patches = [4],
    EM_boundary_conditions = [["periodic"]],
    solve_poisson = False,
    print_every = int(Tsim/dt/20.),
    random_seed = smilei_mpi_rank
)

### Species()
Species(
    name      = "plasma",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = nppc0,
    mass = 1.,
    charge = -1.,
    number_density = n0,
    mean_velocity = [v0,0.,0.],
    temperature = [T0],
    boundary_conditions = [["periodic"]]
)

Species(
    name      = "beam",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = nppcb,
    mass = 1.,
    charge = -1.,
    number_density = nb,
    mean_velocity = [vb,0.,0.],
    temperature = [Tb],
    boundary_conditions = [["periodic"]]
)

### Diagnostics
globalEvery = int(0.2/dt)

DiagScalar(
    every = globalEvery,
    vars=['Utot','Ukin','Uelm','Uelm_Ex']
)

DiagFields(
    every = globalEvery,
    time_average = 1,
    fields = ["Ex", "By", "Bz", "Rho_plasma", "Rho_beam"],
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    time_average = 1,
    species = ["plasma", "beam"],
    axes = [
        ["px", -5.*np.sqrt(T0), vb+5.*np.sqrt(Tb), 200]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    time_average = 1,
    species = ["plasma", "beam"],
    axes = [
        ["x" , 0., Lx, int(Lx/dx/2.)],
        ["px", -5.*np.sqrt(T0), vb+5.*np.sqrt(Tb), 200]
    ]
)
