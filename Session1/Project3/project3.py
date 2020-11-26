import numpy as np

### HERE I FIRST DEFINE MY OWN VARIABLES ###############################################################################

### Definition of the Lamgmuir wave parameters -------------------------------------------------------------------------
klde = 0.4          # k * lambda Debye parameter
eps  = 0.05             # delta n / n0 (amplitude of the density perturbation associated to the Langmuir wave)

### Definition of the plasma parameters
n0   = 1.               # plasma density (we normalize all times to the electron plasma frequency)
T0   = 1.e-4            # electron temperature in units of me c^2

lde  = np.sqrt(T0/n0)   # Debye length in normalized units (electron skin-depth)
vth  = np.sqrt(T0)      # electron thermal velocity in normalized units (speed of light)

k    = klde/lde                # wavenumber the Langmuir wave
w    = np.sqrt(1.+3.*klde**2)  # frequency of the Langmuir wave

# Profile associated to the Langmuir wave
def delta_(x):
    return eps*np.cos(k*x)

def ne_(x):             # electron density profile: carries the perturbation associated to the Langmuir wave
    return n0 *( 1 + delta_(x) )

def vx_(x):
    return -w/k * n0 * delta_(x)


### Simulation box -----------------------------------------------------------------------------------------------------
N    = 1.                     # number of Landmuir wave periods in the simulation box
Lx   = N * 2.*np.pi/k         # length of the simulation box

# # resolution & number of points are computed such that dx ~ lde & Nx is a multiple of 16
# dx   = lde
# Nx   = Lx/dx
# Nx   = int( (Nx+8.)/16.)*16   # final number of cells in the simulation (a multiple of 16 is ensured)
# Nx   = max(Nx,64)
# dx   = Lx/Nx                  # corresponding final resolution
Nx = 64
dx = Lx/Nx

# simulation time & resolution
Tsim = 20.                   # simulation time in normalized units (inverse of electron plasma frequency)
dt   = 0.95*dx               # timestep for the simulation (95% of the CFL)



### NOW I DEFINE SMILEI's VARIABLES (defined in Blocks) ################################################################

### Main() block
Main(
    geometry = "1Dcartesian",
    interpolation_order = 4,
    grid_length = [Lx],
    cell_length = [dx],
    simulation_time = Tsim,
    timestep = dt,
    number_of_patches = [4],
    EM_boundary_conditions = [["periodic"]],
    print_every = int(Tsim/dt/20.),
    random_seed = smilei_mpi_rank
)

### Species('ion')
Species(
    name = "ion",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 1,
    mass = 1.,
    charge = 1.,
    number_density = n0,
    boundary_conditions = [["periodic"]],
    time_frozen = 2.*Tsim,
)

### Species('electron')
Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = 4096,
    mass = 1.,
    charge = -1.,
    number_density = ne_,
    mean_velocity = [vx_,0.,0.],
    temperature = [T0],
    boundary_conditions = [["periodic"]],
    time_frozen = 0.
)

### Diagnostics

diagEvery = int(0.1/dt)

DiagScalar(
    every=diagEvery,
    vars=['Utot','Ukin','Uelm','Uelm_Ex','Uelm_Ey','Uelm_Ez','Uelm_Bx_m','Uelm_By_m','Uelm_Bz_m']
)

DiagFields(
    every = diagEvery,
    fields = ["Ex", "Rho_electron", "Jx_electron"]
)

DiagParticleBinning(
    every = diagEvery,
    deposited_quantity = "weight",
    species = ["electron"],
    axes = [
        ["vx", -6*vth, 6*vth, 200]
    ]
)

DiagParticleBinning(
    every = diagEvery,
    deposited_quantity = "weight",
    species = ["electron"],
    axes = [
        ["x", 0.,      Lx,    int(Nx/2)],
        ["vx", -6*vth, 6*vth, 200]
    ]
)
