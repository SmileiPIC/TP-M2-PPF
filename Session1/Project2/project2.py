# -------------------
# MY PYTHON VARIABLES
# -------------------

from math import pi, sin, sqrt

k  = 2.	            # seeded wave number

# plasma parameters
n0 = 1.					# electron density
v0 = 0.9				# electron flow velocity
g0 = 1./sqrt(1.-v0**2)	# Lorentz factor corresponding to v0

# numerical parameters
Lx    = 2.*2.*pi/k		# simulation length (4 seeded wavelength)
t_sim = 100. 			# simulation times
dx    = 2.*pi/k /128.	# spatial resolution in the x-direction = 500 points in one wavelength
dt    = 0.95*dx         # time-step
nppc  = 200				# nb of particles per cell
diagEvery = max(1,int(0.1/dt)) # frequency diag output


# define the function used to seed the Weibel unstable mode
dBy = 0.0001
def By(x):
	return dBy*sin(k*x)

# --------------------------------------
# SMILEI's VARIABLES (defined in blocks)
# --------------------------------------

# MAIN BLOCK
Main(
    geometry = "1Dcartesian",
    interpolation_order = 4,
    timestep = dt,
    simulation_time = t_sim,
    cell_length = [dx],
    grid_length  = [Lx],
    number_of_patches = [ 4 ],
    EM_boundary_conditions = [ ['periodic'] ] ,
    random_seed = smilei_mpi_rank
)

# APPLY EXTERNAL FIELD (to seeded unstable mode)
ExternalField(
	field='By',
    profile = By
)

# DEFINE ALL SPECIES (ADDING BLOCKS)
Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1836.,
    charge = 1.0,
    number_density = 1.,
    boundary_conditions = [['periodic']],
    time_frozen = 2.*t_sim
)

Species(
    name = 'electron1',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = 0.5,
    mean_velocity = [0.,0.,v0],
    boundary_conditions = [['periodic']]
)

Species(
    name = 'electron2',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = 0.5,
    mean_velocity = [0.,0.,-v0],
    boundary_conditions = [['periodic']]
)

# ---------------------
# DIAGNOSTIC PARAMETERS
# ---------------------

DiagScalar(
    every=diagEvery,
	vars=['Utot','Ukin','Uelm','Uelm_Ex','Uelm_Ey','Uelm_Ez','Uelm_Bx_m','Uelm_By_m','Uelm_Bz_m']
)

DiagFields(
    every = diagEvery,
    fields = ['Ex','Ey','Ez','By_m','Rho','Rho_electron1','Rho_electron2','Jz','Jz_electron1','Jz_electron2']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = diagEvery,
    species = ['electron1', 'electron2'],
    axes = [
        ["x", 0., Lx, int(Lx/dx/2.)],
        ["pz", -4*g0*v0, 4*g0*v0, 100]
    ]
)
