# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import math

l0 = 2.0*math.pi	# laser wavelength
t0 = l0				# optical cicle
Lsim = 60*l0		# length of the simulation
Tsim = 100.*t0		# duration of the simulation
resx = 20 		    # nb of cells in on laser wavelength
rest = resx/0.95	# time of timestep in one optical cycle (0.95 * CFL)

Main(
    geometry = "1Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx],
    grid_length  = [Lsim],
    
    number_of_patches = [ 8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [ ['silver-muller'] ],
     
    random_seed = smilei_mpi_rank,
    clrw = 1
)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = 100,
    mass = 1836.0,
    charge = 1.0,
    number_density = trapezoidal(1./16.,xvacuum=5*l0,xplateau=50*l0),
    temperature = [0.],
    boundary_conditions = [
        ["reflective", "reflective"],
    ],
)
Species(
    name = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'cold',

    particles_per_cell = 100,
    mass = 1.0,
    charge = -1.0,
    number_density = trapezoidal(1./16.,xvacuum=5*l0,xplateau=50*l0),
    temperature = [0.05],
    boundary_conditions = [
        ["reflective", "reflective"],
    ],
)

LaserPlanar1D(
	box_side = 'xmin',
	a0 = ,
    omega = 1.,
    ellipticity = 0.,
    time_envelope =  tgaussian(start=0, fwhm=, order=2),
)


every = int(rest/2.)

DiagFields(
    every = every,
    fields = []
)

DiagScalar(every=every)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = every,
    species = ["eon"],
    axes = [
        ["x",  0.,   Lsim , 200],
        ["px", -10., 20., 200]
    ]
)
