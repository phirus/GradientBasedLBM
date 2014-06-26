# dimensionless numbers
global const REYNOLDS_MAX_INI = 100
global const MORTON = 0.01
global const EOTVOS = 10

# several parameters
global const RHO_L = 1000
global const DIAMETER = 0.01
global const RESOLUTION = 40
global const MACH_MAX = 0.1

# timetracking
global const REFINE_FACTOR = 1.1
global const MAX_ITER = 2e6

# termination criterion
global const RESIDUAL_RE = 1e-6

# output settings
global const TECH_PLOT = 1000
global const RESTART = 10000

# Parameters for relaxation matrix
# bulk viscosity -> s_2
global const MU_RATIO = 2
# s_3 + s_5
global const S_3 = 1.2
global const S_5 = 1.2


println("Parameter wurden geladen")
