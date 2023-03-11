# In this example, we show how to use mesti2s() to compute the scattering
# matrix of a large disordered medium. We want to focus on the performance of
# multithreading (i.e. different # of threads).

# We can tune the # of threads by the environment variable OMP_NUM_THREADS.
# For example, we can set the number of threads to be two by "$export OMP_NUM_THREADS = 2".

# Include mesti2s.jl
include("mesti2s.jl")
# The environmental variable for MUMPS3 (it should be the path to libraries of MUMPS)
ENV["MUMPS_PREFIX"] = "/project/cwhsu_38/hochunli/carc_help/shared_lib/lib"

# Load the essential modules
using MUMPS3 # MUMPS-julia interface
using MPI # MPI module is necessary for MUMPS3
using Random

# System parameters
syst= Syst()
syst.dx = 1
syst.wavelength = 20
syst.epsilon_low = 1
syst.epsilon_high = 1
epsilon_max = 4
epsilon_min = 1
nx = 210; nx_Ey = nx; nx_Ez = nx
ny = 210; ny_Ex = ny; ny_Ez = ny
nz = 20; nz_Ex = nz; nz_Ey = nz
syst.xBC = "periodic"; nx_Ex = nx
syst.yBC = "periodic"; ny_Ey = ny
nz_Ez = nz + 1
syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min
syst.zPML = [PML(10)]

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both s-polarization and p-polarization
input.side = "both"
input.polarization = "both"
# Output to both sides with both s-polarization and p-polarization
output.side = "both"
output.polarization = "both"

opts = Opts()
# Clear variables to reduce peak memory usage
opts.clear_memory = true
opts.clear_BC = true
(S, channels, info)= mesti2s(syst, input, output, opts)
