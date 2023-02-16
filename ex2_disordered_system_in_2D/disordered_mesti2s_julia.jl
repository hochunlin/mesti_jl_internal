# In this example, we show how to use mesti2s() to compute the scattering
# matrix of a strongly scattering disordered medium, analyze the scattering
# matrix for the s-polarization. We can get the same result from mesti2s.m.

# Include mesti2s.jl
include("mesti2s.jl")
# The environmental variable for MUMPS3 (it should be the path to libraries of MUMPS)
ENV["MUMPS_PREFIX"] = "/project/cwhsu_38/hochunli/mumps_5_5_0_par_build/src"

# Load the essential modules 
using MUMPS3 # MUMPS-julia interface
using MPI # MPI module is necessary for MUMPS3
using MAT # MAT module for writing and reading mat filess
using PyPlot # PyPlot module for plotting 
 
# Load and reshape epsilon_xx, epsilon_yy, epsilon_zz from 2D code
file = matopen("epsilon_xx.mat")
epsilon_xx = read(file, "epsilon_xx")
epsilon_xx = reshape(epsilon_xx,size(epsilon_xx,1),1,size(epsilon_xx,2))
close(file)
file = matopen("epsilon_yy.mat")
epsilon_yy = read(file, "epsilon_yy")
epsilon_yy = reshape(epsilon_yy,size(epsilon_yy,1),1,size(epsilon_yy,2))
close(file)
file = matopen("epsilon_zz.mat")
epsilon_zz = read(file, "epsilon_zz")
epsilon_zz=reshape(epsilon_zz,size(epsilon_zz,1),1,size(epsilon_zz,2))
close(file)

# System parameters
syst= Syst()
syst.dx = 1
syst.wavelength = 10
syst.epsilon_low = 1
syst.epsilon_high = 1
syst.epsilon_xx = epsilon_xx
syst.epsilon_yy = epsilon_yy
syst.epsilon_zz = epsilon_zz
syst.xBC = "periodic"
syst.yBC = "periodic"
syst.zPML = [PML(20)]

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with s-polarization
input.side = "both"
input.polarization = "s"
# Output to both sides with s-polarization
output.side = "both"
output.polarization = "s"

(S_julia, channels, info)= mesti2s(syst, input, output);

# Load the result from matlab
file = matopen("S_matlab.mat")
S_matlab = read(file, "S_matlab")
close(file)

# Compare them
pcolormesh(abs.(abs.(S_julia)-abs.(S_matlab))) 
xlabel("Input channels", fontsize=16)
ylabel("Output channels", fontsize=16)
colorbar()
savefig("s_matrix_diff_in_disordered_between_matlab_and_julia")