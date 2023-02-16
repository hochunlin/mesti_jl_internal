# In this example, we show how to use mesti2s() to compute the scattering
# matrix in a homogeneous space and compare the result with the analytic
# prediction.

# Include mesti2s.jl
include("mesti2s.jl")
# The environmental variable for MUMPS3 (it should be the path to libraries of MUMPS)
ENV["MUMPS_PREFIX"] = "/project/cwhsu_38/hochunli/mumps_5_5_0_par_build/src"

# Load the essential modules 
using MUMPS3 # MUMPS-julia interface
using MPI # MPI module is necessary for MUMPS3
using PyPlot # PyPlot module for plotting 

syst= Syst()
syst.dx = 1
syst.wavelength = 5
epsilon_bg = 1
nx = 10; nx_Ey = nx; nx_Ez = nx
ny = 10; ny_Ex = ny; ny_Ez = ny
nz = 11; nz_Ex = nz; nz_Ey = nz
syst.xBC = "periodic"
syst.yBC = "periodic"
nx_Ex = nx # Note that when periodic BC is along x-direction, nx_Ex = nx_Ey = nx_Ez
ny_Ey = ny # Note that when periodic BC is along y-direction, ny_Ey = ny_Ex = ny_Ez
nz_Ez = nz + 1 # Note that when PEC BC is along z-direction, nz_Ez = nz_Ex + 1 = nz_Ey + 1

syst.epsilon_xx = 1*epsilon_bg*ones(nx_Ex,ny_Ex,nz_Ex)
syst.epsilon_yy = 1*epsilon_bg*ones(nx_Ey,ny_Ey,nz_Ey)
syst.epsilon_zz = 1*epsilon_bg*ones(nx_Ez,ny_Ez,nz_Ez)

# Optimized PML parameters for such thickness and resolution
PML_z = PML(60)
PML_z.kappa_max = 7.6489
PML_z.sigma_max_over_omega = 5.3158
PML_z.power_sigma = 7.5019
PML_z.alpha_max_over_omega = 29.3340
PML_z.power_alpha = 5.2944
PML_z.power_kappa = 15.2722
syst.zPML = [PML_z]

syst.epsilon_low = 1
syst.epsilon_high = 1

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both polarizations
input.side = "both"
input.polarization = "both"
# Output to both sides with both polarizations
output.side = "both"
output.polarization = "both"
(S, channels, info)= mesti2s(syst, input, output, opts)

# Compare the scattering matrix to the analytic prediction.
channels_kz = channels.low.kzdx_all[channels.low.ind_prop]
S_ana = zeros(ComplexF64,52,52)
for i=1:13
    S_ana[26+i,i] =  exp(1im*channels_kz[i]*nz)
    S_ana[26+13+i,13+i] =  exp(1im*channels_kz[i]*nz)
    S_ana[i,26+i] =  exp(1im*channels_kz[i]*nz)
    S_ana[13+i,26+13+i] =  exp(1im*channels_kz[i]*nz)
end

pcolormesh(abs.(S))
xlabel("Input channels", fontsize=16)
ylabel("Output channels", fontsize=16)
colorbar()
savefig("s_matrix_in_homo")

pcolormesh(abs.(S-S_ana))
title("\$|S^{mesti.jl}-S^{analytic}|\$")
xlabel("Input channels", fontsize=16)
ylabel("Output channels", fontsize=16)
colorbar()
savefig("s_matrix_diff_in_homo")