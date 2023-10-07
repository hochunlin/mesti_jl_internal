ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0_par_metis/lib"
# Include MESTI module
include("mesti.jl")
using .MESTI
using LinearAlgebra, SparseArrays, Test
using MAT

# Define the size of matrix and the number of RHS
nx, ny, nz = 10, 10, 10
nx_Ex = nx; ny_Ex = ny; nz_Ex = nz -1
nx_Ey = nx; ny_Ey = ny; nz_Ey = nz -1
nx_Ez = nx; ny_Ez = ny; nz_Ez = nz 

syst = Syst()
syst.xBC = "periodic"
syst.yBC = "periodic"
syst.dx = 1
syst.wavelength = 5
syst.epsilon_low = 1
syst.epsilon_high = 1
epsilon_max = 4
epsilon_min = 1
syst.zPML = [PML(100)]

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both s-polarization and p-polarization
input.side = "both"
input.polarization = "both"
# Output to both sides with both s-polarization and p-polarization
output.side = "both"
output.polarization = "both"

# Test the functionality in a test set
@testset "unitarity: " begin
    for i ∈ 1:10
        syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min

	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.

        (S, _, _)= mesti2s(syst, input, output, opts)
        @test maximum(abs.((S'*S) - I(size(S, 1)))) ≤ 1e-2 # Test the unitarity of the system
    end
end
