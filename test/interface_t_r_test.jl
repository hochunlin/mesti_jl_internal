ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0_par_metis/lib"
# Include MESTI module
include("mesti.jl")
using .MESTI
using LinearAlgebra, SparseArrays, Test
using MAT

# Define the size of matrix and the number of RHS
syst = Syst()
syst.yBC = "periodic"
resolution = 10000
syst.wavelength = 1
syst.dx = 1/resolution
syst.zPML = [PML(100)]

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both s-polarization and p-polarization
input.side = "both"
# Output to both sides with both s-polarization and p-polarization
output.side = "both"
k0dx = 2*pi/syst.wavelength*syst.dx

# Test the functionality in a test set
@testset "reflection coefficients and transmission coefficient:" begin
    for i ∈ 1:10
        n1 = 1+3*rand()
        n2 = 1+3*rand()
	syst.epsilon_low = n1^2
	syst.epsilon_high = n2^2
        syst.epsilon_xx = ones(1,0)
 
	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.

        (S, channels, _)= mesti2s(syst, input, output, opts)
	r = (exp(1im*k0dx*n1/2)-exp(1im*k0dx*n2/2))/(exp(1im*k0dx*n2/2)-exp(-1im*k0dx*n1/2))/exp(1im*k0dx*n1/2)
        t = (sqrt(sin(channels.high.kzdx_prop[1]))/sqrt(sin(channels.low.kzdx_prop[1])))*(exp(1im*k0dx*n1/2)-exp(-1im*k0dx*n1/2))/(exp(1im*k0dx*n2/2)-exp(-1im*k0dx*n1/2))/(exp(1im*k0dx*n1/4)*exp(-1im*k0dx*n2/4))
        @test abs(S[1,1]-r) ≤ 1e-3
	@test abs(S[2,1]-t) ≤ 1e-3
    end
end
