# This code test the reflection coefficient r and transmission coefficient t of a single interface system with normal incidence.
# We compare the numerical result and the analytic restuls.

syst = Syst()
syst.yBC = "periodic"
resolution = 1000
syst.wavelength = 1
syst.dx = 1/resolution
syst.zPML = [PML(100)] # Use very thick PML to reduce error

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both s-polarization and p-polarization
input.side = "both"
# Output to both sides with both s-polarization and p-polarization
output.side = "both"
k0dx = 2*pi/syst.wavelength*syst.dx # Dimensionless frequency k0*dx

# Test the functionality in a test set
@testset "reflection coefficients and transmission coefficient:" begin
    for i ∈ 1:50
	# Make an interface between two materials whose refractive indices are n1 and n2
        n1 = 1+3*rand() # n1 is a random number between 1 and 4
        n2 = 1+3*rand() # n2 is a random number between 1 and 4
	syst.epsilon_low = n1^2
	syst.epsilon_high = n2^2
        syst.epsilon_xx = ones(1,0)
 
	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.

        (S, channels, _)= mesti2s(syst, input, output, opts) # S = [r,t;t,r]
	# Analytic expression for reflection coefficient r in normal incidence
	r = (exp(1im*k0dx*n1/2)-exp(1im*k0dx*n2/2))/(exp(1im*k0dx*n2/2)-exp(-1im*k0dx*n1/2))/exp(1im*k0dx*n1/2) 
	# Analytic expression for transmission coefficient t in normal incidence
	t = (sqrt(sin(channels.high.kzdx_prop[1]))/sqrt(sin(channels.low.kzdx_prop[1])))*(exp(1im*k0dx*n1/2)-exp(-1im*k0dx*n1/2))/(exp(1im*k0dx*n2/2)-exp(-1im*k0dx*n1/2))/(exp(1im*k0dx*n1/4)*exp(-1im*k0dx*n2/4))
        @test abs(S[1,1]-r) ≤ 1e-3 # Test the correctness of the solution
	@test abs(S[2,1]-t) ≤ 1e-3 # Test the correctness of the solution
    end
end
