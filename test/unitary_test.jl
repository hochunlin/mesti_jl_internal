# This code test the unitarity of scattering matrices S upon random permittivity profiles.
# We check S'*S ≈ I, an identity matrix.

# Specify parameters of the system
syst = Syst()
syst.xBC = "periodic"
syst.yBC = "periodic"
syst.dx = 1
syst.wavelength = 5
syst.epsilon_low = 1
syst.epsilon_high = 1
epsilon_max = 4
epsilon_min = 1

# Define the size of the scattering region (1 wavelength by 1 wavelength by 1 wavelength)
nx, ny, nz = syst.wavelength, syst.wavelength, syst.wavelength
nx_Ex = nx; ny_Ex = ny; nz_Ex = nz -1
nx_Ey = nx; ny_Ey = ny; nz_Ey = nz -1
nx_Ez = nx; ny_Ez = ny; nz_Ez = nz

# Use optimized PML parameters for this resolution to reduce error
#zpml = PML(10) # 5 PML pixels padded on each side of z-direction
#zpml.sigma_max_over_omega = 4.664171353771144
#zpml.power_sigma = 5.568495004043644
#zpml.kappa_max = 5.260462519440901
#zpml.power_kappa= 11.973331766657591
zpml = get_optimal_PML(5)
zpml.npixels = 10
syst.zPML = [zpml]

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
@testset "unitarity:              " begin
    for i ∈ 1:20
	# Random permittivity profiles, whose value is between 1 and 4
        syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min

	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.

        (S, _, _)= mesti2s(syst, input, output, opts)
        @test maximum(abs.((S'*S) - I(size(S, 1)))) ≤ 1e-2 # Check the unitarity of the scattering matrix
    end
end
