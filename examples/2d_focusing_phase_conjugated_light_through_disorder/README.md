# Focusing Phase Conjugated Light Through Disorder

In this example, we show how to use mesti() to compute the field profile of a point source in a scattering disordered medium, do phase conjugation to determine an incident wavefront that can focus on the the disorder, and then use mesti2s() again to compute the field profile to show its focus.

```julia
# Call necessary packages
using MESTI, Arpack, Plots

# Include the function to build epsilon_xx for the disordered
include("build_epsilon_disorder_wo_subpixel_smoothing.jl")
```

# System parameters

```julia
# dimensions of the system, in units of the wavelength lambda_0
dx      = 1/20  # discretization grid size
W       = 100   # width of the scattering region
L       = 20    # thickness of the scattering region
L_tot   = 40    # full length of the system for plotting
r_min   = 0.2   # minimal radius of the cylindrical scatterers
r_max   = 0.4   # maximal radius of the cylindrical scatterers
min_sep = 0.05  # minimal separation between cylinders
number_density = 1.3  # number density, in units of 1/lambda_0^2
rng_seed = 0   # random number generator seed

# relative permittivity, unitless
epsilon_scat = 2.0^2  # cylindrical scatterers
epsilon_bg   = 1.0^2  # background in the scattering region
epsilon_low  = 1.0^2  # frees space on the low side
epsilon_high = 1.0^2  # frees space on the high side

yBC = "periodic" # boundary condition in y

# generate a random collection of non-overlapping cylinders
# note: subpixel smoothing is not applied for simplicity
build_TM = true
no_scatterered_center = true
(epsilon, y0_list, z0_list, r0_list, y_Ex, z_Ex) = 
build_epsilon_disorder_wo_subpixel_smoothing(W, L, r_min, 
r_max, min_sep, number_density, rng_seed, dx,
epsilon_scat, epsilon_bg, build_TM, no_scatterered_center)
```

# Compute the field profile of a point source 

```julia
syst = Syst()
pml_npixels = 10
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC
# specify the permittivity profile of the simulation domain including the low side, scattering region and the high side.
syst.epsilon_xx = cat(epsilon_low*ones(size(epsilon,1),pml_npixels+1), epsilon, epsilon_high*ones(size(epsilon,1),pml_npixels+1), dims=2)

# specify the input (point source in the middle of the disordered)
Bx = Source_struct()
Bx.pos = [[Int((W/dx)/2),Int((L/dx)/2)+pml_npixels+1,1,1]]
Bx.data = [ones(1,1)]

# put PML along z-direction
pml = PML(pml_npixels)
pml.direction = "z"
syst.PML = [pml]

# field profile: input from a point source
Ex_field, _ = mesti(syst, [Bx])
```
```text:Output
===System size===
ny_Ex = 2000; nz_Ex = 422 for Ex(y,z)
UPML on -z +z sides; ; yBC = periodic; zBC = PEC
Building B,C... elapsed time:   0.202 secs
Building A  ... elapsed time:   3.542 secs
< Method: factorize_and_solve using MUMPS in double precision with AMD ordering >
Analyzing   ... elapsed time:   0.394 secs
Factorizing ... elapsed time:   7.489 secs
Solving     ... elapsed time:   0.799 secs
          Total elapsed time:   9.657 secs
```

# Compute the phase-conjugated focusing profile

```julia
# Specify the system for mesti2s() and mesti_build_channels()
syst = Syst()
syst.epsilon_xx = epsilon
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC
syst.epsilon_low = 1
syst.epsilon_high = 1

# build channels for this system
channels = mesti_build_channels(syst)
N_prop_low = channels.low.N_prop # number of propagating channels on the low side
ind_normal = Int(round((N_prop_low+1)/2)) # index of the normal-incident plane-wave

# build projection matrix C on the low side
C_low = channels.low.sqrt_nu_prop.*exp.((-1im*0.5)*channels.low.kzdx_prop).*convert(Matrix, adjoint(channels.u_x_m(channels.low.kydx_prop)))
proj_coefficient = C_low*Ex_field[:,pml_npixels+1]

# specify two input incident wavefronts:
# (1) normal-incident plane-wave
# (2) phase-conjugated wavefront
input = wavefront()
input.v_low = zeros(ComplexF64, N_prop_low, 2)
input.v_low[ind_normal, 1] = 1

# for the phased conjugated input: 
# conj(coefficient*u) = conj(coefficient)*conj(u) = conj(coefficient)*conj(u) =  conj(coefficient)*perm(u(ky)) = perm(conj(coefficient))*u(ky)
# perm() means permute a vector that switches one propagating channel with one
# having a complex-conjugated transverse profile
# for periodic boundary, this flips the sign of ky. 
input.v_low[:, 2] = conj(proj_coefficient)[channels.low.ind_prop_conj]/norm(proj_coefficient)


# we will also get the field profile in the free spaces on the two sides, for
# plotting purpose.
opts = Opts()
opts.nz_low = round((L_tot-L)/2/dx)
opts.nz_high = opts.nz_low

# for field-profile computations
Ex, _, _ = mesti2s(syst, input, opts)
```
```text:Output
===System size===
ny_Ex = 2000; nz_Ex = 400 => 422 for Ex(y,z)
[N_prop_low, N_prop_high] = [201, 201] per polarization
yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   0.612 secs
            ... elapsed time:   0.848 secs
Building A  ... elapsed time:   0.893 secs
< Method: factorize_and_solve using MUMPS in double precision with AMD ordering >
Analyzing   ... elapsed time:   0.378 secs
Factorizing ... elapsed time:   7.951 secs
Solving     ... elapsed time:   0.734 secs
            ... elapsed time:   1.946 secs
          Total elapsed time:  14.772 secs
```

# Animate the field profiles

```julia
# normalize the field amplitude with respect to the phase-conjugated-input profile
Ex = Ex/maximum(abs.(Ex[:,:,2]))

nframes_per_period = 20

# extend the x coordinate to include free spaces on the two sides
z_Ex = vcat(z_Ex[1] .- (opts.nz_low:-1:1)*dx, z_Ex, z_Ex[end] .+ (1:opts.nz_high)*dx)

# animate the field profile with plane-wave input
anim_pw = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,1]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1), aspect_ratio=:equal))
    theta = range(0, stop=2π, length=100)

    for jj in 1:length(z0_list)
        y0 = y0_list[jj]
        z0 = z0_list[jj]
        r0 = r0_list[jj]

        y0_circle = y0 .+ r0 * sin.(theta)
        z0_circle = z0 .+ r0 * cos.(theta)

        plot!(plt1, z0_circle, y0_circle, lw=0.1, color=:black, legend=false)
    end

    display(plot(plt1))    
end
gif(anim_pw, "plane_wave_input.gif.gif", fps = 10)
```

![plane_wave_input.gif](plane_wave_input.gif)

```julia
# animate the field profile of the phase-conjugated focusing
anim_open_ch = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,2]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1), aspect_ratio=:equal))
    theta = range(0, stop=2π, length=100)

    for jj in 1:length(z0_list)
        y0 = y0_list[jj]
        z0 = z0_list[jj]
        r0 = r0_list[jj]

        y0_circle = y0 .+ r0 * sin.(theta)
        z0_circle = z0 .+ r0 * cos.(theta)

        plot!(plt1, z0_circle, y0_circle, lw=0.1, color=:black, legend=false)
    end

    display(plot(plt1))    
end
gif(anim_open_ch, "phase_conjugated_focusing.gif", fps = 10)
```

![phase_conjugated_focusing](phase_conjugated_focusing.gif)
