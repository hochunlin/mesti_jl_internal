# Open Channel Through Disorder

In this example, we show how to use mesti2s() to compute the transmission matrix of a strongly scattering disordered medium, analyze the transmission matrix to determine an incident wavefront that can penetrate the disorder with almost 100% transmission (called an "open channel"), and then use mesti2s() again to compute the field profile of the open channel while comparing to that of a typical plane-wave input.

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
W       = 30    # width of the scattering region
L       = 12    # thickness of the scattering region
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
(epsilon, y0_list, z0_list, r0_list, y_Ex, z_Ex) = 
build_epsilon_disorder_wo_subpixel_smoothing(W, L, r_min, 
r_max, min_sep, number_density, rng_seed, dx,
epsilon_scat, epsilon_bg, build_TM)
```

# Compute the transmission matrix 

```julia
syst = Syst()
syst.epsilon_xx = epsilon
syst.epsilon_low = epsilon_low
syst.epsilon_high = epsilon_high
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC

# specify the input and output
input = channel_type()
output = channel_type()
# input from the low side
input.side = "low"
# output to the high side
output.side = "high"

# put PML along z-direction
syst.zPML = [PML(10)] 

# transmission matrix: input from the low side, output to the high side
t, channels, _ = mesti2s(syst, input, output)
```
```text:Output
===System size===
ny_Ex = 600; nz_Ex = 240 => 262 for Ex(y,z)
[N_prop_low, N_prop_high] = [61, 61] per polarization
yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   4.042 secs
            ... elapsed time:   0.969 secs
Building A  ... elapsed time:   2.729 secs
< Method: APF using MUMPS in single precision with AMD ordering (symmetric K) >
Building K  ... elapsed time:   0.515 secs
Analyzing   ... elapsed time:   0.076 secs
Factorizing ... elapsed time:   0.688 secs
          Total elapsed time:  14.522 secs
```

# Compare an open channel and a plane-wave input

```julia
# The most-open channels is the singular vector of the transmission matrix with 
# the largest singular value.
(_, sigma_max, v_max), _, _, _, _ = svds(t, nsv=1)

N_prop_low = channels.low.N_prop # number of propagating channels on the low side
ind_normal = Int(round((N_prop_low+1)/2)) # index of the normal-incident plane-wave

# compare the transmission
T_avg = sum(abs.(t).^2)/N_prop_low # average over all channels
T_PW  = sum(abs.(t[:,ind_normal]).^2) # normal-incident plane-wave
T_open = sigma_max[1].^2 # open channel

println(" T_avg  = $(T_avg) \n T_PW   = $(T_PW)\n T_open = $(T_open)")
```
```text:Output
 T_avg  = 0.10750856038630088
 T_PW   = 0.07327834665524958
 T_open = 0.9866000088460921
```

```julia
# specify two input incident wavefronts:
# (1) normal-incident plane-wave
# (2) open channel
input = wavefront()
input.v_low = zeros(ComplexF64, N_prop_low, 2)
input.v_low[ind_normal, 1] = 1
input.v_low[:, 2] = v_max

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
ny_Ex = 600; nz_Ex = 240 => 262 for Ex(y,z)
[N_prop_low, N_prop_high] = [61, 61] per polarization
yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   0.373 secs
            ... elapsed time:   0.405 secs
Building A  ... elapsed time:   0.509 secs
< Method: factorize_and_solve using MUMPS in single precision with AMD ordering >
Analyzing   ... elapsed time:   0.056 secs
Factorizing ... elapsed time:   0.792 secs
Solving     ... elapsed time:   0.377 secs
            ... elapsed time:   1.579 secs
          Total elapsed time:   5.200 secs
```

# Animate the field profiles

```julia
# normalize the field amplitude with respect to the plane-wave-input profile
Ex = Ex/maximum(abs.(Ex[:,:,1]))

nframes_per_period = 20

# extend the x coordinate to include free spaces on the two sides
z_Ex = vcat(z_Ex[1] .- (opts.nz_low:-1:1)*dx, z_Ex, z_Ex[end] .+ (1:opts.nz_high)*dx)

# animate the field profile with plane-wave input
anim_pw = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,1]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1)))
    theta = range(0, stop=2π, length=100)

    for jj in 1:length(z0_list)
        y0 = y0_list[jj]
        z0 = z0_list[jj]
        r0 = r0_list[jj]

        y0_circle = y0 .+ r0 * sin.(theta)
        z0_circle = z0 .+ r0 * cos.(theta)

        plot!(plt1, z0_circle, y0_circle, lw=1, color=:black, legend=false)
    end

    display(plot(plt1))    
end
gif(anim_pw, "disorder_PW_input.gif.gif", fps = 10)
```

![disorder_PW_input.gif](disorder_PW_input.gif)

```julia
# animate the field profile of the open channel
anim_open_ch = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,2]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1)))
    theta = range(0, stop=2π, length=100)

    for jj in 1:length(z0_list)
        y0 = y0_list[jj]
        z0 = z0_list[jj]
        r0 = r0_list[jj]

        y0_circle = y0 .+ r0 * sin.(theta)
        z0_circle = z0 .+ r0 * cos.(theta)

        plot!(plt1, z0_circle, y0_circle, lw=1, color=:black, legend=false)
    end

    display(plot(plt1))    
end
gif(anim_open_ch, "disorder_open_channel.gif", fps = 10)
```

![disorder_open_channel.gif](disorder_open_channel.gif)