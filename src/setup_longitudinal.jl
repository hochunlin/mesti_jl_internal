###### Update on 20220915

mutable struct Side
    N_prop::Integer
    kzdx_all::Vector{ComplexF64}
    ind_prop::Vector{Int64}
    kxdx_prop::Vector{Float64}
    kydx_prop::Vector{Float64}
    kzdx_prop::Vector{Float64}
    sqrt_nu_prop::Vector{Float64}
    ind_prop_conj::Vector{Int64}
    Side() = new()
end

"""
    SETUP_LONGITUDINAL sets up a structure "Side" for one component in the homogeneous space.  
"""
function setup_longitudinal(k0dx::Union{Float64,ComplexF64}, epsilon_bg::Union{Int64,Float64,ComplexF64}, kxdx_all::Union{StepRangeLen{Float64}, Vector{Float64}}, kydx_all::Union{StepRangeLen{Float64}, Vector{Float64}}, kLambda_x::Union{Int64,Float64,ComplexF64,Nothing}=nothing, kLambda_y::Union{Int64,Float64,ComplexF64,Nothing}=nothing, ind_zero_kx::Union{Int64,Nothing}=nothing, ind_zero_ky::Union{Int64,Nothing}=nothing, use_continuous_dispersion::Bool=false)

    k0dx2_epsilon = (k0dx^2)*epsilon_bg
    nx = size(kxdx_all,1)
    ny = size(kydx_all,1)

    side = Side()

    if ~use_continuous_dispersion
        # use the finite-difference dispersion for homogeneous space: k0dx2_epsilon = 4*sin^2(kxdx/2) + 4*sin^2(kydx/2) + 4*sin^2(kzdx/2)
        # sin_kzdx_over_two_sq = sin^2(kzdx/2)
        sin_kzdx_over_two_sq = reshape(0.25*k0dx2_epsilon .- sin.(kxdx_all/2).^2 .- sin.(reshape(vcat(kydx_all),1,:)/2).^2,nx*ny)
        #sin_kzdx_over_two_sq = reshape(0.25*k0dx2_epsilon .- sin.(repeat(kxdx_all,1,ny)/2).^2 .- sin.(repeat(transpose(kydx_all),nx,1)/2).^2, nx*ny)
        # Dimensionless longitudinal wave number
        # asin(sqrt(z)) has two branch points (at z=0 and z=1) and with the branch cuts going outward on the real-z axis; we will address the branch choice below
        # Note kzdx is only defined up to modulo 2*pi (ie, kxdx is equivalent to kzdx + 2*pi, kzdx - 2*pi, etc) because sin(kzdx) and exp(1i*kzdx) are both invariant under 2pi shifts.
        side.kzdx_all = 2*asin.(sqrt.(Complex.((sin_kzdx_over_two_sq))))
        # Indices of the propagating channels
        # When k0dx2_epsilon is real, these are indicies of the channels with real-valued kzdx
        # When k0dx2_epsilon is complex, these are indices of the channels we consider "propagating-like"; they have complex kzdx with 0 < real(kzdx) < pi. When k0dx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real k0dx2_epsilon.
        side.ind_prop = findall(x-> (real(x) > 0 && real(x) < 1), sin_kzdx_over_two_sq)  
        # Here we address the sign choice of kzdx, namely its branch
        # When k0dx2_epsilon is real, we choose the sign of kzdx such that:
        # 1) 0 < kzdx < pi for propagating channels (where kzdx is real)
        # 2) Im(kzdx) >= 0 for evanescent channels
        # Using the correct sign is important when we build the retarded Green's function of the semi-infinite homogeneous space.
        # The default branch choice of asin(sqrt(z)) returns the correct sign for the most part, except when z > 1. We need to flip the sign of those (which can only occur if k0dx2_epsilon > 6).
        # When k0dx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kzdx will be complex-valued, and the sign we "want" for real(kzdx) and the sign we want for imag(kzdx) may be incompatible.
        # What we do with complex k0dx2_epsilon is that we choose the sign for the (complex-valued) kxdx such that when k0dx2_epsilon is tuned to a real number continuously by fixing Re(k0dx2_epsilon) and varying Im(k0dx2_epsilon), the kzdx we choose continuously becomes the "correct" one at real k0dx2_epsilon without crossing any branch cut. To do so, we rotate the two branch cuts of asin(sqrt(z)) by 90 degrees to the lower part of the complex-z plane (ie, the lower part of the complex-k0dx2_epsilon plane), and we pick the branch with the correct sign when k0dx2_epsilon is real. This is implemented by flipping the sign of kzdx for the ones that require flipping.
        # Note that we will get a discontinuity whenever k0dx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
        # The following few lines implement the "flipping".
        if ~(isa(k0dx2_epsilon, Real)) || (isa(k0dx2_epsilon, Real) && k0dx2_epsilon > 6)
            # Note that when imag(sin_kzdx_over_two_sq)=0, flipping is needed for sin_kzdx_over_two_sq>1 but not needed for sin_kzdx_over_two_sq<0.
            ind_flip = findall(x->(real(x)<0 && imag(x)<0) || (real(x)>1 && imag(x)<=0), sin_kzdx_over_two_sq);  
            side.kzdx_all[ind_flip] = -side.kzdx_all[ind_flip]
        end
    else
        # use the continuous dispersion for homogeneous space: k0dx2_epsilon = kxdx^2 + kydx^2 + kzdx^2
        kzdx2 = reshape(k0dx2_epsilon .- kxdx_all.^2 .- reshape(vcat(kydx_all),1,:).^2,nx*ny)
        # kzdx2 = k0dx2_epsilon .- repeat(kxdx_all,1,ny).^2 .- repeat(transpose(kydx_all),nx,1).^2
        
        side.kxdx_all = sqrt.(kzdx2)
        side.ind_prop = findall(real(kzdx2) > 0)
        if ~isa(k0dx2_epsilon, Real)
            ind_flip = findall(x-> (real(x) < 0 && imag(x) < 0), kzdx2)
            side.kzdx_all[ind_flip] = -side.kzdx_all[ind_flip]
        end
    end        
    # Number of propagating channels
    side.N_prop = length(side.ind_prop)

    # Wave numbers of the propagating channels
    side.kzdx_prop = side.kzdx_all[side.ind_prop]
    side.kxdx_prop = kxdx_all[((side.ind_prop).%nx) .+ nx*((side.ind_prop).%nx .== 0)]
    side.kydx_prop = kydx_all[Int.(ceil.((side.ind_prop)./nx))]
    # Square root of the normalized longitudinal group velocity, sqrt(sin(kzdx)), for the propagating channels
    # nu = sin(kxdx)
    # When k0dx2_epsilon is real, sqrt_nu_prop is also real. When k0dx2_epsilon is complex, sqrt_nu_prop is also complex.
    side.sqrt_nu_prop = sqrt.(sin.(side.kzdx_prop))
    
    side.ind_prop_conj = 1:side.N_prop # Just for now.

    return side
end