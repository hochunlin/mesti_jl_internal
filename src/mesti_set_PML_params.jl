"""
    SET_PML_PARAMS sets default values for PML parameters
    Refer to Z. Wang et al (in preparation).
"""
function mesti_set_PML_params(pml::Vector{PML}, k0dx::Union{Real,Complex}, epsilon_bg::Union{Vector{Int64},Vector{Float64}}, direction::String)
    
    wavelength_over_dx = ((2*pi)/k0dx)/sqrt(epsilon_bg)

    if ~(length(pml) == 2)
        throw(ArgumentError("$(direction)PML must be a vector of PML containing two elements."))
    end

    # No PML on the both sides
    if (pml[1].npixels == 0 && pml[2].npixels == 0)
        return pml
    end

    # loop over PML parameters on two sides
    for ii = 1:2
        # Number of PML pixels
        if ~isdefined(pml[ii], :npixels)
            throw(ArgumentError("$(direction)PML[$(ii)] must contain field npixels."))
        else
            temp = pml[ii].npixels
            if temp < 0
                throw(ArgumentError("$(direction)PML[$(ii)].npixels must be a non-negative scalar."))             
            end
            if temp == 0
                continue
            end
        end

        # Power of polynomial grading for the conductivity sigma
        if ~isdefined(pml[ii], :power_sigma)
            pml[ii].power_sigma = 3
        else
            temp = pml[ii].power_sigma
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].power_sigma must be a non-negative scalar."))
            end
        end

        # Conductivity at the end of the PML
        if ~isdefined(pml[ii], :sigma_max_over_omega)
            a = -3.0138; b = 0.9303; c = 1.0128
            pml[ii].sigma_max_over_omega = a + b * wavelength_over_dx^c
        else
            temp = pml[ii].sigma_max_over_omega            
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].sigma_max_over_omega must be a non-negative scalar."))
            end
        end

        # Maximal coordinate stretching factor
        if ~isdefined(pml[ii], :kappa_max)
            a = -2.0944; b = 0.6617; c = 1.0467
            pml[ii].kappa_max = a + b * wavelength_over_dx^c
        else
            temp = pml[ii].kappa_max
            if temp < 1
                throw(ArgumentError("$(direction)pml[$ii].kappa_max must be a real scalar that equals or is larger than 1."))
            end
        end

        # Power of polynomial grading for the coordinate stretching factor kappa
        if ~isdefined(pml[ii], :power_kappa)
            pml[ii].power_kappa = 3
        else
            temp = pml[ii].power_kappa
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].power_kappa must be a non-negative scalar."))    
            end
        end

        # Maximal alpha factor for complex frequency shifting (CFS)
        if ~isdefined(pml[ii], :alpha_max_over_omega)
            # CFS is meant to suppress reflection of low-frequency components for time-domain simulations; it is not necessary for frequency-domain simulations
            pml[ii].alpha_max_over_omega = 0
        else
            temp = pml[ii].alpha_max_over_omega
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].alpha_max_over_omega must be a non-negative scalar."))
            end
        end

        # Power of polynomial grading for the alpha factor
        if ~isdefined(pml[ii], :power_alpha)
            # not relevant when alpha_max_over_omega = 0
            pml[ii].power_alpha = 1
        else
            temp = pml[ii].power_alpha
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].power_alpha must be a non-negative scalar."))   
            end
        end         

    end
    
    return pml       
end

    