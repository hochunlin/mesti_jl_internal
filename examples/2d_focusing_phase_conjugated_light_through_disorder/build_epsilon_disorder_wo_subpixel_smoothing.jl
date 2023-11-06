using Random

function build_epsilon_disorder_wo_subpixel_smoothing(W, L, r_min, r_max, min_sep, number_density, rng_seed, dx, epsilon_scat, epsilon_bg, build_TM, build_TE = false, yBC = "periodic", y1 = 0, y2 = W, z1 = 0, z2 = L; no_scatterered_center = true)
    # Generate a random collection of cylinders and build the relative permittivity profile.
    # The full system spans y in [0, W], z in [0, L].
    # The scattering region spans [z1, z2], [y1, y2].
    # Note that subpixel smoothing is not used.
    #   
    #   === Input Arguments ===
    #   W       = width of the full system
    #   L       = thickness of the full system
    #   r_min   = minimal radius of the cylinders
    #   r_max   = maximal radius of the cylinders
    #   min_sep = minimal separation between cylinders
    #   number_density = number density of the cylinders
    #   rng_seed = random number generator seed
    #   dx       = discretization grid size
    #   epsilon_scat = relative permittivity of the cylinders
    #   epsilon_bg   = relative permittivity of the background
    #   build_TM = whether to build epsilon for TM polarization
    #   build_TE = whether to build 1/epsilon for TE polarization
    #   yBC = boundary condition in y (required only when build_TE = true)
    #   y1 (optional) = bottom edge of the scattering region (defaults to 0)
    #   y2 (optional) = top edge of the scattering region (defaults to W)
    #   z1 (optional) = left edge of the scattering region (defaults to 0)
    #   z2 (optional) = right edge of the scattering region (defaults to L)
    #   
    #   === Output Arguments ===
    #   epsilon_xx (ny_Ex-by-nz_Ex matrix): discretized epsilon_zz
    #   inv_epsilon (2-element cell array):
    #      inv_epsilon_yy (ny_Ex-by-nz_Hx matrix): discretized (1/epsilon)_zz
    #      inv_epsilon_zz (ny_Hx-by-nz_Ex matrix): discretized (1/epsilon)_yy
    #   y0_list (column vector): y coordinates of the cylinder centers
    #   z0_list (column vector): z coordinates of the cylinder centers
    #   r0_list (column vector): radii of the cylinders
    #   y_Ex (row vector): y coordinates of the centers of the Ex pizels
    #   z_Ex (row vector): z coordinates of the centers of the Ex pixels
    #   y_Hx (row vector): y coordinates of the centers of the Hx pixels
    #   z_Hx (row vector): z coordinates of the centers of the Hx pixels

    # round the system size
    ny_Ex = round(Int, W/dx)
    nz_Ex = round(Int, L/dx)
    W = ny_Ex * dx
    L = nz_Ex * dx

    #=
    if nargin() < 12
        error("Not enough input arguments.")
    elseif nargin() == 12 && build_TE
        error("Input argument 'yBC' must be given when build_TE = true.")
    elseif (nargin() == 12 && ~build_TE) || nargin() == 13
        # default scattering region
        z1 = 0
        z2 = L
        y1 = 0
        y2 = W
    elseif nargin() != 17
        error("Number of input arguments must be 12, 13, or 17.")
    elseif z1 < 0 || z2 > L || y1 < 0 || y2 > W
        error("The ranges [z1, z2] and [y1, y2] are invalid.")
    end
    =#

    if y1 < 0 || y2 > W || z1 < 0 || z2 > L
        error("The ranges [y1, y2] and [z1, z2] are invalid.")
    end

    N_scatterers = round(Int, number_density * (y2 - y1) * (z2 - z1))
    if no_scatterered_center
        N_scatterers = N_scatterers + 1
    end
    
    # location of the Ex grid points
    y_Ex = (0.5:ny_Ex) * dx
    z_Ex = (0.5:nz_Ex) * dx

    if build_TM
        epsilon_xx = epsilon_bg*ones(ny_Ex, nz_Ex)
    else
        epsilon_xx = nothing
    end

    if build_TE
        nz_Hx = nz_Ex + 1
        ny_Hx, dm_Hx = 0, 0

        if ismember(lowercase(yBC), ["bloch", "periodic", "pmcpec"])
            ny_Hx = ny_Ex
        elseif lowercase(yBC) == "pec"
            ny_Hx = ny_Ex + 1
            dm_Hx = 1
        elseif lowercase(yBC) == "pmc"
            ny_Hx = ny_Ex - 1
        elseif lowercase(yBC) == "pecpmc"
            ny_Hx = ny_Ex
            dm_Hx = 1
        else
            error("yBC = '$yBC' is not a supported option.")
        end

         = 1/epsilon_bg*ones(ny_Hx, nz_Ex)
        inv_epsilon_zz = 1/epsilon_bg*ones(ny_Ex, nz_Hx)

        # location of the Hx grid points
        y_Hx = ((1:ny_Hx) .- dm_Hx) * dx
        z_Hx = (0:(nz_Hx-1)) * dx
    else
        inv_epsilon_yy = nothing
        inv_epsilon_zz = nothing        
        y_Hx = nothing
        z_Hx = nothing
    end

    #=if N_scatterers > 3000
        # store the indices of the existing cylinders by their approximate locations,
        # to facilitate overlap checking when N_scatterers is large
        use_cell = true
        cell_size = 2 * r_max + min_sep
        nz_cell = ceil((z2 - z1) / cell_size)
        ny_cell = ceil((y2 - y1) / cell_size)
        ind_scatterers = [Vector{Int}[] for _ in 1:ny_cell, _ in 1:nz_cell]
    else
    =#
        use_cell = false
    #end

    # pick the radii and locations of N_scatterers non-overlapping cylinders,
    # and then generate the discretized permittivity profile
    Random.seed!(rng_seed)  # set random number generator seed
    y0_list, z0_list, r0_list = zeros(N_scatterers), zeros(N_scatterers), zeros(N_scatterers)

    for ii in 1:N_scatterers
        # randomly pick the radius of the new cylinder
        r0 = r_min + rand() * (r_max - r_min)

        if use_cell
            # distance for overlap checking
            dist_check = r0 + min_sep + r_max
        end

        # min/max values for the coordinates of the cylinder center
        y_min = y1 + r0
        y_max = y2 - r0
        z_min = z1 + r0
        z_max = z2 - r0

        # keep randomly picking the cylinder location until we get one that doesn't overlap with the existing cylinders
        found = false

        while !found
            global y0
            global z0
            if no_scatterered_center && ii == 1
                y0 = (y1+y2)/2
                z0 = (z1+z2)/2
                found = true
            else
                # center of the new cylinder
                y0 = y_min + rand() * (y_max - y_min)
                z0 = z_min + rand() * (z_max - z_min)

                #=
                if use_cell
                    # indices of cells we need to check for overlap
                    l1 = max(1, round(Int, 0.5 + (z0 - z1 - dist_check) / cell_size))
                    l2 = min(nz_cell, round(Int, 0.5 + (z0 - z1 + dist_check) / cell_size))
                    m1 = max(1, round(Int, 0.5 + (y0 - y1 - dist_check) / cell_size))
                    m2 = min(ny_cell, round(Int, 0.5 + (y0 - y1 + dist_check) / cell_size))

                    # indices of cylinders we need to check for overlap
                    ind_check = vcat(ind_scatterers[m, n] for m in m1:m2, n in l1:l2)

                    # check those cylinders for overlap
                    if isempty(ind_check) || minimum((z0_list[ind] - z0).^2 + (y0_list[ind] - y0).^2 - (r0_list[ind] .+ (r0 + min_sep)).^2) > 0
                        found = true
                    end
                else
                =#
                    # check all existing cylinders for overlap
                    if ii == 1 || minimum((y0_list[1:(ii-1)] .- y0).^2 .+ (z0_list[1:(ii-1)] .- z0).^2 .- (r0_list[1:(ii-1)] .+ (r0 + min_sep)).^2) > 0
                        found = true
                    end
                #end
            end
        end
        
        y0_list[ii] = y0
        z0_list[ii] = z0
        r0_list[ii] = r0

        if use_cell
            # store the index of this cylinder by its approximate location
            m0 = min(ny_cell, max(1, round(Int, 0.5 + (y0 - y1) / cell_size)))
            l0 = min(nz_cell, max(1, round(Int, 0.5 + (z0 - z1) / cell_size)))
            push!(ind_scatterers[m0, l0], ii)
        end

        # set the relative permittivity of the cylinder
        dn = r0 / dx

        if build_TM && ~(no_scatterered_center && ii == 1)
            # identify indices of the interiors of the cylinder within a local box surrounding (z0, y0) with edge length = 2 * r0
            # note: we don't do subpixel smoothing here for simplicity
            # Ex and epsilon_xx are located at (z, y) = (n - 0.5, m - 0.5) * dx with n from 1 to nz_Ex, m from 1 to ny_Ex
            m0 = 0.5 + y0 / dx # location of (y0, z0) in terms of index
            l0 = 0.5 + z0 / dx  
            m1 = max(1, round(Int, m0 - dn))
            m2 = min(ny_Ex, round(Int, m0 + dn))
            l1 = max(1, round(Int, l0 - dn))
            l2 = min(nz_Ex, round(Int, l0 + dn))
            y_local = ((m1 - 0.5):m2) * dx
            z_local = ((l1 - 0.5):l2) * dx
            #Z_local, Y_local = meshgrid(z_local, y_local)
            Y_local = repeat(y_local, 1, length(z_local))                                    
            Z_local = repeat(z_local', length(y_local))
            indices = findall(((Y_local .- y0).^2 + (Z_local .- z0).^2 ) .<= (r0^2))
            m_local = [ci[1] for ci in indices]
            l_local = [ci[2] for ci in indices]
            ind = Base._sub2ind((ny_Ex, nz_Ex), m_local .+ (m1 - 1), l_local .+ (l1 - 1))  # convert to linear indices in the full matrix
            epsilon_xx[ind] .= epsilon_scat
        end

        if build_TE && ~(no_scatterered_center && ii == 1)

            # (1/epsilon)_zz = (1/epsilon)_zz is located at (y, z) = (m - 0.5, l - 1) * dx with n from 1 to nz_Hx = nz_Ex + 1, m from 1 to ny_Ex
            m0 = 0.5 + y0 / dx # location of (z0, y0) in terms of index
            l0 = 1 + z0 / dx 
            m1 = max(1, round(Int, m0 - dn))
            m2 = min(ny_Ex, round(Int, m0 + dn))
            l1 = max(1, round(Int, l0 - dn))
            l2 = min(nz_Hx, round(Int, l0 + dn))
            y_local = ((m1 - 0.5):m2) * dx
            z_local = ((l1:l2) .- 1) * dx
            #Z_local, Y_local = meshgrid(z_local, y_local)
            Y_local = repeat(y_local, 1, length(z_local))                                    
            Z_local = repeat(z_local', length(y_local))
            indices = findall(((Y_local .- y0).^2 + (Z_local .- z0).^2 ) .<= (r0^2))
            m_local = [ci[1] for ci in indices]
            l_local = [ci[2] for ci in indices]
            ind = Base._sub2ind((ny_Ex, nz_Hx), m_local .+ (m1 - 1), l_local .+ (l1 - 1))  # convert to linear indices in the full matrix
            inv_epsilon_zz[ind] .= 1 / epsilon_scat

            # note: inv_epsilon[1][:,1] and inv_epsilon[1][:,end] are half outside the scattering region,
            # so they need to be set separately later based on epsilon_L and epsilon_R.

            # (1/epsilon)_yy is located at z = (l - 0.5) * dx with n from 1 to nz_Ex, and
            # y = m * dx for periodic, Bloch periodic, PMC, and PMCPEC boundary conditions in y
            # y = (m - 1) * dx for PEC and PECPMC boundary conditions in y
            # with m from 1 to ny_Hx
            m0 = dm_Hx + y0 / dx  # location of (z0, y0) in terms of index
            l0 = 0.5 + z0 / dx 
            m1 = max(1, round(Int, m0 - dn))
            m2 = min(ny_Hx, round(Int, m0 + dn))
            l1 = max(1, round(Int, l0 - dn))
            l2 = min(nz_Ex, round(Int, l0 + dn))
            y_local = ((m1:m2) .- dm_Hx) * dx
            z_local = ((l1 - 0.5):l2) * dx
            #Z_local, Y_local = meshgrid(z_local, y_local)
            Y_local = repeat(y_local, 1, length(z_local))                        
            Z_local = repeat(z_local', length(y_local))
            indices = findall(((Y_local .- y0).^2 + (Z_local .- z0).^2) .<= (r0^2))
            m_local = [ci[1] for ci in indices]
            l_local = [ci[2] for ci in indices]
            ind = Base._sub2ind((ny_Hx, nz_Ex), m_local .+ (m1 - 1), l_local .+ (l1 - 1))  # convert to linear indices in the full matrix
            inv_epsilon_yy[ind] .= 1 / epsilon_scat
        end
    end
    if no_scatterered_center
        y0_list = y0_list[2:end]
        z0_list = z0_list[2:end]
        r0_list = r0_list[2:end]    
    end
    
    if build_TM && build_TE
        return epsilon_xx, inv_epsilon_yy, inv_epsilon_zz, y0_list, z0_list, r0_list,  y_Ex, z_Ex, y_Hx, z_Hx
    end

    if build_TM && ~build_TE
        return epsilon_xx, y0_list, z0_list, r0_list, y_Ex, z_Ex
    end

    if ~build_TM && build_TE
        return inv_epsilon_yy, inv_epsilon_zz, y0_list, z0_list, r0_list, y_Hx, z_Hx
    end

end