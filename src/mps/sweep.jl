
#################################################################################

"""
    mutable struct SweepData
        sweepcount::Int
        maxchi::Vector{Int}
        energy::Vector{Float64}
        entropy::Vector{Float64}
        maxtruncerr::Vector{Float64}
        lasteigs::Vector{Vector{Float64}}
    end

Holds historical data after each (full)sweep. Requires for convergence check etc.
 - `sweepcount::Int`: Number of fullsweeps performed.
 - `maxchi::Vector{Int}`: Maximum MPS bond/link dimensions after every sweep.
 - `energy::Vector{Float64}`: Energies after every sweep.
 - `entropy::Vector{Float64}`: Mid-chain entropies after every sweep.
 - `maxtrucerr::Vector{Float64}`: Maximum truncation error after every sweep.
 - `lasteigs::Vector{Vector{Float64}}`: Spectrum of eigenvalues at each bond after
   previous halfsweep.

#### Default constructor:
 - `SweepData()`: Initialize an empty `SweepData` object.
"""
mutable struct SweepData
    sweepcount::Int
    maxchi::Vector{Int}
    energy::Vector{Float64}
    entropy::Vector{Float64}
    maxtruncerr::Vector{Float64}
    lasteigs::Vector{Vector{Float64}}
end

SweepData() = SweepData(0, Int[], Float64[], Float64[], Float64[], Float64[])

#################################################################################
"""
    Base.copy(swdata::SweepData)

Shallow copy of `SweepData`.
"""
Base.copy(swdata::SweepData) = SweepData(swdata.sweepcount,
                                         Base.copy(swdata.maxchi),
                                         Base.copy(swdata.energy),
                                         Base.copy(swdata.entropy),
                                         Base.copy(swdata.mastruncerr),
                                         Base.copy(swdata.lasteigs)
                                         )

#################################################################################

"""
    function fullsweep!(sysenv::StateEnvs, solver, nsite::Int, swdata::SweepData;
                        kwargs...)

Perform a fullsweep (left-to-right and right-to-left) by `solver`.

#### Arguments:
 - `sysenv::StateEnvs`.
 - `solver`: Solver for update. Available ones: `eig_solver` and `exp_solver`.
 - `nsite::Int` of the environment. Either `1` or `2` for one-site or two-site update
    respectively.
 - `swdata::SweepData`.

#### Named arguments and their default values:
 - `time_step::Union{Float64, ComplexF64, Nothing} = nothing`: Time step for TDVP. 
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `noise::Float64 = 0.0`.
 - `reverse_step::Bool = false` if `time_step = nothing`, `true` otherwise.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-14` if `eig_solver`, `1E-12` if `exp_solver`.
 - `solver_krylovdim::Int = 5` if `eig_solver`, `30` if `exp_solver`.
 - `solver_maxiter::Int = 2` if `eig_solver`, `100` if `exp_solver`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false` if `eig_solver`, `true` if `exp_solver`.
 - `solver_check_convergence::Bool = false` if `eig_solver`, `true` if `exp_solver`.

#### Return values:
 - `::Float64`: Change in Energy ΔE
 - `::Float64`: Change in Entropy ΔS

`swdata::SweepData` gets updated.
"""
function fullsweep!(sysenv::StateEnvs, solver, nsite::Int, swdata::SweepData;
                    kwargs...)::Tuple{Float64, Float64}
    
    outputlevel::Int = get(kwargs, :outputlevel, 1)
    noise = get(kwargs, :noise, 0.0)
    
    if !isortho(sysenv.psi) || orthocenter(sysenv.psi) != 1
        orthogonalize!(sysenv.psi, 1)
    end
    
    @assert isortho(sysenv.psi) && orthocenter(sysenv.psi) == 1
    
    energy::Float64 = NaN
    maxtruncerr = 0.0
    swdata.sweepcount += 1   
    N = length(sysenv)
    lasteigs = Vector{Vector{Float64}}(undef, N-1)
    
    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end
    
    sw_time = @elapsed begin
        for bond = 1 : N-1
            energy, err, eigs = update_position!(sysenv, solver, bond, nsite, "left"; kwargs...)
            # DO NOT NEED HERE... 
            # add left eigs
            # lasteigs[bond] = eigs
            #
            maxtruncerr = max(err, maxtruncerr)    
            if (outputlevel > 1)
                @printf("At left sweep %d bond %d => Energy %s, Err %.2E \n",
                        swdata.sweepcount, bond, energy, err)
                flush(stdout)
            end
            if nsite == 1 && bond == N-1
                energy, _, _ = update_position!(sysenv, solver, bond+1, nsite, "left"; kwargs...)
            end 
        end
        # DO NOT NEED HERE... 
        #swdata.lasteigs = lasteigs

        for bond = (N-1) : -1 : 1
            site = nsite == 1 ? bond+1 : bond            
            energy, err, eigs = update_position!(sysenv, solver, site, nsite, "right"; kwargs...)
            # add right eigs
            lasteigs[bond] = eigs
            #
            maxtruncerr = max(err, maxtruncerr)
            if (outputlevel > 1)
                @printf("At right sweep %d bond %d => Energy %s, Err %.2E \n",
                        swdata.sweepcount, bond, energy, err)
                flush(stdout)
            end
            if nsite == 1 && bond == 1
                energy, _, _ = update_position!(sysenv, solver, bond, nsite, "right"; kwargs...)
            end 
        end
        swdata.lasteigs = lasteigs
    end
    
    push!(swdata.maxchi, maxlinkdim(sysenv.psi))
    push!(swdata.energy, energy)
    mideigs = lasteigs[N ÷ 2]
    push!(swdata.entropy, _entropy(mideigs / sum(mideigs)))
    push!(swdata.maxtruncerr, maxtruncerr)
    swdata.lasteigs = lasteigs
    
    if (swdata.sweepcount > 1)
        enerr = swdata.energy[end] - swdata.energy[end-1]
        enterr = swdata.entropy[end] - swdata.entropy[end-1]
    else
        enerr = NaN
        enterr = NaN
    end
    
    if (outputlevel > 0)
        @printf("-----------------------------------------------------------------------------------\n")
        @printf("At sweep %d => E=%s, S=%s, MaxLinkDim=%d, Noise=%0.2E\n",
                swdata.sweepcount, swdata.energy[end], swdata.entropy[end], swdata.maxchi[end],
                noise)
        @printf("At sweep %d => ΔE=%.2E, ΔS=%.2E, MaxErr=%.2E, Time=%0.3f\n",
                swdata.sweepcount, enerr, enterr, maxtruncerr, sw_time)
        @printf("-----------------------------------------------------------------------------------\n")
        flush(stdout)
    end
    
    
    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end

    return enerr, enterr
end

#################################################################################

"""
    function dynamic_fullsweep!(sysenv::StateEnvs, solver, swdata::SweepData;
                                kwargs...)

Perform a dynamic fullsweep (left-to-right and right-to-left) by `solver`.
The very first sweep, as dictated by `swdata.sweepcount=0`, Global Subspace Expansion
(see below) is performed followed by a pure one-site sweep if
`typeof(sysenv) == StateEnvs{ProjMPO}`, else performs a full two-site sweep. At each bond,
if the lowest eigenvalue is below `eigthreshold` or the bond dimension at that bond
has reached `maxdim` at a particular halfsweep, performs single-site update across
that bond in the subsequent halfsweep, otherwise performs two-site update.

#### Arguments:
 - `sysenv::StateEnvs`.
 - `solver`: Solver for update. Available ones: `eig_solver` and `exp_solver`.
 - `swdata::SweepData`.

#### Named arguments and their default values:
 - `time_step::Union{Float64, ComplexF64, Nothing} = nothing`: Time step for TDVP. 
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `noise::Float64 = 0.0`.
 - `reverse_step::Bool = false` if `time_step = nothing`, `true` otherwise.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.
 - `eigthreshold::Float64 = 1E-12`.
 - `extendat::Union{Nothing, Int} = nothing`: If specified, at every `extendat`th sweep,
   Global Subspace Expansion is performed followed by a pure one-site sweep if
   `typeof(sysenv) == StateEnvs{ProjMPO}`, else performs a full two-site sweep.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-14` if `eig_solver`, `1E-12` if `exp_solver`.
 - `solver_krylovdim::Int = 3` if `eig_solver`, `30` if `exp_solver`.
 - `solver_maxiter::Int = 1` if `eig_solver`, `100` if `exp_solver`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false` if `eig_solver`, `true` if `exp_solver`.
 - `solver_check_convergence::Bool = false` if `eig_solver`, `true` if `exp_solver`.


#### Arguments for Global Subspace Expansion and their default values:
 - `extension_krylovdim::Int = 3`: Number of Krylov vectors used for GSE.
 - `extension_applyH_cutoff::Float64 = 0.0`: Cutoff for the application of the MPO to the MPS.
 - `extension_applyH_maxdim::Int = maxlinkdim(psi) + 1`: Maximum bond/link dimension
   for the application of the MPO to the MPS.
 - `extension_cutoff::Float64 = 1E-10`: Cutoff for the basis extension step in GSE.

#### Return values:
 - `::Float64`: Change in Energy ΔE
 - `::Float64`: Change in Entropy ΔS

`swdata::SweepData` gets updated.
"""
function dynamic_fullsweep!(sysenv::StateEnvs, solver, swdata::SweepData;
                            eigthreshold::Float64 = 1E-12,
                            extendat::Union{Nothing, Int} = nothing,
                            kwargs...)::Tuple{Float64, Float64}

    maxdim::Int = get(kwargs, :maxdim, typemax(Int))    
    outputlevel::Int = get(kwargs, :outputlevel, 1)
    noise = get(kwargs, :noise, 0.0)
        
    if (swdata.sweepcount == 0)
        if typeof(sysenv) == StateEnvs{ProjMPO}
            krylov_extend!(sysenv; kwargs...)
            return fullsweep!(sysenv, solver, 1, swdata; kwargs...)
        else
            return fullsweep!(sysenv, solver, 2, swdata; kwargs...)
        end
    end
    
    if (!isnothing(extendat) && (swdata.sweepcount+1) % extendat == 0)
        if typeof(sysenv) == StateEnvs{ProjMPO}
            krylov_extend!(sysenv; kwargs...)
            return fullsweep!(sysenv, solver, 1, swdata; kwargs...)
        else
            return fullsweep!(sysenv, solver, 2, swdata; kwargs...)
        end
    end

    if !isortho(sysenv.psi) || orthocenter(sysenv.psi) != 1
        orthogonalize!(sysenv.psi, 1)
    end
    
    @assert isortho(sysenv.psi) && orthocenter(sysenv.psi) == 1
    
    energy = 0.0
    maxtruncerr = 0.0
    swdata.sweepcount += 1   
    N = length(sysenv)
    
    lasteigs = Vector{Vector{Float64}}(undef, N-1)
    
    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end
    
    sw_time = @elapsed begin
        for bond = 1 : N-1
            nsite = (swdata.lasteigs[bond][end] < eigthreshold
                     ||  dim(linkind(sysenv.psi, bond)) >= maxdim) ? 1 : 2
            energy, err, eigs = update_position!(sysenv, solver, bond, nsite, "left"; kwargs...)
            # add left eigs
            lasteigs[bond] = eigs
            #
            maxtruncerr = max(err, maxtruncerr)    
            if (outputlevel > 1)
                @printf("At left sweep %d bond %d => Energy %s, Err %.2E \n",
                        swdata.sweepcount, bond, energy, err)
                flush(stdout)
            end
            if nsite == 1 && bond == N-1
                energy, _, _ = update_position!(sysenv, solver, bond+1, nsite, "left"; kwargs...)
            end            
        end
        swdata.lasteigs = lasteigs
        
        for bond = (N-1) : -1 : 1
            nsite = (swdata.lasteigs[bond][end] < eigthreshold
                     ||  dim(linkind(sysenv.psi, bond)) >= maxdim) ? 1 : 2
            site = nsite == 1 ? bond+1 : bond            
            energy, err, eigs = update_position!(sysenv, solver, site, nsite, "right"; kwargs...)
            # add right eigs
            lasteigs[bond] = eigs
            #
            maxtruncerr = max(err, maxtruncerr)
            if (outputlevel > 1)
                @printf("At right sweep %d bond %d => Energy %s, Err %.2E \n",
                        swdata.sweepcount, bond, energy, err)
                flush(stdout)
            end
            if nsite == 1 && bond == 1
                energy, _, _ = update_position!(sysenv, solver, bond, nsite, "right"; kwargs...)
            end
        end
        swdata.lasteigs = lasteigs
    end
    
    push!(swdata.maxchi, maxlinkdim(sysenv.psi))
    push!(swdata.energy, energy)
    mideigs = lasteigs[N ÷ 2]
    push!(swdata.entropy, _entropy(mideigs / sum(mideigs)))
    push!(swdata.maxtruncerr, maxtruncerr)
    swdata.lasteigs = lasteigs
    
    if (swdata.sweepcount > 1)
        enerr = swdata.energy[end] - swdata.energy[end-1]
        enterr = swdata.entropy[end] - swdata.entropy[end-1]
    else
        enerr = NaN
        enterr = NaN
    end
    
    if (outputlevel > 0)
        @printf("-----------------------------------------------------------------------------------\n")
        @printf("At sweep %d => E=%s, S=%s, MaxLinkDim=%d, Noise=%0.2E\n",
                swdata.sweepcount, swdata.energy[end], swdata.entropy[end], swdata.maxchi[end],
                noise)
        @printf("At sweep %d => ΔE=%.2E, ΔS=%.2E, MaxErr=%.2E, Time=%0.3f\n",
                swdata.sweepcount, enerr, enterr, maxtruncerr, sw_time)
        @printf("-----------------------------------------------------------------------------------\n")
        flush(stdout)
    end
    
    
    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end

    return enerr, enterr
end

#################################################################################

"""
    function krylov_extend!(psi::MPS, H::MPO; kwargs...)

Performs Global Subspace Expansion.

#### Named arguments and their default values:
 - `extension_krylovdim::Int = 3`: Number of Krylov vectors used for GSE.
 - `extension_applyH_cutoff::Float64 = Float64_threshold()`: Cutoff for the application
   the MPO to the MPS.
 - `extension_applyH_maxdim::Int = maxlinkdim(psi) + 1`: Maximum bond/link
   dimension of the resulting MPS after the application of the MPO to the MPS.
 - `extension_cutoff::Float64 = 1E-10`: Cutoff for the basis extension step in GSE.
"""
function krylov_extend!(psi::MPS, H::MPO; kwargs...)

    extension_krylovdim::Int = get(kwargs, :extension_krylovdim, 3)
    extension_applyH_cutoff::Float64 = get(kwargs, :extension_applyH_cutoff,
                                           Float64_threshold())
    extension_applyH_maxdim::Int = get(kwargs, :extension_applyH_maxdim,
                                       maxlinkdim(psi) + 1)
    extension_cutoff::Float64 = get(kwargs, :extension_cutoff,
                                    1E-10)
    
    phis = Vector{MPS}(undef, extension_krylovdim)
        
    for k in 1 : extension_krylovdim
        prev = k == 1 ? psi : phis[k - 1]
        phis[k] = apply(H, prev;
                        maxdim = extension_applyH_maxdim,
                        cutoff = extension_applyH_cutoff)
        normalize!(phis[k])
    end
    _krylov_addbasis!(psi, phis, extension_cutoff)
    return nothing
end
    
#################################################################################

"""
    function krylov_extend!(sysenv::StateEnvs{ProjMPO}; kwargs...)

Performs Global Subspace Expansion. The `StateEnvs` must be created by a single MPO.

#### Named arguments and their default values:
 - `extension_krylovdim::Int = 3`: Number of Krylov vectors used for GSE.
 - `extension_applyH_cutoff::Float64 = Float64_threshold()`: Cutoff for the application
   the MPO to the MPS.
 - `extension_applyH_maxdim::Int = maxlinkdim(psi) + 1`: Maximum bond/link
   dimension of the resulting MPS after the application of the MPO to the MPS.
 - `extension_cutoff::Float64 = 1E-10`: Cutoff for the basis extension step in GSE.
"""
function krylov_extend!(sysenv::StateEnvs{ProjMPO}; kwargs...)::Nothing

    extension_krylovdim::Int = get(kwargs, :extension_krylovdim, 3)
    extension_applyH_cutoff::Float64 = get(kwargs, :extension_applyH_cutoff,
                                           Float64_threshold())
    extension_applyH_maxdim::Int = get(kwargs, :extension_applyH_maxdim,
                                       maxlinkdim(sysenv.psi) + 1)
    extension_cutoff::Float64 = get(kwargs, :extension_cutoff,
                                    1E-10)
    
    kr_time = @elapsed begin

        krylov_extend!(sysenv.psi, sysenv.PH.H; kwargs...)
        
        sysenv.PH.lpos = 0
        sysenv.PH.rpos = length(sysenv.PH.H) + 1
        sysenv.PH.nsite = 2
        sysenv.PH.LR = Vector{ITensor}(undef, length(sysenv.PH.H))
    end
    
    outputlevel::Int = get(kwargs, :outputlevel, 1)
    if (outputlevel > 0)
        @printf("-----------------------------------------------------------------------------------\n")
        @printf("Global Subspace Expansion: KrylovDim=%d, applyH Cutoff=%.2E, applyH MaxDim=%d\n",
                extension_krylovdim, extension_applyH_cutoff, extension_applyH_maxdim)
        @printf("Global Subspace Expansion: Cutoff=%.2E, MaxLinkDim=%d, Time=%0.3f\n",
                extension_cutoff, maxlinkdim(sysenv.psi), kr_time)
        @printf("-----------------------------------------------------------------------------------\n")
        flush(stdout)
    end
    
    return nothing
end

#################################################################################
"""
Copied from https://github.com/ITensor/ITensorTDVP.jl/pull/24/files#diff-d7e5dafa1fb67b61b1e780001c8066ef3f01409e70877cd8647d62ab692c1169.
"""
function _krylov_addbasis!(psi::MPS, phis::Vector{MPS}, extension_cutoff::Float64)::MPS
    N = length(psi)

    orthogonalize!(psi, N)
    for phi in phis
        orthogonalize!(phi, N)
    end

    s = siteinds(psi)

    for j in reverse(2:N)
        # SVD psi[j] to compute B
        linds = (s[j - 1], linkind(psi, j - 1))
        _, S, B = svd(psi[j], linds...; righttags="Link,l=$j")
        rinds = uniqueinds(B, S)

        # Make projector
        Id = ITensor(1.0)
        for r in rinds
            Id *= denseblocks(delta(r', dag(r)))
        end
        P = Id - prime(B, rinds) * dag(B)

        # Sum phi density matrices
        rho = ITensor()
        for phi in phis
            rho += prime(phi[j], rinds) * dag(phi[j])
        end
        rho /= tr(rho)

        # Apply projector
        PrhoP = apply(apply(P, rho), P)

        if norm(PrhoP) > 1E-14
            # Diagonalize projected density matrix PrhoP
            # to compute Bphi, which spans part of right basis 
            # of phis which is orthogonal to right basis of psi
            D, Bphi = eigen(PrhoP;
                            cutoff=extension_cutoff,
                            ishermitian=true,
                            righttags="bϕ_$j,Link")

            ## Test Bphi is ortho to B
            #O = Bphi*B
            #if norm(O) > 1E-10
            #  @show norm(O)
            #  error("Non-zero overlap of extended basis with original basis")
            #end

            # Form direct sum of B and Bphi over left index
            bψ = commonind(B, S)
            bϕ = commonind(Bphi, D)
            
            #if !hasqns(bψ)
            #    bx = Index(dim(bψ) + dim(bϕ), "Link,l=$(j-1)")
            #else
            #    bx = Index(vcat(space(bψ), space(bϕ)), "Link,l=$(j-1)")
            #end
            #D1, D2 = directsum_itensors(bψ, bϕ, dag(bx))
            #Bx = D1 * B + D2 * Bphi

            # TODO: test
            Bx, _ = directsum(B => bψ, Bphi => bϕ; tags = "Link,l=$(j-1)")
        else
            Bx = B
        end

        # Shift ortho center one site left using dag(Bx)
        # and replace tensor at site j with Bx
        psi[j - 1] = psi[j - 1] * (psi[j] * dag(Bx))
        psi[j] = Bx
        for phi in phis
            phi[j - 1] = phi[j - 1] * (phi[j] * dag(Bx))
            phi[j] = Bx
        end
    end
    orthogonalize!(psi, 1)    
    return psi
end

#################################################################################
