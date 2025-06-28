
#################################################################################

"""
    halfsweep_done(N::Int, pos::Int, nsite::Int, ortho::String)::Bool

Check wheather half sweep (left-to-right or right-to-left) is completed.
 - `N::Int`: Length of the MPS.
 - `pos::Int`: Current position.
 - `nsite` of the environment. Either `1` or `2` for one-site or two-site update respectively.
 - `ortho::String`: Direction of the sweep. Either `"left"` or `"right"`.
"""
function halfsweep_done(N::Int, pos::Int, nsite::Int, ortho::String)::Bool
    if pos == 1 && ortho == "right"
        return true
    elseif pos == N && ortho == "left" && nsite == 1
        return true
    elseif pos == N-1 && ortho == "left" && nsite == 2
        return true
    else
        return false
    end
end

#################################################################################

function _update_two_site!(sysenv::StateEnvs, solver, pos::Int, ortho::String,
                           time_step::Union{Float64, ComplexF64, Nothing},
                           normalize::Bool,                
                           maxdim::Int,
                           mindim::Int,
                           cutoff::Float64,
                           svd_alg::String,
                           noise::Float64,
                           reverse_step::Bool;                           
                           kwargs...)::Tuple{Float64, Float64, Vector{Float64}}
    # Sanity check
    @assert pos > 0 && pos < length(sysenv)
    @assert (orthocenter(sysenv.psi) == pos && ortho == "left") || 
        (orthocenter(sysenv.psi) == pos+1 && ortho == "right")
    
    nsite = 2
    set_nsite!(sysenv, nsite)
    phi = sysenv.psi[pos] * sysenv.psi[pos+1]
    position!(sysenv, pos)        
    energy, phi = solver(sysenv, phi, time_step; kwargs...)
    normalize && normalize!(phi)
    isnan(energy) && (energy = real(scalar(dag(phi) * sysenv.PH(phi))))

    drho = nothing
    if abs(noise) > Float64_threshold()
        drho = noise * noiseterm(sysenv.PH, phi, ortho)
    end
    
    spec = replacebond!(
        sysenv.psi,
        pos,
        phi;
        maxdim=maxdim,
        mindim=mindim,
        cutoff=cutoff,
        eigen_perturbation=drho,
        ortho=ortho,
        normalize=normalize,
        which_decomp=nothing,
        svd_alg=svd_alg
    )
    
    if reverse_step && !halfsweep_done(length(sysenv), pos, nsite, ortho)
        pos1 = ortho == "left" ? pos+1 : pos
        phi0 = sysenv.psi[pos1]       
        set_nsite!(sysenv, nsite - 1)
        position!(sysenv, pos1)
        energy, phi0 = solver(sysenv, phi0, -time_step; kwargs...)
        normalize && normalize!(phi0)
        isnan(energy) && (energy = real(scalar(dag(phi0) * sysenv.PH(phi0))))
        sysenv.psi[pos1] = phi0
    end
    
    return energy, spec.truncerr, spec.eigs
end

#################################################################################

function _update_one_site!(sysenv::StateEnvs, solver, pos::Int, ortho::String,
                           time_step::Union{Float64, ComplexF64, Nothing},
                           normalize::Bool,                
                           maxdim::Int,
                           mindim::Int,
                           cutoff::Float64,
                           svd_alg::String,
                           noise::Float64,
                           reverse_step::Bool;
                           kwargs...)::Tuple{Float64, Float64, Vector{Float64}}

    # Sanity check
    @assert pos > 0 && pos <= length(sysenv)
    @assert orthocenter(sysenv.psi) == pos
    
    nsite = 1
    set_nsite!(sysenv, nsite)
    phi = sysenv.psi[pos]
    position!(sysenv, pos)        
    energy, phi = solver(sysenv, phi, time_step; kwargs...)
    normalize && normalize!(phi)
    isnan(energy) && (energy = real(scalar(dag(phi) * sysenv.PH(phi))))

    truncerr = 0.0
    eigs::Vector{Float64} = Float64[]         

    if halfsweep_done(length(sysenv), pos, nsite, ortho)
        sysenv.psi[pos] = phi
    else   
        posnext = ortho == "left" ? pos+1 : pos-1
        pos0 = ortho == "left" ? pos : pos-1

        if abs(noise) > Float64_threshold()
            phi *= sysenv.psi[posnext]
            set_nsite!(sysenv, nsite+1)
            position!(sysenv, pos0)
            drho = noise * noiseterm(sysenv.PH, phi, ortho)

            spec = replacebond!(
                sysenv.psi,
                pos0,
                phi;
                maxdim=maxdim,
                mindim=mindim,
                cutoff=cutoff,
                eigen_perturbation=drho,
                ortho=ortho,
                normalize=normalize,
                which_decomp=nothing,
                svd_alg=svd_alg
            )
            truncerr = spec.truncerr
            eigs = spec.eigs
            set_nsite!(sysenv, nsite)

        else            
            uinds = uniqueinds(phi, sysenv.psi[posnext])
            origtag = tags(commonind(sysenv.psi[pos], sysenv.psi[posnext]; tags = "Link"))    

            U, S, V, spec = svd(
                phi, uinds; 
                maxdim=maxdim,
                mindim=mindim,
                cutoff=cutoff,
                alg = svd_alg,
                lefttags = origtag)
            
            normalize && normalize!(S)
            sysenv.psi[pos] = U
            phi0 = S*V        
            ortho == "left" && setleftlim!(sysenv.psi, pos)
            ortho == "right" && setrightlim!(sysenv.psi, pos)            
            truncerr = spec.truncerr
            eigs = spec.eigs
            
            if reverse_step
                pos1 = ortho == "left" ? pos + 1 : pos
                set_nsite!(sysenv, nsite - 1)
                position!(sysenv, pos1)
                energy, phi0 = solver(sysenv, phi0, -time_step; kwargs...)
                normalize && normalize!(phi0)
                isnan(energy) && (energy = real(scalar(dag(phi0) * sysenv.PH(phi0))))
            end
            sysenv.psi[posnext] = phi0 * sysenv.psi[posnext]
        end
    end
    return energy, truncerr, eigs
end

#################################################################################

"""
    update_position!(sysenv::StateEnvs, solver, pos::Int, nsite::Int, ortho::String; kwargs...)

Updates StateEnvs at position `pos` by `solver`.

#### Arguments:
 - `sysenv::StateEnvs`
 - `solver`: Solver for update. Available ones: `eig_solver` and `exp_solver`.
 - `pos::Int`: Position of the bond (`nsite=2`) or site (`nsite=1`).
 - `nsite` of the environment. Either `1` or `2` for one-site or two-site update respectively.
 - `ortho::String`: Direction of the sweep. Either `"left"` or `"right"`.

#### Named arguments and their default values:
 - `time_step::Union{Float64, ComplexF64, Nothing} = nothing`: Time step for TDVP. 
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `noise::Float64 = 0.0`.
 - `reverse_step::Bool = false` if `time_step = nothing`, `true` otherwise.

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
 - `::Float64`: Energy.
 - `::Float64`: Truncation Error.
 - `::Vector{Float64}`: SVD spectrum.
"""
function update_position!(sysenv::StateEnvs, solver,
                          pos::Int, nsite::Int, ortho::String; 
                          kwargs...)::Tuple{Float64, Float64, Vector{Float64}}

    time_step::Union{Float64, ComplexF64, Nothing} = get(kwargs, :time_step, nothing)
    normalize::Bool = get(kwargs, :normalize, true)
    
    # SVD/Decomp Algs    
    maxdim::Int = get(kwargs, :maxdim, typemax(Int))
    mindim::Int = get(kwargs, :mindim, 1)
    cutoff::Float64 = get(kwargs, :cutoff, Float64_threshold())
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")    
    
    noise::Float64 = get(kwargs, :noise, 0.0)
    reverse_step::Bool = get(kwargs, :reverse_step, isnothing(time_step) ? false : true)
    
    if noise > 0 && reverse_step
        error(string("`updatePosition()` :: `noise=$noise` cannot be greater than zero",
                     " for `reverse_step=$reverse_step` !!"))
    end
    
    if nsite == 2
        return _update_two_site!(sysenv, solver, pos, ortho,
                                 time_step,
                                 normalize,
                                 maxdim,
                                 mindim,
                                 cutoff,
                                 svd_alg,
                                 noise,
                                 reverse_step;
                                 kwargs...)
    elseif nsite == 1
        return _update_one_site!(sysenv, solver, pos, ortho,
                                 time_step,
                                 normalize,
                                 maxdim,
                                 mindim,
                                 cutoff,
                                 svd_alg,
                                 noise,
                                 reverse_step;
                                 kwargs...)
    else
        error("`update_position()` with `nsite=$nsite` not implemented !!")
    end    
end

#################################################################################

#=

function _expandterm(P::ProjMPO, phi::ITensor, 
                     ortho::String)::Tuple{ITensor, Union{Index, Nothing},
                                           Int, Union{Index, Nothing}}
    
    itensor_map = Union{ITensor,OneITensor}[]
    push!(itensor_map, ortho == "left" ? lproj(P) : rproj(P))
    append!(itensor_map, P.H[site_range(P)])    
    Mix = phi
    for it in itensor_map
        Mix *= it
    end
    noprime!(Mix)
    
    b = ortho == "left" ? site_range(P)[end] : site_range(P)[begin]
    
    if (b == 1 && ortho == "right") || (b == length(P) && ortho == "left")
        return Mix, nothing, b, nothing
    end
    
    enlinkind = commonind(phi, ortho == "left" ? rproj(P) : lproj(P); tags = "Link")
    enlinkindH = commonind(P.H[b], ortho == "left" ? rproj(P) : lproj(P); tags = "Link")
    tag = addtags(tags(enlinkind), "Subspace")
    comb = combiner(enlinkind, enlinkindH; tags=tag)
    mixind = combinedind(comb)
    Mix *= comb
    return Mix, mixind, b, enlinkind
end

#################################################################################

_expandterm(P::ProjMPO_MPS2, phi::ITensor, ortho::String) = _expandterm(P.PH, phi, ortho)

#################################################################################

function _expandterm(P::ProjMPOSum, phi::ITensor, 
                     ortho::String)::Tuple{ITensor, Union{Index, Nothing},
                                           Int, Union{Index, Nothing}}
    
    Mixers = [_expandterm(x, phi, ortho) for x in P.terms]    
    
    if isnothing(Mixers[begin][2])
        return sum(x -> x[begin], Mixers), nothing, Mixers[begin][3], nothing
    end
    
    tag = tags(Mixers[begin][2])
    pairmix = [x[1] => x[2] for x in Mixers]
    S, s = directsum(pairmix...; tags=tag)
    return S, s, Mixers[begin][3], Mixers[begin][4]
end

#################################################################################

"""
    subspace_expand!(sysenv::StateEnvs, phi::ITensor, ortho::String, noise::Float64)

Returns subspace expanded `phi`. Also pads next site with zeros.
"""
function subspace_expand!(sysenv::StateEnvs, 
                          phi::ITensor,
                          ortho::String, 
                          noise::Float64)::ITensor
    
    Mixer, mixind, b,  enlinkind = _expandterm(sysenv.PH, phi, ortho)
    Mixer *= noise
    
    if isnothing(mixind)
        return phi + Mixer
    end
    
    bnext = ortho == "left" ? b+1 : b-1
    
    indzero = [dag(mixind), siteind(sysenv.psi, bnext)]
    linkindzero = uniqueind(sysenv.psi[bnext], sysenv.psi[b]; tags="Link")
    !isnothing(linkindzero) && push!(indzero, linkindzero)
    zerotensor = diagITensor(indzero)
    
    phienlarged, sumind = directsum(phi => enlinkind, Mixer => mixind, tags=tags(enlinkind))
    sysenv.psi[bnext], sumindtemp = directsum(sysenv.psi[bnext] => enlinkind, 
                                              zerotensor => mixind)

    ## ADD COMBINERS
    replaceind!(sysenv.psi[bnext], sumindtemp, dag(sumind))
    ortho == "left" && setrightlim!(sysenv.psi, b)
    ortho == "right" && setleftlim!(sysenv.psi, b)
    return phienlarged
end

#################################################################################

function _updateOneSite456!(sysenv::StateEnvs, solver, pos::Int, ortho::String,
                            time_step::Union{Float64, ComplexF64, Nothing},
                            normalize::Bool,                
                            maxdim::Int,
                            mindim::Int,
                            cutoff::Float64,
                            svd_alg::String,
                            noise::Float64,
                            reverse_step::Bool;
                            kwargs...)::Tuple{Float64, Float64, Vector{Float64}}

    # Sanity check
    @assert pos > 0 && pos <= length(sysenv)
    @assert orthocenter(sysenv.psi) == pos
    
    nsite = 1
    set_nsite!(sysenv, nsite)
    phi = sysenv.psi[pos]
    position!(sysenv, pos)        
    ts = @elapsed begin
        energy, phi = solver(sysenv.PH, phi, time_step; kwargs...)
    end
    global Tsol += ts
    normalize && normalize!(phi)
    isnan(energy) && (energy = real(scalar(dag(phi) * sysenv.PH(phi))))

    truncerr = 0.0
    eigs::Vector{Float64} = []         

    if halfsweep_done(length(sysenv), pos, nsite, ortho)
        sysenv.psi[pos] = phi
    else   
        posnext = ortho == "left" ? pos+1 : pos-1
        uinds = uniqueinds(phi, sysenv.psi[posnext])
        origtag = tags(commonind(sysenv.psi[pos], sysenv.psi[posnext]; tags = "Link"))
        
  	if noise > 0.0
            tp = @elapsed begin
       	        phi = subspace_expand!(sysenv, phi, ortho, noise::Float64)
            end
            global Tpad += tp
        end

        tsvd = @elapsed begin
            U, S, V, spec = svd(
                phi, uinds; 
                maxdim=maxdim,
                mindim=mindim,
                cutoff=cutoff,
                alg = svd_alg,
                lefttags = origtag)
        end
        global Tsvd += tsvd
        
        normalize && normalize!(S)
        sysenv.psi[pos] = U
        phi0 = S*V        
        ortho == "left" && setleftlim!(sysenv.psi, pos)
        ortho == "right" && setrightlim!(sysenv.psi, pos)            
        truncerr = spec.truncerr
        eigs = spec.eigs
        
        if reverse_step
            # Required as next site gets updated by noise
            noise > 0.0 && position!(sysenv, posnext)
            #
            pos1 = ortho == "left" ? pos + 1 : pos
            set_nsite!(sysenv, nsite - 1)
            position!(sysenv, pos1)
            energy, phi0 = solver(sysenv.PH, phi0, -time_step; kwargs...)
            normalize && normalize!(phi0)
            isnan(energy) && (energy = real(scalar(dag(phi0) * sysenv.PH(phi0))))
        end
        sysenv.psi[posnext] = phi0 * sysenv.psi[posnext]
    end

    if pos == 1
        println(Tsol, " ", Tpad, " ", Tsvd)
        global Tsol = 0
        global Tpad = 0
        global Tsvd = 0
    end
    
    return energy, truncerr, eigs
end

#################################################################################

=#

