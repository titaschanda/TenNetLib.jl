
#################################################################################

"""
    mutable struct DMRGParams
        maxdim::Vector{Int}
        nsweeps::Vector{Int}
        cutoff::Vector{Float64}
        noise::Vector{Float64}
        noisedecay::Vector{Float64}
        disable_noise_after::Vector{Int}
    end

Holds parameters to control DMRG sweeps.
 - `maxdim::Vector{Int}`: Maximum allowed MPS bond/link dimensions at each stages of DMRG.
 - `nsweeps::Vector{Int}`: Number of sweeps to be performed at each statges of DMRG.
 - `cutoff::Vector{Float64}`: Cutoff for SVD truncation at each stages of DMRG.
 - `noise::Vector{Float64}`: Noise level at each stages of DMRG.
 - `noisedecay::Vector{Float64}`: Decay of noise level at each states of DMRG.
   Noise is divided by `noisedecay` after each sweep.
 - `disable_noise_after::Vector{Int}`: Switch of noise after this many sweeps at each
states of DMRG.

All these `Vector`s must have same dimension.
"""
mutable struct DMRGParams
    maxdim::Vector{Int}
    nsweeps::Vector{Int}
    cutoff::Vector{Float64}
    noise::Vector{Float64}
    noisedecay::Vector{Float64}
    disable_noise_after::Vector{Int}
end

#################################################################################

"""
    Base.copy(params::DMRGParams)

Shallow copy DMRGParams.
"""
Base.copy(params::DMRGParams) = DMRGParams(Base.copy(params.maxdim),
                                           Base.copy(params.nsweeps),
                                           Base.copy(params.cutoff),
                                           Base.copy(params.noise),
                                           Base.copy(params.noisedecay),
                                           Base.copy(params.disable_noise_after))

#################################################################################

"""
    function DMRGParams(;maxdim::Vector{Int}, nsweeps::Vector{Int}, 
                        cutoff::Union{Vector{Float64}, Float64, Int} = _Float64_Threshold,
                        noise::Union{Vector{Float64}, Float64, Int} = 0.0,
                        noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0,
                        disable_noise_after::Union{Vector{Int}, Int} = typemax(Int))

Constructor for `DMRGParams`. Takes only named arguments.
 - `maxdim::Vector{Int}`: Maximum allowed MPS bond/link dimensions at each stages of DMRG.
 - `nsweeps::Vector{Int}`: Number of sweeps to be performed at each statges of DMRG.
 - `cutoff::Union{Int, Float64, Vector{Float64}} = Float64_threshold()`: Cutoff for SVD truncation
   at each stages of DMRG. If `Float64`, `cutoff` remains same throughout the DMRG simulation. 
 - `noise::Union{Float64, Int, Vector{Float64}} = 0.0`: Noise level at each stages of DMRG.
   If `Float64` or `Int`, initial `noise` remains same throughout the DMRG simulation.
 - `noisedecay::Union{Float64, Int, Vector{Float64}} = 1.0`: Decay of noise level at each states
   of DMRG. Noise is divided by `noisedecay` after each sweep. If `Float64` or `Int`, `noisedecay`
   remains same throughout the DMRG simulation.
 - `disable_noise_after::Union{Int, Vector{Int}} = typemax(Int)`: Switch of noise after this
   many sweeps at each states of DMRG. If `Int`, `disable_noise_after` remains same throughout the
   DMRG simulation.
"""
function DMRGParams(;maxdim::Vector{Int}, nsweeps::Vector{Int}, 
                    cutoff::Union{Vector{Float64}, Float64, Int} = Float64_threshold(),
                    noise::Union{Vector{Float64}, Float64, Int} = 0.0,
                    noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0,
                    disable_noise_after::Union{Vector{Int}, Int} = -1)
    
    cutoffVec = typeof(cutoff) != Vector{Float64} ?
        fill(Float64(cutoff), length(nsweeps)) : cutoff
    
    noiseVec = typeof(noise) != Vector{Float64} ?
        fill(Float64(noise), length(nsweeps)) : noise
    
    noisedecayVec = typeof(noisedecay) != Vector{Float64} ?
        fill(Float64(noisedecay), length(nsweeps)) : noisedecay
    
    disable_noise_afterVec = typeof(disable_noise_after) == Int ?
        fill(disable_noise_after, length(nsweeps)) : disable_noise_after
        
    sanitycheck = length(nsweeps) == length(maxdim) ==
        length(cutoffVec) == length(noiseVec) == length(noisedecayVec) ==
        length(disable_noise_afterVec)
    
    if !sanitycheck
        error(string("`DMRGParams()` :: Size mismatch in input vectors !! \n",
                     "Lengths of `maxdim` and `nsweeps` must be same !!\n",            
                     "`cutoff`, `noise`, `noisedecay`, `disable_noise_after` ",
                     "can be either scalar or vector having length = `length(maxdim)` !!"))
    end
    return DMRGParams(maxdim, nsweeps, cutoffVec,
                      noiseVec, noisedecayVec, disable_noise_afterVec)
end

#################################################################################

"""
    function dmrg!(sysenv::StateEnvs, params::DMRGParams, nsite::Int; kwargs...)

Performs DMRG.

#### Arguments:
 - `sysenv::StateEnvs`.
 - `params::DMRGParams`.
 - `nsite::Int` of the environment. Either `1` or `2` for one-site or two-site update
   respectively.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize after update.
 - `svd_alg::String = "divide_and_conquer"`.
 - `weight::Float64 = -1.0`: Weight for the excited state DMRG. Must be set to greater than `0`
   for the excited state DMRG.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.

#### Convergence criteria:
 - `energyErrGoal`: DMRG (at a particluar stage) stops when energy difference between
   two consecutive sweeps falls below this threshold and DMRG moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.
 - `entropyErrGoal`: DMRG (at a particluar stage) stops when mid-chain entropy difference
   between two consecutive sweeps falls below this threshold and DMRG moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.
When both `energyErrGoal` and `entropyErrGoal` are given, both conditions must be satisfied
to trigger this early stopping.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 5`.
 - `solver_maxiter::Int = 2`.
 - `solver_outputlevel::Int = 0`.: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `SweepData`
"""
function dmrg!(sysenv::StateEnvs, params::DMRGParams, nsite::Int;
               kwargs...)::SweepData
    
    outputlevel = get(kwargs, :outputlevel, 1)
    enerrgoal = get(kwargs, :energyErrGoal, nothing)
    enterrgoal = get(kwargs, :entropyErrGoal, nothing)
    
    swdata = SweepData()
    
    for ii=1:length(params.nsweeps) 
        errGoalMet = false
        maxdim = params.maxdim[ii]
        cutoff = params.cutoff[ii]
        errGoalMet = false
        noise = params.noise[ii]
        noisedecay = params.noisedecay[ii]
        disable_noise_after = params.disable_noise_after[ii]

        if outputlevel > 0
            @printf("-----------------------------------------------------------------------------------\n")
            @printf("DMRG level=%d => maxdim=%d, nsweeps=%d, cutoff=%0.2E\n",
                     ii, maxdim, params.nsweeps[ii], cutoff)   
            @printf("DMRG level=%d => noise=%0.2E, noisedecay=%0.3f, disable_noise_after=%d\n",
                    ii, noise, noisedecay, disable_noise_after)
            @printf("-----------------------------------------------------------------------------------\n")
            flush(stdout)
        end

        
        for jj=1:params.nsweeps[ii]  
            
            enerr, enterr = fullsweep!(sysenv, eig_solver, nsite, swdata; 
                                       maxdim=maxdim,
                                       cutoff=cutoff,
                                       noise=noise,
                                       kwargs...)
            
            if !isnothing(enerrgoal) ||  !isnothing(enterrgoal)            
                if !isnothing(enerrgoal) &&  !isnothing(enterrgoal)
                    errGoalMet = (abs(enerr) < abs(enerrgoal)) && (abs(enterr) < abs(enterrgoal))
                elseif !isnothing(enerrgoal)
                    errGoalMet = (abs(enerr) < abs(enerrgoal))
                elseif !isnothing(enterrgoal)
                    (abs(enterr) < abs(enterrgoal))
                end
            end
            
            if errGoalMet && abs(noise) < Float64_threshold()
                if outputlevel > 0
                    @printf("-----------------------------------------------------------------------------------\n")
                    @printf("Error goal for maxdim=%d met at sweep %d \n",
                            params.maxdim[ii], swdata.sweepcount)                    
                    if ii != length(params.nsweeps)
                        @printf("Continuing with maxdim=%d \n",
                                params.maxdim[ii+1])
                    else
                        @printf("DMRG Finished !! \n")
                    end
                    @printf("-----------------------------------------------------------------------------------\n")
                    flush(stdout)
                end
                break
            end
            
            if jj == disable_noise_after            
                if outputlevel > 0
                    @printf("-----------------------------------------------------------------------------------\n")
                    @printf("Disabling noise after %d sweeps. Last noise level = %0.2E\n",
                            swdata.sweepcount, noise)
                    @printf("-----------------------------------------------------------------------------------\n")
                    flush(stdout)
                end
                noise = 0.0
            end
            noise /= noisedecay
            noise < 100 * Float64_threshold() && (noise = 0)
        end
    end
    return swdata
end  

#################################################################################

function dmrg(psi0::MPS, H::T, params::DMRGParams, nsite::Int;
              kwargs...)::Tuple{Union{Float64, ComplexF64}, MPS} where T <: Union{MPO, Vector{MPO}, CouplingModel}
    sysenv = StateEnvs(psi0, H)
    swdata = dmrg!(sysenv, params, nsite; kwargs...)
    return swdata.energy[end], sysenv.psi
end

function dmrg(psi0::MPS, H::T, Ms::Vector{MPS}, params::DMRGParams, nsite::Int;
              kwargs...)::Tuple{Union{Float64, ComplexF64}, MPS} where T <: Union{MPO, Vector{MPO}, CouplingModel}
    weight::Float64 = get(kwargs, :weight, -1.0)
    sysenv = StateEnvs(psi0, H, Ms; weight)
    swdata = dmrg!(sysenv, params, nsite; kwargs...)
    return swdata.energy[end], sysenv.psi
end

#################################################################################

"""
    function dmrg2(psi0::MPS, H::MPO, params::DMRGParams; kwargs...)
    function dmrg2(psi0::MPS, H::CouplingModel, params::DMRGParams; kwargs...)
    function dmrg2(psi0::MPS, Hs::Vector{MPO}, params::DMRGParams; kwargs...)
    function dmrg2(psi0::MPS, H::MPO, Ms::Vector{MPS}, params::DMRGParams; kwargs...)
    function dmrg2(psi0::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}, params::DMRGParams; kwargs...)
    function dmrg2(psi0::MPS, H::CouplingModel, Ms::Vector{MPS}, params::DMRGParams; kwargs...)

Performs two-site DMRG.

#### Arguments:
 - `psi0::MPS`: Initial MPS.
 - `H::MPO`, `H::CouplingModel`, `Hs::Vector{MPO}`, `Ms::Vector{MPS}`.
 - `params::DMRGParams`.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize after update.
 - `svd_alg::String = "divide_and_conquer"`.
 - `weight::Float64 = -1.0`: Weight for the excited state DMRG. Must be set to greater than `0`
   for the excited state DMRG.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after every
   fullsweep, if `2` prints at every update step.

#### Convergence criteria:
 - `energyErrGoal`: DMRG (at a particluar stage) stops when energy difference between
   two consecutive sweeps falls below this threshold and DMRG moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.
 - `entropyErrGoal`: DMRG (at a particluar stage) stops when mid-chain entropy difference
   between two consecutive sweeps falls below this threshold and DMRG moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.
When both `energyErrGoal` and `entropyErrGoal` are given, both conditions must be satisfied
to trigger this early stopping.

#### Named arguments for `eig_solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 5`.
 - `solver_maxiter::Int = 2`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `::Float64`: Energy.
 - `::MPS`: The state psi.
"""
dmrg2(psi0::MPS, H::T, params::DMRGParams;
      kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
          dmrg(psi0::MPS, H::T, params::DMRGParams, 2; kwargs...)

dmrg2(psi0::MPS, H::T, Ms::Vector{MPS}, params::DMRGParams;
      kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
          dmrg(psi0::MPS, H::T, Ms::Vector{MPS}, params::DMRGParams, 2; kwargs...)


"""
    function dmrg1(psi0::MPS, H::MPO, params::DMRGParams; kwargs...)
    function dmrg1(psi0::MPS, H::CouplingModel, params::DMRGParams; kwargs...)
    function dmrg1(psi0::MPS, Hs::Vector{MPO}, params::DMRGParams; kwargs...)
    function dmrg1(psi0::MPS, H::MPO, Ms::Vector{MPS}, params::DMRGParams; kwargs...)
    function dmrg1(psi0::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}, params::DMRGParams; kwargs...)
    function dmrg1(psi0::MPS, H::CouplingModel, Ms::Vector{MPS}, params::DMRGParams; kwargs...)

Performs single-site DMRG. All other details are same as in [`dmrg2`](@ref).
"""
dmrg1(psi0::MPS, H::T, params::DMRGParams;
      kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
          dmrg(psi0::MPS, H::T, params::DMRGParams, 1; kwargs...)

dmrg1(psi0::MPS, H::T, Ms::Vector{MPS}, params::DMRGParams;
      kwargs...) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
          dmrg(psi0::MPS, H::T, Ms::Vector{MPS}, params::DMRGParams, 1; kwargs...)

#################################################################################
