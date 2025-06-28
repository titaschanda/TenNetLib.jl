
#################################################################################

"""
    mutable struct OptimizeParamsTTN
        maxdim::Vector{Int}
        nsweeps::Vector{Int}
        cutoff::Vector{Float64}
        noise::Vector{Float64}
        noisedecay::Vector{Float64}
        disable_noise_after::Vector{Int}
    end

Holds parameters to control the optimzation sweeps of the TTN.
 - `maxdim::Vector{Int}`: Maximum allowed TTN bond/link dimensions at each stages of optimization.
 - `nsweeps::Vector{Int}`: Number of sweeps to be performed at each statges of optimization.
 - `cutoff::Vector{Float64}`: Cutoff for SVD truncation at each stages of optimization.
 - `noise::Vector{Float64}`: Noise level at each stages of optimization.
 - `noisedecay::Vector{Float64}`: Decay of noise level at each states of optimization.
Noise is divided by `noisedecay` after each sweep.
 - `disable_noise_after::Vector{Int}`: Switch of noise after this many sweeps at each
states of optimization.

All these `Vector`s must have same dimension.
"""
mutable struct OptimizeParamsTTN
    maxdim::Vector{Int}
    nsweeps::Vector{Int}
    cutoff::Vector{Float64}
    noise::Vector{Float64}
    noisedecay::Vector{Float64}
    disable_noise_after::Vector{Int}
end

#################################################################################

"""
    Base.copy(params::OptimizeParamsTTN)

Shallow copy OptimizeParamsTTN.
"""
Base.copy(params::OptimizeParamsTTN) = OptimizeParamsTTN(Base.copy(params.maxdim),
                                                         Base.copy(params.nsweeps),
                                                         Base.copy(params.cutoff),
                                                         Base.copy(params.noise),
                                                         Base.copy(params.noisedecay),
                                                         Base.copy(params.disable_noise_after))

#################################################################################

"""
    function OptimizeParamsTTN(;maxdim::Vector{Int}, nsweeps::Vector{Int}, 
                               cutoff::Union{Vector{Float64}, Float64} = Float64_threshold(),
                               noise::Union{Vector{Float64}, Float64, Int} = 0.0,
                               noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0,
                               disable_noise_after::Union{Vector{Int}, Int} = typemax(Int))

Constructor for `OptimizeParamsTTN`. Takes named arguments.
 - `maxdim::Vector{Int}`: Maximum allowed TTN bond/link dimensions at each stages of optimization.
 - `nsweeps::Vector{Int}`: Number of sweeps to be performed at each statges of optimization.
 - `cutoff::Union{Float64, Vector{Float64}} = Float64_threshold()`: Cutoff for SVD truncation
   at each stages of optimization. If `Float64`, `cutoff` remains same throughout the
   optimization simulation. 
 - `noise::Union{Float64, Int, Vector{Float64}} = 0.0`: Noise level at each stages of
   optimization. If `Float64` / `Int`, initial `noise` remains same throughout the optimization.
 - `noisedecay::Union{Float64, Int, Vector{Float64}} = 1.0`: Decay of noise level at each states
   of optimization. Noise is divided by `noisedecay` after each sweep.
   If `Float64` / `Int`, `noisedecay` remains same throughout the DMRG simulation.
 - `disable_noise_after::Union{Int, Vector{Int}} = typemax(Int)`: Switch of noise after this
   many sweeps at each states of optimization. If `Int`, `disable_noise_after` remains same
   throughout the optimization.
"""
function OptimizeParamsTTN(;maxdim::Vector{Int}, nsweeps::Vector{Int}, 
                           cutoff::Union{Vector{Float64}, Float64} = Float64_threshold(),
                           noise::Union{Vector{Float64}, Float64, Int} = 0.0,
                           noisedecay::Union{Vector{Float64}, Float64, Int} = 1.0,
                           disable_noise_after::Union{Vector{Int}, Int} = typemax(Int))
    
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
        error(string("`OptimizeParamsTTN()` :: Size mismatch in input vectors !! \n",
                     "Lengths of `maxdim` and `nsweeps` must be same !!\n",            
                     "`cutoff`, `noise`, `noisedecay`, `disable_noise_after` ",
                     "can be either scalar or vector having length = `length(maxdim)` !!"))
    end
    return OptimizeParamsTTN(maxdim, nsweeps, cutoffVec,
                             noiseVec, noisedecayVec, disable_noise_afterVec)
end

#################################################################################

"""
    function optimize!(sysenv::StateEnvsTTN,
                       params::OptimizeParamsTTN,
                       sweeppath::Vector{Int2};
                       kwargs...)

Performs optimization of the TTN.

#### Arguments:
 - `sysenv::StateEnvsTTN`.
 - `params::OptimizationParamsTTN`.
 - `sweeppath::Vector{Int2}`: The path to be followed during optimization sweep. A vector
   that must contain all the nodes atleast once.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize after update.
 - `svd_alg::String = "divide_and_conquer"`.
 - `weight::Float64 = -1.0`: Weight for the excited state calculation.
   Must be set to greater than `0` for the excited state optimization.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.

#### Convergence criteria:
 - `energyErrGoal`: Optimization (at a particluar stage) stops when energy difference between
   two consecutive sweeps falls below this threshold and the optimzation moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.

#### Named arguments for `eig_solver` and their default values:
See documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 3`.
 - `solver_maxiter::Int = 1`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `SweepDataTTN`
"""
function optimize!(sysenv::StateEnvsTTN,
                   params::OptimizeParamsTTN,
                   sweeppath::Vector{Int2};
                   kwargs...)::SweepDataTTN
    
    outputlevel = get(kwargs, :outputlevel, 1)
    enerrgoal = get(kwargs, :energyErrGoal, nothing)
    
    swdata = SweepDataTTN()
    
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
            @printf("TTN Optimization level=%d => maxdim=%d, nsweeps=%d, cutoff=%0.2E\n",
                    ii, maxdim, params.nsweeps[ii], cutoff)   
            @printf("TTN Optimization level=%d => noise=%0.2E, noisedecay=%0.3f, disable_noise_after=%d\n",
                    ii, noise, noisedecay, disable_noise_after)
            @printf("-----------------------------------------------------------------------------------\n")
            flush(stdout)
        end

        
        
        for jj=1:params.nsweeps[ii]  

            enerr = fullsweep!(sysenv, sweeppath, eig_solver, swdata;
                               maxdim=maxdim,
                               cutoff=cutoff,
                               noise=noise,
                               kwargs...)
            
            if !isnothing(enerrgoal)  
                errGoalMet = (abs(enerr) < abs(enerrgoal))
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
                        @printf("TTN Optimization Finished !! \n")
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

"""
    function optimize(psi0::TTN, H::CouplingModel,
                      params::OptimizeParamsTTN,
                      sweeppath::Vector{Int2};
                      kwargs...)

    function optimize(psi0::TTN, H::CouplingModel, Ms::Vector{TTN},
                      params::OptimizeParamsTTN,
                      sweeppath::Vector{Int2};
                      kwargs...)

Performs optimization of the TTN.

#### Arguments:
 - `psi0::TTN`: Initial TTN.
 - `H::CouplingModel`, `Ms::Vector{MPS}`.
 - `params::OptimizationParamsTTN`.
 - `sweeppath::Vector{Int2}`: The path to be followed during optimization sweep. A vector
   that must contain all the nodes atleast once.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize after update.
 - `svd_alg::String = "divide_and_conquer"`.
 - `weight::Float64 = -1.0`: Weight for the excited state calculation.
   Must be set to greater than `0` for the excited state optimization.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.

#### Convergence criteria:
 - `energyErrGoal`: Optimization (at a particluar stage) stops when energy difference between
   two consecutive sweeps falls below this threshold and the optimzation moves to the next stage.
   `noise` must be below `Float64_threshold()` to trigger this early stopping.

#### Named arguments for `eig_solver` and their default values:
See documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 3`.
 - `solver_maxiter::Int = 1`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `::Float64`: Energy.
 - `::TTN`: The state psi.
"""
function optimize(psi0::TTN, H::CouplingModel,
                  params::OptimizeParamsTTN,
                  sweeppath::Vector{Int2};
                  kwargs...)::Tuple{Float64, TTN}

    sysenv = StateEnvsTTN(psi0, H)
    swdata = optimize!(sysenv, params, sweeppath; kwargs...)
    return swdata.energy[end], sysenv.psi
end

function optimize(psi0::TTN, H::CouplingModel, Ms::Vector{TTN},
                  params::OptimizeParamsTTN,
                  sweeppath::Vector{Int2};
                  kwargs...)::Tuple{Float64, TTN}
    weight::Float64 = get(kwargs, :weight, -1.0)
    sysenv = StateEnvsTTN(psi0, H, Ms; weight)
    swdata = optimize!(sysenv, params, sweeppath; kwargs...)
    return swdata.energy[end], sysenv.psi
end

function optimize(psi0::TTN, H::CouplingModel, Ms::Vector{TTN{T}},
                  params::OptimizeParamsTTN,
                  sweeppath::Vector{Int2};
                  kwargs...)::Tuple{Float64, TTN} where T
    weight::Float64 = get(kwargs, :weight, -1.0)
    sysenv = StateEnvsTTN(psi0, H, Ms; weight)
    swdata = optimize!(sysenv, params, sweeppath; kwargs...)
    return swdata.energy[end], sysenv.psi
end

#################################################################################
