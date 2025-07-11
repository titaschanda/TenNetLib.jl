
#################################################################################

"""
    mutable struct TDVPEngine{T <: Union{ProjMPO,
                                         ProjMPOSum2,
                                         ProjMPO_MPS2,
                                         ProjMPOSum_MPS,
                                         ProjCouplingModel,
                                         ProjCouplingModel_MPS}}

Holds the `MPS` state, `SweepData`, and absolute elpased time for TDVP simulations.
 - `sysenv::StateEnvs`: Holds the state psi and environments PH.
 - `swdata::SweepData`: Holds the historical data after each (full)sweep.
 - `abstime::Float64`: Absolute elapsed time.
"""
mutable struct TDVPEngine{T <: Union{ProjMPO,
                                     ProjMPOSum2,
                                     ProjMPO_MPS2,
                                     ProjMPOSum_MPS,
                                     ProjCouplingModel,
                                     ProjCouplingModel_MPS}}
    sysenv::StateEnvs{T}
    swdata::SweepData
    abstime::Float64
end

#################################################################################

"""
    Base.copy(engine::TDVPEngine)

Shallow copy of `TDVPEngine`.
"""
Base.copy(engine::TDVPEngine) = TDVPEngine(Base.copy(engine.sysenv),
                                           Base.copy(swdata), engine.abstime)

#################################################################################

"""
    function TDVPEngine(psi::MPS, H::{T}) where T <: Union{MPO, Vector{MPO}, CouplingModel}

Constructor of the `TDVPEngine` from different forms of Hamiltonians.
"""
function TDVPEngine(psi::MPS, H::T) where T <: Union{MPO, Vector{MPO}, CouplingModel}
    sysenv = StateEnvs(psi, H)
    swdata = SweepData()
    return TDVPEngine(sysenv, swdata, 0.0)
end

#################################################################################

"""
    function TDVPEngine(psi::MPS, H::T, Ms::Vector{MPS};
                        weight::Float64) where T <: Union{MPO, Vector{MPO}, CouplingModel}

Constructor of `TDVPEngine` from different forms of Hamiltonians and a vector of MPS.
"""
function TDVPEngine(psi::MPS, H::T, Ms::Vector{MPS};
                    weight::Float64) where T <: Union{MPO, Vector{MPO}, CouplingModel}
    sysenv = StateEnvs(psi, H, Ms; weight=weight)
    swdata = SweepData()
    return TDVPEngine(sysenv, swdata, 0.0)
end

#################################################################################

"""
    function getpsi(engine::TDVPEngine)

Returns (shallow copy of) the state psi.
"""
function getpsi(engine::TDVPEngine)
    return Base.copy(engine.sysenv.psi)
end

#################################################################################

"""
    function sweepcount(engine::TDVPEngine)

Returns the number of sweeps performed.
"""
sweepcount(engine::TDVPEngine) = engine.swdata.sweepcount

#################################################################################

"""
    function getenergy(engine::TDVPEngine)

Returns the energy of the state psi.
"""
getenergy(engine::TDVPEngine) = engine.swdata.energy[end]

#################################################################################

"""
    function getentropy(engine::TDVPEngine)

Returns the mid-chain entropy of the state psi.
"""
getentropy(engine::TDVPEngine) = engine.swdata.entropy[end]

#################################################################################

"""
    function maxchi(engine::TDVPEngine)

Returns the maximum bond/link dimension of the state psi.
"""
maxchi(engine::TDVPEngine) = engine.swdata.maxchi[end]

#################################################################################

"""
    function totalerror(engine::TDVPEngine)

Returns the sum of the truncation errors of all the sweeps performed.
"""
totalerror(engine::TDVPEngine) = sum(engine.swdata.maxtruncerr)

#################################################################################

"""
    function sweeperror(engine::TDVPEngine)

Returns the truncation error in the last sweep.
"""
sweeperror(engine::TDVPEngine) = engine.swdata.maxtruncerr[end]

#################################################################################

"""
    function krylov_extend!(engine::TDVPEngine{ProjMPO}; kwargs...)

Performs Global Subspace Expansion.

#### Arguments and their default values:
 - `extension_krylovdim::Int = 3`: Number of Krylov vectors used for GSE.
 - `extension_applyH_cutoff::Float64 = Float64_threshold()`: Cutoff for the application
   the MPO to the MPS.
 - `extension_applyH_maxdim::Int = maxlinkdim(psi) + 1`: Maximum bond/link
   dimension of the resulting MPS after the application of the MPO to the MPS.
 - `extension_cutoff::Float64 = 1E-10`: Cutoff for the basis extension step in GSE.
"""
krylov_extend!(engine::TDVPEngine{ProjMPO}; kwargs...) =
    krylov_extend!(engine.sysenv; kwargs...)

#################################################################################

"""
    function sweepdata(engine::TDVPEngine)

Returns the (shallow copy of) `SweepData`.
"""
sweepdata(engine::TDVPEngine) = engine.swdata

#################################################################################

"""
    function abstime(engine::TDVPEngine)

Returns the absolute elapsed time.
"""
abstime(engine::TDVPEngine) = engine.abstime

#################################################################################

"""
    function updateH!(engine::TDVPEngine, H::T;
                      recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel}

Update Hamiltonian `H` in `engine::TDVPEngine`. If `recalcEnv = false`,
it reuses previous environments. `recalcEnv = false` is only defined when the
`TDVPEngine` is created from a single `MPO`.
"""
updateH!(engine::TDVPEngine, H::T;
         recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
    updateH!(engine.sysenv, H; recalcEnv)

#################################################################################

"""
    function updateH!(engine::TDVPEngine, H::T, Ms::Vector{MPS};
                      weight::Float64,
                      recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel}

Update Hamiltonian `H` in `engine::TDVPEngine`. `recalcEnv = false` is not supported.
"""
updateH!(engine::TDVPEngine, H::T, Ms::Vector{MPS};
         weight::Float64,
         recalcEnv::Bool = true) where T <: Union{MPO, Vector{MPO}, CouplingModel} =
             updateH!(engine.sysenv, H, Ms; weight, recalcEnv)

#################################################################################

"""
    tdvpsweep!(engine::TDVPEngine, time_step::Union{Float64, ComplexF64}; 
               nsite::Union{Int, String} = "dynamic", 
               solver = exp_solver,
               kwargs...)::Nothing

Performs one TDVP sweep. Computes `ψ' = exp(time_step * H) * ψ`.
Therefore, for real-time dynamics with step `dt`, `time_step` should be `-im * dt`. 

#### Arguments:
 - `engine::TDVPEngine`.
 - `time_step::Union{Float64, ComplexF64}`.
 - `nsite::Union{Int, String} = "dynamic"`: For two or one site sweeps, `nsite=2` or
   `nsite=1` respectively. For `nsite="dynamic"`, [`dynamic_fullsweep!`](@ref) is performed.
 - `solver = exp_solver`: Only `exp_solver` is supported now.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after every
   fullsweep, if `2` prints at every update step.
 - `eigthreshold::Float64 = 1E-12`. Only applicable for `nsite = "dynamic"`
   (see [`dynamic_fullsweep!`](@ref)).
 - `extendat::Union{Nothing, Int} = nothing`: If specified, at every `extendat`th sweep,
   Global Subspace Expansion is performed followed by a pure one-site sweep if
   `typeof(sysenv) == StateEnvs{ProjMPO}`, else performs a full two-site sweep.
   Only applicable for `nsite = "dynamic"` (see [`dynamic_fullsweep!`](@ref)).

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-12`.
 - `solver_krylovdim::Int = 30`.
 - `solver_maxiter::Int = 100`.
 - `solver_outputlevel::Int = 0`.: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = true`.
 - `solver_check_convergence::Bool = true`.

#### Arguments for Global Subspace Expansion and their default values:
Only applicable for `nsite = "dynamic"` (see [`dynamic_fullsweep!`](@ref)).
 - `extension_krylovdim::Int = 3`: Number of Krylov vectors used for GSE.
 - `extension_applyH_cutoff::Float64 = Float64_threshold()`: Cutoff for the application
   the MPO to the MPS.
 - `extension_applyH_maxdim::Int = maxlinkdim(psi) + 1`: Maximum bond/link
   dimension of the resulting MPS after the application of the MPO to the MPS.
 - `extension_cutoff::Float64 = 1E-10`: Cutoff for the basis extension step in GSE.
"""
function tdvpsweep!(engine::TDVPEngine, time_step::Union{Float64, ComplexF64}; 
                    nsite::Union{Int, String} = "dynamic", 
                    solver = exp_solver,
                    kwargs...)
    
    if solver != exp_solver
        error("`solver` must be `exp_solver` for `tdvpsweep` !!")
    end
    
    if nsite == "dynamic"
        dynamic_fullsweep!(engine.sysenv, solver, engine.swdata;
                           time_step = 0.5 * time_step,
                           reverse_step=true,
                           kwargs...)
    elseif nsite == 2 || nsite == 1
        fullsweep!(engine.sysenv, solver, nsite, engine.swdata;
                   time_step = 0.5 * time_step,
                   reverse_step=true,
                   kwargs...)
    else
        error("""`nsite` must be `"dynamic"`, `2`, or `1` for `tdvpsweep` !!""")
    end
    engine.abstime += abs(time_step)
    return nothing
end

#################################################################################

