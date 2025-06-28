
#################################################################################

"""
    mutable struct StateEnvsTTN{T1, T2 <: Union{EnvCouplingModelTTN,
                                                EnvCouplingModelProjTTN}}
        psi::TTN{T1}
        env::T{T2}
    end

Holds the TTN and the environment.
"""
mutable struct StateEnvsTTN{T1, T2 <: Union{EnvCouplingModelTTN,
                                            EnvCouplingModelProjTTN}}
    psi::TTN{T1}
    env::T2
end

#################################################################################

"""
    function StateEnvsTTN(psi::TTN, H::CouplingModel)

Constructor of the `StateEnvsTTN` from a `CouplingModel`. Each terms in the `CouplingModel`
are contracted in parallel.
"""
function StateEnvsTTN(psi::TTN, M::CouplingModel)
    psi = copy(psi)
    env = EnvCouplingModelTTN(psi, M)
    return StateEnvsTTN(psi, env)
end

#################################################################################

"""
    function StateEnvsTTN(psi::TTN, H::CouplingModel, Ms::Vector{TTN}; weight::Float64)

Constructor of StateEnvsTTN from a `CouplingModel` and a vector of TTN used
for excited state search. Each terms in the `CouplingModel` are contracted in parallel.
"""
function StateEnvsTTN(psi::TTN,
                      M::CouplingModel,
                      psis::Vector{TTN};
                      weight::Float64)
    psi = copy(psi)
    if weight <= 0.0
        error(string("`weight` parameter should be > 0.0",
                     " (value passed was `weight=$weight`)")
              )
    end
    env = EnvCouplingModelProjTTN(psi, M, psis; weight = weight)
    return StateEnvsTTN(psi, env)
end

function StateEnvsTTN(psi::TTN{T},
                      M::CouplingModel,
                      psis::Vector{TTN{T}};
                      weight::Float64) where T
    psi = copy(psi)
    if weight <= 0.0
        error(string("`weight` parameter should be > 0.0",
                     " (value passed was `weight=$weight`)")
              )
    end
    env = EnvCouplingModelProjTTN(psi, M, psis; weight = weight)
    return StateEnvsTTN(psi, env)
end

StateEnvsTTN(psi::TTN, M::CouplingModel, psis::TTN...; weight::Float64) =
    StateEnvsTTN(psi, M, TTN[psis...]; weight = weight)

#################################################################################

"""
    function Base.copy(sysenv::StateEnvsTTN)

Shallow copy of `StateEnvsTTN`.
"""
Base.copy(sysenv::StateEnvsTTN) =
    StateEnvsTTN(Base.copy(sysenv.psi), Base.copy(sysenv.env))

#################################################################################

"""
    function getpsi(sysenv::StateEnvsTTN)

Returns (shallow copy of) the state `psi`.
"""
getpsi(sysenv::StateEnvsTTN) = Base.copy(sysenv.psi)

#################################################################################

"""
    function getenv(sysenv::StateEnvsTTN)

Returns (shallow copy of) the environment `env`.
"""
getenv(sysenv::StateEnvsTTN) = base.copy(sysenv.env)

#################################################################################

"""
    function position!(sysenv::StateEnvsTTN, node::Int2;
                       maxdim::Int = typemax(Int),
                       mindim::Int = 1,
                       cutoff::Float64 = Float64_threshold(),
                       svd_alg::String = "divide_and_conquer",
                       normalize::Bool = true)

Moves the orthogonality center of the TTN to `node` and updates the effective Hamiltonian of
the same.

#### Named arguments and their default values:
 - `normalize::Bool = true`: Whether to normalize afterwards.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
"""
function position!(sysenv::StateEnvsTTN, node::Int2;
                   maxdim::Int = typemax(Int),
                   mindim::Int = 1,
                   cutoff::Float64 = Float64_threshold(),
                   svd_alg::String = "divide_and_conquer",
                   normalize::Bool = true,
                   node_to_skip::Union{Int2, Nothing} = nothing
                   )::Nothing
    oc = orthocenter(sysenv.psi)

    isometrize!(sysenv.psi, node;
                normalize=normalize,
                maxdim=maxdim,
                mindim=mindim,
                cutoff=cutoff,
                svd_alg=svd_alg)
    
    move_environment!(sysenv.env, sysenv.psi, oc, node;
                      node_to_skip=node_to_skip)
    return nothing
end

#################################################################################

Base.eltype(sysenv::StateEnvsTTN) = Base.eltype(sysenv.env, sysenv.psi)

#################################################################################

"""
    function product(sysenv::StateEnvsTTN, v::ITensor)

Returns the `Matrix-Vector` product between the environment and input ITensor `v`.
"""
product(sysenv::StateEnvsTTN, v::ITensor) = product(sysenv.env, sysenv.psi, v)

(sysenv::StateEnvsTTN)(v::ITensor) = product(sysenv, v)

#################################################################################

