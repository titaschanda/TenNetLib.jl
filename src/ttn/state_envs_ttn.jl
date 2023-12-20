
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

## TO ADD MORE

function StateEnvsTTN(psi::TTN, M::CouplingModel)
    psi = copy(psi)
    env = EnvCouplingModelTTN(psi, M)
    return StateEnvsTTN(psi, env)
end

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

Base.copy(sysenv::StateEnvsTTN) =
    StateEnvsTTN(Base.copy(sysenv.psi), Base.copy(sysenv.env))

#################################################################################

getpsi(sysenv::StateEnvsTTN) = Base.copy(sysenv.psi)

#################################################################################

function position!(sysenv::StateEnvsTTN, node::Int2;
                   maxdim::Int = typemax(Int),
                   mindim::Int = 1,
                   cutoff::Float64 = _Float64_Threshold,
                   svd_alg::String = "divide_and_conquer",
                   normalize::Bool = false,
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

product(sysenv::StateEnvsTTN, v::ITensor) = product(sysenv.env, sysenv.psi, v)
   
#################################################################################

(sysenv::StateEnvsTTN)(v::ITensor) = product(sysenv, v)

#################################################################################
