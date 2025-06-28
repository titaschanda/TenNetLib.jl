
#################################################################################

"""
    mutable struct EnvCouplingModelTTN
        env_tensors::LinkTensorsTTN
    end

A wrapper around `LinkTensorsTTN`.
"""
mutable struct EnvCouplingModelTTN
    env_tensors::LinkTensorsTTN
end

#################################################################################

EnvCouplingModelTTN(psi::TTN, M::CouplingModel) =
    EnvCouplingModelTTN(LinkTensorsTTN(psi::TTN, M::CouplingModel))

#################################################################################

move_environment!(env::EnvCouplingModelTTN, psi::TTN,
                  source_node::Int2,
                  destination_node::Int2;
                  kwargs...) = move_linktensors!(env.env_tensors, psi,
                                                 source_node,
                                                 destination_node;
                                                 kwargs...)

#################################################################################

Base.eltype(env::EnvCouplingModelTTN, psi::TTN) = Base.eltype(env.env_tensors, psi)

#################################################################################

product(env::EnvCouplingModelTTN, psi::TTN, v::ITensor) =
    product(env.env_tensors, psi, v)

#################################################################################

"""
    mutable struct EnvCouplingModelProjTTN
        env_tensors::LinkTensorsTTN
        proj_tensors::Vector{LinkProjTTN}
        weight::Float64
    end

A wrapper to hold `LinkTensorsTTN` and a vector of `LinkProjTTN`.
"""
mutable struct EnvCouplingModelProjTTN
    env_tensors::LinkTensorsTTN
    proj_tensors::Vector{LinkProjTTN}
    weight::Float64
end

#################################################################################

EnvCouplingModelProjTTN(psi::TTN, M::CouplingModel, psis::Vector{TTN}; weight = 1.0) =
    EnvCouplingModelProjTTN(LinkTensorsTTN(psi, M),
                            LinkProjTTN[LinkProjTTN(psi, x) for x in psis],
                            weight)

EnvCouplingModelProjTTN(psi::TTN{T}, M::CouplingModel,
                        psis::Vector{TTN{T}}; weight = 1.0) where T =
    EnvCouplingModelProjTTN(LinkTensorsTTN(psi, M),
                            LinkProjTTN[LinkProjTTN(psi, x) for x in psis],
                            weight)

EnvCouplingModelProjTTN(psi::TTN, M::CouplingModel, psis::TTN...; weight = 1.0) =
    EnvCouplingModelProjTTN(psi, M, TTN[psis...]; weight = weight)

#################################################################################

function move_environment!(env::EnvCouplingModelProjTTN, psi::TTN,
                           source_node::Int2,
                           destination_node::Int2;
                           kwargs...)::Nothing
    
    move_linktensors!(env.env_tensors, psi,
                      source_node,
                      destination_node;
                      kwargs...)

    for proj in env.proj_tensors
        move_linkproj!(proj, psi,
                       source_node,
                       destination_node;
                       kwargs...)
    end
end

#################################################################################

function Base.eltype(env::EnvCouplingModelProjTTN, psi::TTN)
    elt = eltype(env.env_tensors, psi)
    for proj in env.proj_tensors
        elt = promote_type(elt, eltype(proj, psi))
    end
    return elt
end

#################################################################################

function product(env::EnvCouplingModelProjTTN, psi::TTN, v::ITensor)::ITensor
    Pv = product(env.env_tensors, psi, v)
    for proj in env.proj_tensors
        Pv += env.weight * product(proj, psi, v)
    end
    return Pv
end

#################################################################################

