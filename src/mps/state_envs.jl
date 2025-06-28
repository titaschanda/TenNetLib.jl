
#################################################################################

"""
    mutable struct StateEnvs{T <: Union{ProjMPO,
                                  ProjMPO_MPS2,
                                  ProjMPOSum2,
                                  ProjMPOSum_MPS,
                                  ProjCouplingModel,
                                  ProjCouplingModel_MPS}
                            }
        psi::MPS
        PH::T
    end

Holds the MPS state `psi` and its environments `PH`.
"""
mutable struct StateEnvs{T <: Union{ProjMPO,
                                    ProjMPO_MPS2,
                                    ProjMPOSum2,
                                    ProjMPOSum_MPS,
                                    ProjCouplingModel,
                                    ProjCouplingModel_MPS}
                         }
    psi::MPS
    PH::T
end

#################################################################################

"""
    function getpsi(sysenv::StateEnvs)

Returns (shallow copy of) the state `psi`.
"""
getpsi(sysenv::StateEnvs) = Base.copy(sysenv.psi)

#################################################################################

"""
    function getenv(sysenv::StateEnvs)

Returns (shallow copy of) the environment `PH`.
"""
getenv(sysenv::StateEnvs) = Base.copy(sysenv.PH)

#################################################################################

"""
    function StateEnvs(psi::MPS, H::MPO)

Constructor of the `StateEnvs` from a single `MPO`.
"""
function StateEnvs(psi::MPS, H::MPO)
    ITensorMPS.check_hascommoninds(siteinds, H, psi)
    ITensorMPS.check_hascommoninds(siteinds, H, psi')
    H = permute(H, (linkind, siteinds, linkind))
    PH = ProjMPO(H)
    return StateEnvs(copy(psi), PH)
end

#################################################################################

"""
    function StateEnvs(psi::MPS, Hs::Vector{MPO})

Constructor of the  `StateEnvs` from a vector of `MPO`.
Environments are contracted in parallel.
"""
function StateEnvs(psi::MPS, Hs::Vector{MPO})
    for H in Hs
        ITensorMPS.check_hascommoninds(siteinds, H, psi)
        ITensorMPS.check_hascommoninds(siteinds, H, psi')
    end
    Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))
    PHS = ProjMPOSum2(Hs)
    return StateEnvs(copy(psi), PHS)
end

#################################################################################

"""
    function StateEnvs(psi::MPS, H::MPO, Ms::Vector{MPS}; weight::Float64)

Constructor of the `StateEnvs` from a single `MPO` and a vector of `MPS` used
for excited state DMRG.
"""
function StateEnvs(psi::MPS, H::MPO, Ms::Vector{MPS}; weight::Float64)
    
    ITensorMPS.check_hascommoninds(siteinds, H, psi)
    ITensorMPS.check_hascommoninds(siteinds, H, psi')
    for M in Ms
        ITensorMPS.check_hascommoninds(siteinds, M, psi)
    end
    H = permute(H, (linkind, siteinds, linkind))
    Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
    if weight <= 0.0
        error(string("`weight` parameter should be > 0.0",
                     " (value passed was `weight=$weight`)")
              )
    end
    PMM = ProjMPO_MPS2(H, Ms; weight=weight)
    return StateEnvs(copy(psi), PMM)
end

#################################################################################

"""
    function StateEnvs(psi::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}; weight::Float64)

Constructor of the `StateEnvs` from a vector of `MPO` and a vector of `MPS` used
for excited state DMRG. Environments for `MPO`s are contracted in parallel.
"""
function StateEnvs(psi::MPS, Hs::Vector{MPO}, Ms::Vector{MPS}; weight::Float64)

    for H in Hs
        ITensorMPS.check_hascommoninds(siteinds, H, psi)
        ITensorMPS.check_hascommoninds(siteinds, H, psi')
    end
    Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))

    for M in Ms
        ITensorMPS.check_hascommoninds(siteinds, M, psi)
    end
    H = permute(H, (linkind, siteinds, linkind))
    Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
    if weight <= 0.0
        error(string("`weight` parameter should be > 0.0",
                     " (value passed was `weight=$weight`)")
              )
    end
    PMM = ProjMPOSum_MPS(Hs, Ms; weight=weight)
    return StateEnvs(copy(psi), PMM)
end

#################################################################################

"""
    function StateEnvs(psi::MPS, H::CouplingModel)

Constructor of the `StateEnvs` from a `CouplingModel`. Each terms in the `CouplingModel`
are contracted in parallel.
"""
function StateEnvs(psi::MPS, H::CouplingModel)

    @assert siteinds(psi) == siteinds(H)
    PH = ProjCouplingModel(H)
    return StateEnvs(copy(psi), PH)
end

#################################################################################

"""
    function StateEnvs(psi::MPS, H::CouplingModel, Ms::Vector{MPS}; weight::Float64)

Constructor of StateEnvs from a `CouplingModel` and a vector of MPS used
for excited state DMRG. Each terms in the `CouplingModel` are contracted in parallel.
"""
function StateEnvs(psi::MPS, H::CouplingModel, Ms::Vector{MPS}; weight::Float64)
    
    @assert siteinds(psi) == siteinds(H)

    Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
    if weight <= 0.0
        error(string("`weight` parameter should be > 0.0",
                     " (value passed was `weight=$weight`)")
              )
    end
    PMM = ProjCouplingModel_MPS(H, Ms; weight=weight)
    return StateEnvs(copy(psi), PMM)
end

#################################################################################

"""
    function updateH!(sysenv::StateEnvs{ProjMPO}, H::MPO; recalcEnv::Bool = true)

Update Hamiltonian `H` in `sysenv::StateEnvs`. If `recalcEnv = false`,
it reuses previous environments.
"""
function updateH!(sysenv::StateEnvs{ProjMPO}, H::MPO; recalcEnv::Bool = true)
    
    if recalcEnv
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi)
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi')
        H = permute(H, (linkind, siteinds, linkind))
        PH = ProjMPO(H)
        sysenv.PH = PH
    else
        @assert sysenv.PH.lpos == 0
        @assert typeof(sysenv.PH) == ProjMPO

        if !isortho(sysenv.psi) || orthocenter(sysenv.psi) != 1
            orthogonalize!(sysenv.psi, 1)
        end
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi)
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi')
        H = permute(H, (linkind, siteinds, linkind))
        
        for ii = 1:length(H)    
            linkindNew = filter(x -> hastags(x, "Link"), inds(sysenv.PH.H[ii]))
            linkindOld = filter(x -> hastags(x, "Link"), inds(H[ii]))
            replaceinds!(H[ii], linkindOld, linkindNew)
        end            
        
        sysenv.PH.H = H
    end
end

#################################################################################

"""
    function updateH!(sysenv::StateEnvs{ProjMPOSum2}, Hs::Vector{MPO}; recalcEnv::Bool = true)

Update Hamiltonian `H` in `sysenv::StateEnvs`. `recalcEnv = false` is not supported.
"""
function updateH!(sysenv::StateEnvs{ProjMPOSum2}, Hs::Vector{MPO}; recalcEnv::Bool = true)
    
    if recalcEnv
        for H in Hs
            ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi)
            ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi')
        end
        Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))
        PHS = ProjMPOSum2(Hs)
        sysenv.PH = PHS
    else
        error("`updateH!()` :: Not implemented for `recalcEnv=$recalcEnv` with `ProjMPOSum2` !!")
    end
end

#################################################################################

"""
    function updateH!(sysenv::StateEnvs{ProjMPO_MPS2}, H::MPO, Ms::Vector{MPS};
                      weight::Float64, recalcEnv::Bool = true)

Update Hamiltonian `H` in `sysenv::StateEnvs`. `recalcEnv = false` is not supported.
 - `weight::Union{Nothing, Float64} = nothing`.
"""
function updateH!(sysenv::StateEnvs{ProjMPO_MPS2}, H::MPO, Ms::Vector{MPS};
                  weight::Float64, recalcEnv::Bool = true)
    
    if recalcEnv
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi)
        ITensorMPS.check_hascommoninds(siteinds, H, sysenv.psi')
        for M in Ms
            ITensorMPS.check_hascommoninds(siteinds, M, psi)
        end
        H = permute(H, (linkind, siteinds, linkind))
        Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
        if weight <= 0.0
            error(string("`weight` parameter should be > 0.0",
                         " (value passed was `weight=$weight`)")
                  )
        end
        PMM = ProjMPO_MPS2(H, Ms; weight=weight)
        sysenv.PH = PMM
    else
        error(string("`updateH!()` :: Not implemented for `recalcEnv=$recalcEnv` ",
                     "with `ProjMPO_MPS2` !!"))
    end
end

#################################################################################

"""
    function updateH!(sysenv::StateEnvs{ProjMPOSum_MPS}, Hs::Vector{MPO}, Ms::Vector{MPS};
                      weight::Float64, recalcEnv::Bool = true)

Update Hamiltonian `H` in `sysenv::StateEnvs`. `recalcEnv = false` is not supported.
 - `weight::Union{Nothing, Float64} = nothing`.
"""
function updateH!(sysenv::StateEnvs{ProjMPOSum_MPS}, H::Vector{MPO}, Ms::Vector{MPS};
                  weight::Float64, recalcEnv::Bool = true)
    
    if recalcEnv

        for H in Hs
            ITensorMPS.check_hascommoninds(siteinds, H, psi)
            ITensorMPS.check_hascommoninds(siteinds, H, psi')
        end
        Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))
                
        for M in Ms
            ITensorMPS.check_hascommoninds(siteinds, M, psi)
        end
        H = permute(H, (linkind, siteinds, linkind))
        Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
        if weight <= 0.0
            error(string("`weight` parameter should be > 0.0",
                         " (value passed was `weight=$weight`)")
                  )
        end
        PMM = ProjMPOSum_MPS(Hs, Ms; weight=weight)
        sysenv.PH = PMM
    else
        error(string("`updateH!()` :: Not implemented for `recalcEnv=$recalcEnv` ",
                     "with `ProjMPOSum_MPS` !!"))
    end
end

#################################################################################

"""
    function updateH!(sysenv::StateEnvs{ProjCouplingModel_MPS}, H::CouplingModel, Ms::Vector{MPS};
                      weight::Float64, recalcEnv::Bool = true)

Update Hamiltonian `H` in `sysenv::StateEnvs`. `recalcEnv = false` is not supported.
- `weight::Union{Nothing, Float64} = nothing`.
"""
function updateH!(sysenv::StateEnvs{ProjCouplingModel_MPS}, H::CouplingModel, Ms::Vector{MPS};
                  weight::Float64, recalcEnv::Bool = true)
    
    if recalcEnv
        @assert siteinds(psi) == siteinds(H)

        Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
        if weight <= 0.0
            error(string("`weight` parameter should be > 0.0",
                         " (value passed was `weight=$weight`)")
                  )
        end
        PMM = ProjCouplingModel_MPS(H, Ms; weight=weight)
        sysenv.PH = PMM
    else
        error(string("`updateH!()` :: Not implemented for `recalcEnv=$recalcEnv` ",
                     "with `ProjCouplingModel_MPS` !!"))
    end
end

#################################################################################

"""
    function nsite(sysenv::StateEnvs)

Returns the `nsite` of the environment. `nsite = 1` for single-site environment,
`nsite = 2` for two-site environmrnt, and so on.
 Currently, only uses `nsite = 0`, `1`, or `2`.
"""
nsite(sysenv::StateEnvs) = nsite(sysenv.PH)

#################################################################################

"""
    function set_nsite!(sysenv::StateEnvs, nsite::Int)

Set `nsite` of the environment. `nsite = 1` for single-site environment,
`nsite = 2` for two-site environmrnt, and so on.
 Currently, only uses `nsite = 0`, `1`, or `2`.
""" 
function set_nsite!(sysenv::StateEnvs, nsite::Int)
    set_nsite!(sysenv.PH, nsite)
    return sysenv
end

#################################################################################

"""
    function position!(sysenv::StateEnvs, pos::Int)

Compute the left and right environments at position `pos`.
"""
function position!(sysenv::StateEnvs, pos::Int)
    position!(sysenv.PH, sysenv.psi, pos)
    return sysenv
end

#################################################################################

"""
    function product(sysenv::StateEnvs, v::ITensor)

Returns the `Matrix-Vector` product between the environment and input ITensor `v`.
"""
product(sysenv::StateEnvs, v::ITensor)::ITensor = sysenv.PH(v)

(sysenv::StateEnvs)(v::ITensor) = product(sysenv, v)

#################################################################################

"""
    function Base.copy(sysenv::StateEnvs)

Shallow copy of `StateEnvs`.
"""
Base.copy(sysenv::StateEnvs) = StateEnvs(Base.copy(sysenv.psi), Base.copy(sysenv.PH))

#################################################################################

"""
    function Base.length(sysenv::StateEnvs)

Returns the length of the underlying MPS/Environment.
"""
Base.length(sysenv::StateEnvs) = Base.length(sysenv.psi)

#################################################################################


