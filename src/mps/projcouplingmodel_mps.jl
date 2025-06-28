
#################################################################################

"""
    mutable struct ProjCouplingModel_MPS
        PH::ProjCouplingModel
        pm::Vector{ProjMPS2}
        weight::Float64
    end

Holds Environments from CouplingModel and MPS.
"""
mutable struct ProjCouplingModel_MPS
    PH::ProjCouplingModel
    pm::Vector{ProjMPS2}
    weight::Float64
end

#################################################################################

"""
    ProjCouplingModel_MPS(H::CouplingModel, mpsv::Vector{MPS}; weight=1.0)

Constructor of ProjCouplingModel_MPS.
"""
function ProjCouplingModel_MPS(H::CouplingModel, mpsv::Vector{MPS}; weight=1.0)
    return ProjCouplingModel_MPS(ProjCouplingModel(H), [ProjMPS2(m) for m in mpsv], weight)
end

#################################################################################

"""
    Base.length(P::ProjCouplingModel_MPS)

Returns the length of the underlying Environment.
"""
Base.length(P::ProjCouplingModel_MPS) = Base.length(P.PH)

#################################################################################

"""
    Base.copy(P::ProjCouplingModel_MPS)

Shallow copy ProjCouplingModel_MPS.
"""
Base.copy(P::ProjCouplingModel_MPS) =
    ProjCouplingModel_MPS(Base.copy(P.PH), Base.copy.(P.pm), P.weight)

#################################################################################

"""
    nsite(P::ProjCouplingModel_MPS)

Returns `nsite` of ProjCouplingModel_MPS
"""
nsite(P::ProjCouplingModel_MPS) = nsite(P.PH)

#################################################################################

"""
    set_nsite!(P::ProjCouplingModel_MPS, nsite::Int)::ProjCouplingModel_MPS

Set `nsite` of the `ProjCouplingModel_MPS` object.
""" 
function set_nsite!(Ps::ProjCouplingModel_MPS, nsite::Int)::ProjCouplingModel_MPS
    set_nsite!(Ps.PH, nsite)
    for P in Ps.pm
        set_nsite!(P, nsite)
    end
    return Ps
end

#################################################################################

"""
    Base.eltype(Ps::ProjCouplingModel_MPS)

Returns the element type (e.g., `Float64` or `ComplexF64`) of `ProjCouplingModel_MPS`.
"""
function Base.eltype(Ps::ProjCouplingModel_MPS)
    elt = eltype(Ps.PH)
    for p in Ps.pm
        elt = promote_type(elt, eltype(p))
    end
    return elt
end

#################################################################################

"""
    product(P::ProjCouplingModel_MPS, v::ITensor)

Compute the `Matrix-Vector` product between ProjCouplingModel_MPS and input ITensor `v`.
"""
function product(P::ProjCouplingModel_MPS, v::ITensor)::ITensor
    Pv = P.PH(v)
    for p in P.pm
        Pv += P.weight * product(p, v)
    end
    return Pv
end

"""
    (P::ProjCouplingModel_MPS)(v::ITensor)

Compute the `Matrix-Vector` product between ProjCouplingModel_MPS and input ITensor `v`.
"""
(P::ProjCouplingModel_MPS)(v::ITensor) = product(P, v)

#################################################################################

"""
    position!(P::ProjCouplingModel_MPS, psi::MPS, k::Int)::ProjCouplingModel_MPS

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjCouplingModel_MPS, psi::MPS, pos::Int)::ProjCouplingModel_MPS
    position!(P.PH, psi, pos)
    for p in P.pm
        position!(p, psi, pos)
    end
    return P
end

#################################################################################

ITensorMPS.noiseterm(P::ProjCouplingModel_MPS, phi::ITensor, ortho::String) =
    noiseterm(P.PH, phi, ortho)

#################################################################################

