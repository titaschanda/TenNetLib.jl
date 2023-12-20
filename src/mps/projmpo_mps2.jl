
#################################################################################

"""
    mutable struct ProjMPO_MPS2
        PH::ProjMPO
        pm::Vector{ProjMPS2}
        weight::Float64
    end

Holds Environments from MPO and MPS.
"""
mutable struct ProjMPO_MPS2
    PH::ProjMPO
    pm::Vector{ProjMPS2}
    weight::Float64
end

#################################################################################

"""
    ProjMPO_MPS2(H::MPO, mpsv::Vector{MPS}; weight=1.0)

Constructor of the `ProjMPO_MPS2`.
"""
function ProjMPO_MPS2(H::MPO, mpsv::Vector{MPS}; weight=1.0)
    return ProjMPO_MPS2(ProjMPO(H), ProjMPS2[ProjMPS2(m) for m in mpsv], weight)
end

#################################################################################

"""
    Base.length(P::ProjMPO_MPS2)

Returns the length of the underlying Environment.
"""
Base.length(P::ProjMPO_MPS2) = Base.length(P.PH)

#################################################################################

"""
    Base.copy(P::ProjMPO_MPS2)

Shallow copy `ProjMPO_MPS2`.
"""
Base.copy(P::ProjMPO_MPS2) = ProjMPO_MPS2(Base.copy(P.PH), Base.copy.(P.pm), P.weight)

#################################################################################

"""
    nsite(P::ProjMPO_MPS2)

Returns `nsite` of ProjMPO_MPS2
"""
nsite(P::ProjMPO_MPS2) = nsite(P.PH)

#################################################################################

"""
    set_nsite!(P::ProjMPO_MPS2, nsite::Int)::ProjMPO_MPS2

Set `nsite` of the `ProjMPO_MPS2` object.
""" 
function set_nsite!(Ps::ProjMPO_MPS2, nsite::Int)::ProjMPO_MPS2
    set_nsite!(Ps.PH, nsite)
    for P in Ps.pm
        set_nsite!(P, nsite)
    end
    return Ps
end

#################################################################################

"""
    Base.eltype(Ps::ProjMPO_MPS2)

Returns the element type (e.g., `Float64` or `ComplexF64`) of `ProjMPO_MPS2`.
"""
function Base.eltype(Ps::ProjMPO_MPS2)
  elt = eltype(Ps.PH)
  for p in Ps.pm
    elt = promote_type(elt, eltype(p))
  end
  return elt
end

#################################################################################

"""
    product(P::ProjMPO_MPS2, v::ITensor)::ITensor

Compute the `Matrix-Vector` product between ProjMPO_MPS2 and input ITensor `v`.
"""
function product(P::ProjMPO_MPS2, v::ITensor)::ITensor
    Pv = P.PH(v)
    for p in P.pm
        Pv += P.weight * product(p, v)
    end
    return Pv
end

"""
    (P::ProjMPO_MPS2)(v::ITensor)

Compute the `Matrix-Vector` product between ProjMPO_MPS2 and input ITensor `v`.
"""
(P::ProjMPO_MPS2)(v::ITensor) = product(P, v)

#################################################################################

Base.size(P::ProjMPO_MPS2) = Base.size(P.PH)

#################################################################################

"""
    position!(P::ProjMPO_MPS2, psi::MPS, k::Int)::ProjMPO_MPS2

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjMPO_MPS2, psi::MPS, pos::Int)::ProjMPO_MPS2
    position!(P.PH, psi, pos)
    for p in P.pm
        position!(p, psi, pos)
    end
    return P
end

#################################################################################

ITensors.noiseterm(P::ProjMPO_MPS2, phi::ITensor, ortho::String) =
    noiseterm(P.PH, phi, ortho)

#################################################################################
