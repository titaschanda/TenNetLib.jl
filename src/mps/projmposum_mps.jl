
#################################################################################

"""
    mutable struct ProjMPOSum_MPS
        PH::ProjMPOSum2
        pm::Vector{ProjMPS2}
        weight::Float64
    end

Holds Environments from MPOs and MPSs.
"""
mutable struct ProjMPOSum_MPS
    PH::ProjMPOSum2
    pm::Vector{ProjMPS2}
    weight::Float64
end

#################################################################################

"""
    ProjMPOSum_MPS(Hs::Vector{MPO}, mpsv::Vector{MPS}; weight=1.0)

Constructor of the `ProjMPOSum_MPS`.
"""
function ProjMPOSum_MPS(Hs::Vector{MPO}, mpsv::Vector{MPS}; weight=1.0)
    return ProjMPOSum_MPS(ProjMPOSum2(Hs), ProjMPS2[ProjMPS2(m) for m in mpsv], weight)
end

#################################################################################

"""
    Base.length(P::ProjMPOSum_MPS)

Returns the length of the underlying Environment.
"""
Base.length(P::ProjMPOSum_MPS) = Base.length(P.PH)

#################################################################################

"""
    Base.copy(P::ProjMPOSum_MPS)

Shallow copy `ProjMPOSum_MPS`.
"""
Base.copy(P::ProjMPOSum_MPS) = ProjMPOSum_MPS(Base.copy(P.PH), Base.copy.(P.pm), P.weight)

#################################################################################

"""
    nsite(P::ProjMPOSum_MPS)

Returns `nsite` of ProjMPOSum_MPS
"""
nsite(P::ProjMPOSum_MPS) = nsite(P.PH)

#################################################################################

"""
    set_nsite!(P::ProjMPOSum_MPS, nsite::Int)::ProjMPOSum_MPS

Set `nsite` of the `ProjMPOSum_MPS` object.
""" 
function set_nsite!(Ps::ProjMPOSum_MPS, nsite::Int)::ProjMPOSum_MPS
    set_nsite!(Ps.PH, nsite)
    for P in Ps.pm
        set_nsite!(P, nsite)
    end
    return Ps
end

#################################################################################

"""
    Base.eltype(Ps::ProjMPOSum_MPS)

Returns the element type (e.g., `Float64` or `ComplexF64`) of `ProjMPOSum_MPS`.
"""
function Base.eltype(Ps::ProjMPOSum_MPS)
  elt = eltype(Ps.PH)
  for p in Ps.pm
    elt = promote_type(elt, eltype(p))
  end
  return elt
end

#################################################################################

"""
    product(P::ProjMPOSum_MPS, v::ITensor)::ITensor

Compute the `Matrix-Vector` product between ProjMPOSum_MPS and input ITensor `v`.
"""
function product(P::ProjMPOSum_MPS, v::ITensor)::ITensor
    Pv = P.PH(v)
    for p in P.pm
        Pv += P.weight * product(p, v)
    end
    return Pv
end

"""
    (P::ProjMPOSum_MPS)(v::ITensor)

Compute the `Matrix-Vector` product between ProjMPOSum_MPS and input ITensor `v`.
"""
(P::ProjMPOSum_MPS)(v::ITensor) = product(P, v)

#################################################################################

Base.size(P::ProjMPOSum_MPS) = Base.size(P.PH)

#################################################################################

"""
    position!(P::ProjMPOSum_MPS, psi::MPS, k::Int)::ProjMPOSum_MPS

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjMPOSum_MPS, psi::MPS, pos::Int)::ProjMPOSum_MPS
    position!(P.PH, psi, pos)
    for p in P.pm
        position!(p, psi, pos)
    end
    return P
end

#################################################################################

ITensorMPS.noiseterm(P::ProjMPOSum_MPS, phi::ITensor, ortho::String) =
    noiseterm(P.PH, phi, ortho)

#################################################################################

