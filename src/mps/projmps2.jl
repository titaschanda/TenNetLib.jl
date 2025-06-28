    
#################################################################################

"""
    mutable struct ProjMPS2 <: AbstractProjMPO
        lpos::Int
        rpos::Int
        nsite::Int
        M::MPS
        LR::Vector{ITensor}
    end


Holds the following data where `psi`
is the MPS being optimized and `M` is the 
MPS held constant by the ProjMPS.

Inherited from `AbstractProjMPO` from `ITesnors.jl`.
```
     o--o--o--o--o--o--o--o--o--o--o <M|
     |  |  |  |  |  |  |  |  |  |  |
     *--*--*-      -*--*--*--*--*--* |psi>
```
Required for calculating projectors for a MPS.
"""
mutable struct ProjMPS2 <: AbstractProjMPO
    lpos::Int
    rpos::Int
    nsite::Int
    M::MPS
    LR::Vector{ITensor}
end

#################################################################################

"""
    ProjMPS2(M::MPS)::ProjMPS2

Constructor of `ProjMPS2`.
"""
function ProjMPS2(M::MPS)::ProjMPS2
    return ProjMPS2(0, length(M) + 1, 2, M, Vector{ITensor}(undef, length(M)))
end

#################################################################################

"""
    Base.length(P::ProjMPS2)

Returns the length of the underlying MPS/Environment.
"""
Base.length(P::ProjMPS2) = Base.length(P.M)

#################################################################################

"""
    Base.copy(P::ProjMPS2)

Shallow copy `ProjMPS2`.
"""
Base.copy(P::ProjMPS2) = ProjMPS2(P.lpos, P.rpos, P.nsite, Base.copy(P.M), Base.copy(P.LR))

#################################################################################

"""
    set_nsite!(P::ProjMPS2, nsite::Int)::ProjMPS2

Set `nsite` of the `ProjMPS2` object.
""" 
function set_nsite!(P::ProjMPS2, nsite::Int)::ProjMPS2
    P.nsite = nsite
    return P
end

#################################################################################

function _makeL!(P::ProjMPS2, psi::MPS, k::Int)::Union{ITensor, Nothing}
    ll = P.lpos
    if ll ≥ k
        P.lpos = k
        return nothing
    end
    ll = max(ll, 0)
    L = lproj(P)
    while ll < k
        L = L * psi[ll + 1] * dag(prime(P.M[ll + 1], "Link"))
        P.LR[ll + 1] = L
        ll += 1
    end
    P.lpos = k
    return L
end

#################################################################################

"""
    makeL!(P::ProjMPS2, psi::MPS, k::Int)::ProjMPS2

Compute Left Environments for the MPS `psi` at position `k`.
"""
function makeL!(P::ProjMPS2, psi::MPS, k::Int)::ProjMPS2
    _makeL!(P, psi, k)
    return P
end

#################################################################################
    
function _makeR!(P::ProjMPS2, psi::MPS, k::Int)::Union{ITensor, Nothing}
    rl = P.rpos
    if rl ≤ k
        P.rpos = k
        return nothing
    end
    N = length(P.M)
    rl = min(rl, N + 1)
    R = rproj(P)
    while rl > k
        R = R * psi[rl - 1] * dag(prime(P.M[rl - 1], "Link"))
        P.LR[rl - 1] = R
        rl -= 1
    end
    P.rpos = k
    return R
end

#################################################################################

"""
    makeR!(P::ProjMPS2, psi::MPS, k::Int)::ProjMPS2

Compute Right Environments for the MPS `psi` at position `k`.
"""
function makeR!(P::ProjMPS2, psi::MPS, k::Int)::ProjMPS2
    _makeR!(P, psi, k)
    return P
end

#################################################################################

"""
    position!(P::ProjMPS2, psi::MPS, k::Int)::ProjMPS2

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjMPS2, psi::MPS, pos::Int)::ProjMPS2
    makeL!(P, psi, pos - 1)
    makeR!(P, psi, pos + nsite(P))
    return P
end

#################################################################################

"""
    contract(P::ProjMPS2, v::ITensor)::ITensor

Tensor contraction with input ITensor `v` with the ProjMPS2.
Note: This is different from `product`.
"""
function contract(P::ProjMPS2, v::ITensor)::ITensor
    itensor_map = Union{ITensor,OneITensor}[lproj(P)]
    append!(itensor_map, [dag(prime(t, "Link")) for t in P.M[site_range(P)]])
    push!(itensor_map, rproj(P))

    if dim(first(itensor_map)) == 1
        reverse!(itensor_map)
    end

    Mv = v
    for it in itensor_map
        Mv *= it
    end
    return Mv
end

#################################################################################

"""
    proj_mps(P::ProjMPS2)::ITensor

Compute the projector of `ProjMPS2`.
"""
function proj_mps(P::ProjMPS2)::ITensor
    itensor_map = Union{ITensor,OneITensor}[lproj(P)]
    append!(itensor_map, [dag(prime(t, "Link")) for t in P.M[site_range(P)]])
    push!(itensor_map, rproj(P))

    if dim(first(itensor_map)) == 1
        reverse!(itensor_map)
    end

    m = ITensor(true)
    for it in itensor_map
        m *= it
    end
    return m
end

#################################################################################

"""
    product(P::ProjMPS2, v::ITensor)::ITensor

Compute the `Matrix-Vector` product between ProjMPS2 and input ITensor `v`.
"""
function product(P::ProjMPS2, v::ITensor)::ITensor

    Pv = contract(P, v) * dag(proj_mps(P))
    if order(Pv) != order(v)
        error(
            string(
                "The order of the ProjMPS2-ITensor product P*v is not equal to the ",
                "order of the ITensor v, this is probably due to an index mismatch.\n",
                "Common reasons for this error: \n",
                "(1) You are trying to multiply the ProjMPS2 with the $(nsite(P))-site ",
                "wave-function at the wrong position.\n",
                "(2) `orthogonalize!` was called, changing the MPS without updating ",
                "the ProjMPS.\n\n",
                "P*v inds: $(inds(Pv)) \n\n",
                "v inds: $(inds(v))",
            ),
        )
    end
    return Pv 
end

#################################################################################

"""
    (P::ProjMPS2)(v::ITensor)

Compute the `Matrix-Vector` product between ProjMPS2 and input ITensor `v`.
"""
(P::ProjMPS2)(v::ITensor) = product(P, v)

#################################################################################
