
#################################################################################

"""
    mutable struct ProjMPOSum2
        PHs::Vector{ProjMPO}
    end

Holds Environments from a sum of MPOs. 
"""
mutable struct ProjMPOSum2
    PHs::Vector{ProjMPO}
end

#################################################################################

"""
    ProjMPOSum2(Hs::Vector{MPO})

Constructor of the `ProjMPOSum2`.
"""
function ProjMPOSum2(Hs::Vector{MPO})
    return ProjMPOSum2(ProjMPO[ProjMPO(ph) for ph in Hs])
end

#################################################################################

"""
    Base.length(P::ProjMPOSum2)

Returns the length of the underlying Environment.
"""
Base.length(P::ProjMPOSum2) = Base.length(P.PHs[begin])

#################################################################################

"""
    Base.copy(P::ProjMPOSum2)

Shallow copy `ProjMPOSum2`.
"""
Base.copy(P::ProjMPOSum2) = ProjMPOSum2(Base.copy.(P.PHs))

#################################################################################

"""
    nsite(P::ProjMPOSum2)

Returns `nsite` of ProjMPOSum2
"""
nsite(P::ProjMPOSum2) = nsite(P.PHs[begin])

#################################################################################

"""
    set_nsite!(P::ProjMPOSum2, nsite::Int)::ProjMPOSum2

Set `nsite` of the `ProjMPOSum2` object.
""" 
function set_nsite!(Ps::ProjMPOSum2, nsite::Int)::ProjMPOSum2
    for P in Ps.PHs
        set_nsite!(P, nsite)
    end
    return Ps
end

#################################################################################

"""
    Base.eltype(Ps::ProjMPOSum2)

Returns the element type (e.g., `Float64` or `ComplexF64`) of `ProjMPOSum2`.
"""
function Base.eltype(Ps::ProjMPOSum2)
    elts = eltype.(Ps.PHs)
    return promote_type(elts...)
end

#################################################################################

"""
    product(P::ProjMPOSum2, v::ITensor)::ITensor

Compute the `Matrix-Vector` product between ProjMPOSum2 and input ITensor `v`.
"""
function product(P::ProjMPOSum2, v::ITensor)::ITensor

    Pv = ITensor()
    using_threaded_loop() && (mutex = Threads.SpinLock())
    @threaded_loop for p in P.PHs
        temp_pv = p(v)
        if using_threaded_loop()
            lock(mutex) do
                Pv += temp_pv
            end
        else
            Pv += temp_pv
        end
    end
    return Pv
end

"""
    (P::ProjMPOSum2)(v::ITensor)

Compute the `Matrix-Vector` product between ProjMPOSum2 and input ITensor `v`.
"""
(P::ProjMPOSum2)(v::ITensor) = product(P, v)

#################################################################################

Base.size(P::ProjMPOSum2) = Base.size(P.PHs[begin])

#################################################################################

"""
    position!(P::ProjMPOSum2, psi::MPS, k::Int)::ProjMPOSum2

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjMPOSum2, psi::MPS, pos::Int)::ProjMPOSum2
    @threaded_loop for p in P.PHs
        position!(p, psi, pos)
    end
    return P
end

#################################################################################

function ITensorMPS.noiseterm(P::ProjMPOSum2, phi::ITensor, ortho::String)
    nt = ITensor()
    using_threaded_loop() && (mutex = Threads.SpinLock())
    @threaded_loop for p in P.PHs
        temp_nt = noiseterm(p, phi, ortho)

        if using_threaded_loop()
            lock(mutex) do
                nt += temp_nt
            end
        else
            nt += temp_nt
        end
    end
    return nt
end

#################################################################################

