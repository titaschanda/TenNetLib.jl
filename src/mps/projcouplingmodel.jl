
#################################################################################

"""
    mutable struct ProjCouplingModel
        lpos::Int
        rpos::Int
        nsite::Int
        M::CouplingModel
        LR::Vector{IDTensors}
    end

Holds Environments from CouplingModel.
"""
mutable struct ProjCouplingModel
    lpos::Int
    rpos::Int
    nsite::Int
    M::CouplingModel
    LR::Vector{IDTensors}
end

#################################################################################

"""
    ProjCouplingModel(M::CouplingModel)

Constructor of the `ProjCouplingModel`.
"""
ProjCouplingModel(M::CouplingModel) = ProjCouplingModel(0, length(M) + 1, 2, M, 
                                                        Vector{IDTensors}(undef, length(M)))

#################################################################################

"""
    Base.copy(M::ProjCouplingModel)

Shallow copy of `ProjCouplingModel`.
"""
Base.copy(M::ProjCouplingModel) = ProjCouplingModel(M.lpos, M.rpos, M.nsite,
                                                    Base.copy(M), Base.copy(LR))

#################################################################################

"""
    set_nsite!(P::ProjCouplingModel, nsite::Int)::ProjCouplingModel

Set `nsite` of the `ProjCouplingModel` object.
""" 
function set_nsite!(P::ProjCouplingModel, nsite::Int)::ProjCouplingModel
    P.nsite = nsite
    return P
end

#################################################################################

"""
    nsite(P::ProjCouplingModel)

Returns `nsite` of ProjCouplingModel
"""
nsite(P::ProjCouplingModel) = P.nsite

#################################################################################

"""
    site_range(P::ProjCouplingModel)

Returns the range of sites where the Environments are acting. 
"""
site_range(P::ProjCouplingModel) = (P.lpos + 1):(P.rpos - 1)

#################################################################################

"""
    Base.length(P::ProjCouplingModel)

Returns the length of the underlying Environment.
"""
Base.length(P::ProjCouplingModel) = length(P.M)

#################################################################################

"""
    lproj(P::ProjCouplingModel)::IDTensors

Returns the Left Environment `IDTenors`.
"""
function lproj(P::ProjCouplingModel)::IDTensors
    (P.lpos <= 0) && return IDTensors()
    return P.LR[P.lpos]
end

#################################################################################

"""
    rproj(P::ProjCouplingModel)::IDTensors

Returns the Right Environment `IDTenors`.
"""
function rproj(P::ProjCouplingModel)::IDTensors
    (P.rpos >= length(P) + 1) && return IDTensors()
    return P.LR[P.rpos]
end

#################################################################################

"""
    Base.eltype(P::ProjCouplingModel)

Returns the element type (e.g., `Float64` or `ComplexF64`) of `ProjCouplingModel`.
"""
function Base.eltype(P::ProjCouplingModel)
    elt = eltype(lproj(P))
    for j in site_range(P)
        elt = promote_type(elt, eltype(P.M[j]))
    end
    return promote_type(elt, eltype(rproj(P)))
end

#################################################################################

function _makeL!(P::ProjCouplingModel, psi::MPS, k::Int)::Union{IDTensors, Nothing}
    ll = P.lpos
    if ll ≥ k
        P.lpos = k
        return nothing
    end
    
    ll = max(ll, 0)
    L = lproj(P)
    while ll < k
        
        P.LR[ll + 1] = IDTensors()
        idkeys = collect(union(keys(L), keys(P.M[ll+1])))

        local_tensor = ITensor()

        using_threaded_loop() && (mutex = Threads.SpinLock())        
        @threaded_loop for id in idkeys
            
            if haskey(L, id) && haskey(P.M[ll+1], id) 
                phi = psi[ll+1]
                phidag = dag(prime(phi))
                phidag *= L[id]
                phidag *= P.M[ll+1][id]
                phidag *= phi
                
            elseif haskey(L, id) || haskey(P.M[ll+1], id)
                phi = psi[ll+1]
                uncommon_tensor = haskey(L, id) ? L[id] : P.M[ll+1][id]
                ind1 = commonind(phi, uncommon_tensor)
                ind2 = commonind(phi, psi[ll+2])

                # Workaround for 0-site environment
                if isnothing(ind2)
                    ind2 = ll == 0 ? uniqueind(phi, uncommon_tensor; tags = "Link") :
                        uniqueind(phi, psi[ll]; tags = "Link")
                end
                
                phidag = dag(prime(phi, [ind1, ind2]))
                phidag *= uncommon_tensor;
                phidag *= phi;
            else
                error("SOMETHING IS WRONG")            
            end

            if using_threaded_loop()
                lock(mutex) do
                    if order(phidag) > 2
                        P.LR[ll+1][id] = phidag
                    else
                        local_tensor += phidag
                    end                    
                end
            else
                if order(phidag) > 2
                    P.LR[ll+1][id] = phidag
                else
                    local_tensor += phidag
                end
            end
            
        end

        if order(local_tensor) != 0
            newid = gen_rand_id()
            P.LR[ll + 1][newid] = local_tensor
        end
        L = P.LR[ll + 1]
        ll += 1
    end

    P.lpos = k
    return L
end

#################################################################################

"""
    makeL!(P::ProjCouplingModel, psi::MPS, k::Int)::ProjCouplingModel

Compute Left Environments for the MPS `psi` at position `k`.
"""
function makeL!(P::ProjCouplingModel, psi::MPS, k::Int)::ProjCouplingModel
    _makeL!(P, psi, k)
    return P
end

#################################################################################

function _makeR!(P::ProjCouplingModel, psi::MPS, k::Int)::Union{IDTensors, Nothing}
    rl = P.rpos
    if rl ≤ k
        P.rpos = k
        return nothing
    end
    
    N = length(P)
    rl = min(rl, N + 1)
    R = rproj(P)
    while rl > k
        
        P.LR[rl - 1] = IDTensors()
        idkeys = collect(union(keys(R), keys(P.M[rl-1])))
        
        local_tensor = ITensor()

        using_threaded_loop() && (mutex = Threads.SpinLock())
        @threaded_loop for id in idkeys
            
            if haskey(R, id) && haskey(P.M[rl-1], id) 
                phi = psi[rl-1]
                phidag = dag(prime(phi))
                phidag *= R[id]
                phidag *= P.M[rl-1][id]
                phidag *= phi
                
            elseif haskey(R, id) || haskey(P.M[rl-1], id)
                phi = psi[rl-1]
                uncommon_tensor = haskey(R, id) ? R[id] : P.M[rl-1][id]
                ind1 = commonind(phi, uncommon_tensor)
                ind2 = commonind(phi, psi[rl-2])

                # Workaround for 0-site environment
                if isnothing(ind2)
                    ind2 = rl == N+1 ? uniqueind(phi, uncommon_tensor; tags = "Link") :
                        uniqueind(phi, psi[rl]; tags = "Link")
                end
                
                phidag = dag(prime(phi, [ind1, ind2]))     
                phidag *= uncommon_tensor;
                phidag *= phi;
                
            else
                error("SOMETHING IS WRONG")            
            end

            if using_threaded_loop()
                lock(mutex) do
                    if order(phidag) > 2
                        P.LR[rl-1][id] = phidag
                    else
                        local_tensor += phidag
                    end
                end
            else
                if order(phidag) > 2
                    P.LR[rl-1][id] = phidag
                else
                    local_tensor += phidag
                end
            end
        end

        
        if order(local_tensor) != 0
            newid = gen_rand_id()
            P.LR[rl-1][newid] = local_tensor
        end
        R = P.LR[rl-1]
        rl -= 1
    end
    P.rpos = k
    return R
end

#################################################################################

"""
    makeR!(P::ProjCouplingModel, psi::MPS, k::Int)::ProjCouplingModel

Compute Right Environments for the MPS `psi` at position `k`.
"""
function makeR!(P::ProjCouplingModel, psi::MPS, k::Int)::ProjCouplingModel
    _makeR!(P, psi, k)
    return P
end

#################################################################################

"""
    position!(P::ProjCouplingModel, psi::MPS, pos::Int)::ProjCouplingModel

Compute Left and Right Environments for the MPS `psi` at position `pos`.
"""
function position!(P::ProjCouplingModel, psi::MPS, pos::Int)::ProjCouplingModel
    makeL!(P, psi, pos - 1)
    makeR!(P, psi, pos + nsite(P))
    return P
end

#################################################################################

function _contract(P::ProjCouplingModel, v::ITensor)::ITensor
    
    idtensor_map = IDTensors[lproj(P)]
    append!(idtensor_map, P.M[site_range(P)])
    push!(idtensor_map, rproj(P))

    if isnothing(first(idtensor_map))
        reverse!(idtensor_map)
    end
    
    idtens = Dict{IDType, Vector{ITensor}}()
    for it in idtensor_map
        for term in it
            if !haskey(idtens, term.first)
                idtens[term.first] = ITensor[term.second]
            else
                push!(idtens[term.first], term.second)
            end
        end
    end
    
    sum_tensor = ITensor()    
    idkeys = collect(keys(idtens))

    using_threaded_loop() && (mutex = Threads.SpinLock())
    @threaded_loop for id in idkeys
        
        Hv = ITensors.contract(v, idtens[id]...;
                               sequence = ITensors.default_sequence())
        
        if using_threaded_loop()
            lock(mutex) do
                sum_tensor += noprime(Hv)
            end
        else
            sum_tensor += noprime(Hv)
        end
    end

    
    return sum_tensor
end

#################################################################################

"""
    product(P::ProjCouplingModel, v::ITensor)::ITensor

Compute the `Matrix-Vector` product between ProjCouplingModel and input ITensor `v`.
"""
function product(P::ProjCouplingModel, v::ITensor)::ITensor
    Pv = _contract(P, v)
    if order(Pv) != order(v)
        error(
            string(
                "The order of the ProjCouplingModel-ITensor product P*v ", 
                "is not equal to the order of the ITensor v, ",
                "this is probably due to an index mismatch.\nCommon reasons for this error: \n",
                "(1) You are trying to multiply the ProjCouplingModel with the $(nsite(P))-site ",
                "wave-function at the wrong position.\n",
                "(2) `orthogonalize!` was called, changing the MPS without updating ",
                "the ProjMPO.\n\n",
                "P*v inds: $(inds(Pv)) \n\n",
                "v inds: $(inds(v))",
            ),
        )
    end
    return Pv
end

#################################################################################

(P::ProjCouplingModel)(v::ITensor) = product(P, v)

#################################################################################

function ITensorMPS.noiseterm(P::ProjCouplingModel, phi::ITensor, ortho::String)::ITensor
    if nsite(P) != 2
        error("noise term only defined for 2-site ProjMPO")
    end
    
    site_range_P = site_range(P)
    if ortho == "left"        
        idtensor_map = IDTensors[lproj(P)]
        push!(idtensor_map, P.M[first(site_range_P)])
        if isnothing(first(idtensor_map))
            reverse!(idtensor_map)
        end
        
        idtens = Dict{IDType, Vector{ITensor}}()
        for it in idtensor_map
            for term in it
                if !haskey(idtens, term.first)
                    idtens[term.first] = ITensor[term.second]
                else
                    push!(idtens[term.first], term.second)
                end
            end
        end
        
        nt = ITensor()        
        idkeys = collect(keys(idtens))

        using_threaded_loop() && (mutex = Threads.SpinLock())
        @threaded_loop for id in idkeys
            
            phitemp = ITensors.contract(phi, idtens[id]...;
                                        sequence = ITensors.default_sequence())

            oplinkinds = inds(phitemp; tags = "OpLink")
            (length(oplinkinds) != 0) && setprime!(phitemp, 0, oplinkinds)
            linds = inds(phitemp; tags = "Link,l=$(first(site_range_P)-1)")
            (length(linds) != 0) && setprime!(phitemp, 1, linds)
            sinds = inds(phitemp; tags = "Site,n=$(first(site_range_P))")
            (length(sinds) != 0) && setprime!(phitemp, 1, sinds)
            
            phitemp *= dag(noprime(phitemp))          

            if using_threaded_loop()
                lock(mutex) do
                    nt +=  phitemp
                end
            else
                nt +=  phitemp
            end
            
        end
        
    elseif ortho == "right"
        idtensor_map = IDTensors[rproj(P)]
        push!(idtensor_map, P.M[last(site_range_P)])
        if isnothing(first(idtensor_map))
            reverse!(idtensor_map)
        end
        
        idtens = Dict{IDType, Vector{ITensor}}()
        for it in idtensor_map
            for term in it
                if !haskey(idtens, term.first)
                    idtens[term.first] = ITensor[term.second]
                else
                    push!(idtens[term.first], term.second)
                end
            end
        end
        
        nt = ITensor()
        idkeys = collect(keys(idtens))

        using_threaded_loop() && (mutex = Threads.SpinLock())
        @threaded_loop for id in idkeys

            phitemp = ITensors.contract(phi, idtens[id]...;
                                        sequence = ITensors.default_sequence())

            oplinkinds = inds(phitemp; tags = "OpLink")
            (length(oplinkinds) != 0) && setprime!(phitemp, 0, oplinkinds)
            rinds = inds(phitemp; tags = "Link,l=$(last(site_range_P))")
            (length(rinds) != 0) && setprime!(phitemp, 1, rinds)
            sinds = inds(phitemp; tags = "Site,n=$(last(site_range_P))")
            (length(sinds) != 0) && setprime!(phitemp, 1, sinds)
            
            phitemp *= dag(noprime(phitemp))

            if using_threaded_loop()
                lock(mutex) do
                    nt +=  phitemp
                end
            else
                nt +=  phitemp
            end
        end
    else
        error("In noiseterm, got ortho = $ortho, only supports `left` and `right`")
    end

    return nt
end

#################################################################################
