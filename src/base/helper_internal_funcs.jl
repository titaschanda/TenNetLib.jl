
################################################################################

"""
    function gen_rand_id()

Generates a random id.
"""
gen_rand_id() = rand(ITensors.index_id_rng(), IDType)

#################################################################################

"""
    _divide_by_chunksize(vec::Vector{T}, size::Int)

Helper function to divide a vector by chunksize `size`. After the division,
remaining elements are put in another vector. Returns a Vector2{T}.
"""
function _divide_by_chunksize(vec::Vector{T}, size::Int)::Vector2{T} where T
    n = length(vec)
    numchunks, leftover = divrem(n, size)
    subvectors = Vector{T}[]
    for i in 1 : numchunks
        push!(subvectors, vec[(i-1)*size + 1 : i * size])
    end

    leftover == 1 && push!(subvectors, T[vec[end]])
    leftover > 1 && push!(subvectors, vec[numchunks * size + 1 : end])
    return subvectors
end

#################################################################################

"""
    _add_oplinks!(tensors::Vector{ITensor})

Adds virtual links to the operator strings.
"""
function _add_oplinks!(tensors::Vector{ITensor})::Nothing

    N = length(tensors)    
    if hasqns(tensors[1])
        lflux = QN()
        for j in 1:(N - 1)
            lflux += flux(tensors[j])
        end
        links = Vector{Index{QNBlocks}}(undef, N - 1)
        for j in (N - 1):-1:1
            links[j] = dag(Index(lflux => 1; tags="OpLink"))
            lflux -= flux(tensors[j])
        end
    else
        links = Index{Int}[Index(1; tags="OpLink") for n in 1:N]
    end

    for b = 1 : N - 1
        tensors[b] *= onehot(links[b] => 1)
        tensors[b+1] *= onehot(dag(links[b]) => 1)
    end
    return nothing
end

#################################################################################

"""
    _directsum(ψ⃗::Vector2{ITensor})

Internal function to perform directsum of vector of tensors (like summing MPS / MPO).
Neighboring tensors must have one shared index. Copied from ITensor.jl.
"""
function _directsum(ψ⃗::Vector2{ITensor})::Vector{ITensor}
    
    n = length(first(ψ⃗))
    @assert all(ψᵢ -> length(first(ψ⃗)) == length(ψᵢ), ψ⃗)

    # Output tensor
    ϕ = Vector{ITensor}(undef, n)

    # Direct sum first tensor
    j = 1
    l⃗j = map(ψᵢ -> commonind(ψᵢ[j], ψᵢ[j+1]), ψ⃗)
    ϕj, (lj,) = directsum((ψ⃗[i][j] => (l⃗j[i],) for i in 1:length(ψ⃗))...;
                          tags=[tags(first(l⃗j))])
    ljm_prev = lj
    ϕ[j] = ϕj
    
    for j in 2:(n - 1)
        l⃗jm = map(ψᵢ -> commonind(ψᵢ[j], ψᵢ[j-1]), ψ⃗)
        l⃗j = map(ψᵢ -> commonind(ψᵢ[j], ψᵢ[j+1]), ψ⃗)
        ϕj, (ljm, lj) = directsum((ψ⃗[i][j] => (l⃗jm[i], l⃗j[i]) for i in 1:length(ψ⃗))...;
                                  tags=[tags(first(l⃗jm)), tags(first(l⃗j))])
        ϕj = replaceind(ϕj, ljm => dag(ljm_prev))
        ljm_prev = lj
        ϕ[j] = ϕj
    end
    
    j = n
    l⃗jm = map(ψᵢ -> commonind(ψᵢ[j], ψᵢ[j-1]), ψ⃗)
    ϕj, (ljm,) = directsum((ψ⃗[i][j] => (l⃗jm[i],) for i in 1:length(ψ⃗))...;
                           tags=[tags(first(l⃗jm))])
    ϕj = replaceind(ϕj, ljm => dag(ljm_prev))
    ϕ[j] = ϕj
    return ϕ
end

#################################################################################

"""
    function combineinds(inds::Vector{Index};
                         maxdim::Union{Nothing, Int} = nothing, 
                         maxqnblocks::Union{Nothing, Int} = nothing,
                         kwargs...)

Combine a vector of `Index` into one (like ITensors.jl's `combiner`). `maxdim` is be the maximum
dimension of the output `Index`, `maxqnblocks` represents maximum number of QN blocks to retain
in the output `Index`.
"""
function combineinds(inds::Vector{Index{Int}};
                     maxdim::Union{Nothing, Int} = nothing, 
                     kwargs...)::Index{Int}
    sumdim = sum(x -> dim(x), inds)     
    !isnothing(maxdim) && sumdim > maxdim && (sumdim = maxdim)
    return Index(sumdim; kwargs...)
end

function combineinds(inds::Vector{Index{QNBlocks}};
                     maxdim::Union{Nothing, Int} = nothing, 
                     maxqnblocks::Union{Nothing, Int} = nothing,
                     kwargs...)::Index{QNBlocks}    

    combind = combinedind(combiner(inds; kwargs...))

    if isnothing(maxdim) && isnothing(maxqnblocks)
        return combind
    else
        combblocks = space(combind)

        if !isnothing(maxqnblocks)
            sort!(combblocks, lt = (x, y) -> x.second > y.second)        
            maxqnblocks < length(combblocks) && (combblocks = combblocks[1:maxqnblocks])
        end

        if !isnothing(maxdim)
            blockdims::Vector{Int} = Int[x.second for x in combblocks]        
            if sum(blockdims) > maxdim
                blockdims = Int.(round.(maxdim * (blockdims / sum(blockdims))))            
                for ii = 1 : length(combblocks)           
                    combblocks[ii] = combblocks[ii].first => blockdims[ii]
                end
                combblocks = filter(x -> x.second > 0, combblocks)
            end
        end
        return Index(combblocks; dir=dir(combind), kwargs...)
    end
end

#################################################################################

"""
    function indexintersection(inds1::Vector{Index}, inds2::Vector{Index};
                               maxdim::Union{Nothing, Int} = nothing,
                               maxqnblocks::Union{Nothing, Int} = nothing,
                               kwargs...)

Performs set intersection of two vectors of `Index`. `maxdim` is be the maximum
dimension of the output `Index`, `maxqnblocks` represents maximum number of QN blocks to retain
in the output `Index`.
"""
function indexintersection(inds1::Vector{Index{Int}}, inds2::Vector{Index{Int}};
                           maxdim::Union{Nothing, Int} = nothing,
                           kwargs...)::Index{Int}
    sumdim = min(sum(x -> dim(x), inds1), sum(x -> dim(x), inds2))
    !isnothing(sumdim) && sumdim > maxdim && (sumdim = maxdim)
    return Index(sumdim; kwargs...)
end

function indexintersection(inds1::Vector{Index{QNBlocks}},
                           inds2::Vector{Index{QNBlocks}};
                           maxdim::Union{Nothing, Int} = nothing,
                           maxqnblocks::Union{Nothing, Int} = nothing,
                           kwargs...)::Index{QNBlocks}    
    
    combind1 = combinedind(combiner(inds1))
    combind2 = combinedind(combiner(inds2))

    combblocks1 = space(combind1)
    combblocks2 = space(combind2)

    sameQN(qn1::QN, qn2::QN) = dir(combind1) == dir(combind2) ?
        qn1 == qn2 : qn1 + qn2 == QN()
            
    newblocks = QNBlocks()
    for block in combblocks1
        loc = findall(x -> sameQN(block.first, x.first), combblocks2)        
        if length(loc) > 1
            error("`indexintersection()`: Multiple QN blocks with same QN !!")
        elseif length(loc) == 1
            blockdim = min(block.second, combblocks2[loc[1]].second)
            !isnothing(maxdim) && (blockdim = min(blockdim, maxdim))
            push!(newblocks, block.first => blockdim)
        end
    end
    if length(newblocks) == 0
        error("`indexintersection()`: No common QN blocks present !!")
    end

    if !isnothing(maxqnblocks)
        sort!(newblocks, lt = (x, y) -> x.second > y.second)      
        maxqnblocks < length(newblocks) && (newblocks = newblocks[1:maxqnblocks])
    end

    if !isnothing(maxdim)
        blockdims::Vector{Int} = Int[x.second for x in newblocks]
        
        if sum(blockdims) > maxdim
            blockdims = Int.(round.(maxdim * (blockdims / sum(blockdims))))
            for ii = 1 : length(newblocks)           
                newblocks[ii] = newblocks[ii].first => blockdims[ii]
            end
            newblocks = filter(x -> x.second > 0, newblocks)
        end
    end
    return Index(newblocks; dir = dir(combind1), kwargs...)
end

#################################################################################

