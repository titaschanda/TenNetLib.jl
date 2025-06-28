
#################################################################################

"""
    struct CouplingModel{T}
        sites::Vector{Index{T}}
        terms::Vector{IDTensors}
    end

`CouplingModel` for a given `OpStrings` Hamiltonian terms and `sites::Vector{Index}`.
 - `sites::Vector{Index{T}}`: Site `Index`s.
 - `terms::Vector{IDTensors}`: Collection of Hamiltonian terms. 
"""
struct CouplingModel{T}
    sites::Vector{Index{T}}
    terms::Vector{IDTensors}
end

#################################################################################

"""
    function Base.length(model::CouplingModel)

Number of sites in the system.
"""
Base.length(model::CouplingModel) = length(model.sites)

#################################################################################

"""
    function Base.copy(model::CouplingModel)

Shallow copy of `CouplingModel`.
"""
Base.copy(model::CouplingModel) = CouplingModel(Base.copy(model.sites), Base.copy(terms))

#################################################################################

"""
    function Base.getindex(model::CouplingModel, n)

Returns Hamiltonian terms for `n`th site(s) in the `CouplingModel`.
"""
Base.getindex(model::CouplingModel, n) = model.terms[n]

#################################################################################

"""
    function ITensors.siteinds(model::CouplingModel)

Returns the site `Index`s of the `CouplingModel`.
"""
ITensors.siteinds(model::CouplingModel) = model.sites

#################################################################################
"""
    _chunksum_cm_terms(cmterms::Vector2{ITensor},
                       chunksize::Int,
                       maxdim::Int,
                       mindim::Int,
                       cutoff::Float64,
                       svd_alg::String)::Vector{ITensor}

Internal helper function to iteratively sum and compress CouplingModel terms.
"""
function _chunksum_cm_terms(cmterms::Vector2{ITensor},
                            chunksize::Int,
                            maxdim::Int,
                            mindim::Int,
                            cutoff::Float64,
                            svd_alg::String)::Vector{ITensor}
    
    if length(cmterms) == 1
        return cmterms[begin]
        
    elseif length(cmterms) <= chunksize
        summed_cmterms = _directsum(cmterms)
        numsites = length(summed_cmterms)
        
        for dir = ["left", "right"]
            locs = dir == "left" ? (1 : numsites - 1) : reverse(1 : numsites - 1) 
            for b in locs                

                phi = summed_cmterms[b] * summed_cmterms[b+1]
                uinds = uniqueinds(summed_cmterms[b],
                                   summed_cmterms[b+1])

                U, S, V = svd(phi, uinds;
                              maxdim=maxdim,
                              mindim=mindim,
                              cutoff=cutoff,
                              alg = svd_alg,
                              lefttags = "OpLink",
                              righttags = "OpLink")

                norm = ITensors.norm(S)
                S /= norm
                if dir == "left"
                    V *= S
                else
                    U *= S
                end                
                summed_cmterms[b] = U * norm^0.5
                summed_cmterms[b+1] = V * norm^0.5
            end
        end
        return summed_cmterms
    else
        chunks = _divide_by_chunksize(collect(1 : length(cmterms)),
                                      chunksize)
        new_cmterms = Vector2{ITensor}(undef, length(chunks))

        @threaded_loop for ii = 1 : length(chunks)          
            new_cmterms[ii] = _chunksum_cmterms(cmterms[chunks[ii]],
                                                chunksize,
                                                maxdim,
                                                mindim,
                                                cutoff,
                                                svd_alg)
        end

        return _chunksum_cmterms(new_cmterms,
                                 chunksize,
                                 maxdim,
                                 mindim,
                                 cutoff,
                                 svd_alg)
    end
end

#################################################################################

"""
    _initCouplingModel(os::OpStrings{T1},
                       sites::Vector{Index{T2}};
                       merge::Bool = true,
                       maxdim::Int = typemax(Int),
                       mindim::Int = 1,
                       cutoff::Float64 = Float64_threashold(),
                       svd_alg::String = "divide_and_conquer",
                       chunksize::Int = 12) where {T1 <: Number, T2}

Internal helper function for the construction of `CouplingModel`.

#### Named arguments and their default values:
 - `merge::Bool = true`.
 - `maxdim::Int = typemax(Int)`.
 - `mindim::Int = 1`.
 - `cutoff::Float64 = Float64_threashold()`.
 - `svd_alg::String = "divide_and_conquer"`.
 - `chunksize::Int = 12`.
"""
function _initCouplingModel(os::OpStrings{T1},
                            sites::Vector{Index{T2}};
                            merge::Bool = true,
                            maxdim::Int = typemax(Int),
                            mindim::Int = 1,
                            cutoff::Float64 = Float64_threashold(),
                            svd_alg::String = "divide_and_conquer",
                            chunksize::Int = 12
                            )::Vector{IDTensors} where {T1 <: Number, T2}
    
    numsites::Int = length(sites)
    collected_tensors = Dict{Vector{Int}, Vector2{ITensor}}()

    os = removeIdsZeros(os)
    os = bosonize(os, sites)
    os = mergeterms(os)
    
    for term in os
        
        if minsite(term) < 1 || maxsite(term) > numsites
            error("`CouplingModel()`: Out-of-bounds error in OpStrings !!")
        end
        
        coeff = coefficient(term)
        ops = operators(term)
        
        positions = Int[x.second for x in ops]
        tensors = ITensor[op(x.first, sites[x.second]) for x in ops]    
        signcoeff = sign(coeff)
        abscoeff = abs(coeff)^(1.0/length(term))
        tensors .*= abscoeff
        tensors[begin] *= signcoeff
        
        _add_oplinks!(tensors)        
        
        if !haskey(collected_tensors, positions)
            collected_tensors[positions] = Vector2{ITensor}()
            push!(collected_tensors[positions], tensors)
        else
            if length(positions) == 1
                collected_tensors[positions][begin] .+= tensors                
            else
                push!(collected_tensors[positions], tensors)
            end
        end
    end

    terms = Vector{IDTensors}(undef, numsites)
    for ii = 1 : numsites
        terms[ii] = IDTensors()
    end
    
    if merge
        termdict = Dict{Vector{Int}, Vector{ITensor}}()    
        for (positions, tensors) in collected_tensors
            
            if length(tensors) == 1
                termdict[positions] = tensors[begin]
            else
                termdict[positions] = _chunksum_cm_terms(tensors,
                                                         chunksize,
                                                         maxdim,
                                                         mindim,
                                                         cutoff,
                                                         svd_alg)
            end
        end
        
        for (positions, tensors) in termdict
            id = gen_rand_id()            
            for ii = 1 : length(positions)
                b = positions[ii]
                terms[b][id] = tensors[ii]
            end
        end
    else
        for positions in keys(collected_tensors)
            for tensors in collected_tensors[positions]
                id = gen_rand_id()            
                for ii = 1 : length(positions)
                    b = positions[ii]
                    terms[b][id] = tensors[ii]
                end
            end
        end
    end
    
    return terms
end

#################################################################################

"""
    function CouplingModel(os::OpStrings{T1},
                           sites::Vector{Index{T2}};
                           merge::Bool = true,
                           maxdim::Int = typemax(Int),
                           mindim::Int = 1,
                           cutoff::Float64 = Float64_threashold(),
                           svd_alg::String = "divide_and_conquer",
                           chunksize::Int = 12) where {T1 <: Number, T2}

Constructor of the `CouplingModel` from `os::OpStrings` and `sites::Vector{Index}`.

#### Named arguments and their default values:
 - `merge::Bool = true`. If `true`, merges all the terms having same spatial support resulting
   in larger virtual dimension. Otherwise, all the terms have virtual dimension one.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_Threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `chunksize::Int = 12`.

#### Example:

    os = OpStrings()    
    for j=1:N-1
        os += 1, "Sz" => j, "Sz" => j+1
        os += 0.5, "S+" => j, "S-" => j+1
        os += 0.5, "S-" => j, "S+" => j+1
    end    
    H = CouplingModel(os,sites)
"""
function CouplingModel(os::OpStrings{T1},
                       sites::Vector{Index{T2}};
                       merge::Bool = true,
                       maxdim::Int = typemax(Int),
                       mindim::Int = 1,
                       cutoff::Float64 = Float64_threshold(),
                       svd_alg::String = "divide_and_conquer",
                       chunksize::Int = 12
                       )::CouplingModel{T2} where {T1 <: Number, T2}

    terms = _initCouplingModel(os, sites;
                               merge,
                               maxdim,
                               mindim,
                               cutoff,
                               svd_alg,
                               chunksize
                               )
    return CouplingModel(sites, terms)
end

#################################################################################

"""
    function CouplingModel(os::OpStrings{T1},
                           mpo::MPO;
                           merge::Bool = true,
                           maxdim::Int = typemax(Int),
                           mindim::Int = 1,
                           cutoff::Float64 = Float64_threshold(),
                           svd_alg::String = "divide_and_conquer",
                           chunksize::Int = 12) where {T1 <: Number}

Constructor of the `CouplingModel` from `os::OpStrings` and `mpo::MPO`.
Resultant `CouplingModel` is the sum of `os` and the `mpo`.

#### Named arguments and their default values:
 - `merge::Bool = true`. If `true`, merges all the terms having same spatial support resulting
   in larger virtual dimension. Otherwise, all the terms have virtual dimension one.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_Threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `chunksize::Int = 12`.
"""
function CouplingModel(os::OpStrings{T1},
                       mpo::MPO;
                       merge::Bool = true,
                       maxdim::Int = typemax(Int),
                       mindim::Int = 1,
                       cutoff::Float64 = Float64_threshold(),
                       svd_alg::String = "divide_and_conquer",
                       chunksize::Int = 12
                       )::CouplingModel where {T1 <: Number}
    
    sites = dag.(siteinds(first, mpo; plev = 0))
    numsites = length(sites)
    terms = _initCouplingModel(os, sites;
                               merge,
                               maxdim,
                               mindim,
                               cutoff,
                               svd_alg,
                               chunksize
                               )

    replacetags!(linkinds, mpo, "Link", "OpLink")
    id = gen_rand_id()
    for b = 1 : numsites
        terms[b][id] = mpo[b]
    end
    
    return CouplingModel(sites, terms)
end

#################################################################################

"""
    function CouplingModel(mpos::MPO...)

Constructs `CouplingModel` from a collection of `MPO`s. Different MPO terms are contracted
in parallel. Useful for TTN codes.
"""
function CouplingModel(mpos::MPO...)::CouplingModel

    sites = dag.(siteinds(first, first(mpos); plev = 0))
    @assert all(mpo -> dag.(siteinds(first, mpo; plev = 0)) == sites, mpos)
    
    numsites = length(sites)
    terms::Vector{IDTensors} = Vector{IDTensors}(undef, numsites)
    for ii = 1 : numsites
        terms[ii] = IDTensors()
    end

    for mpo in mpos
        mpo = replacetags(linkinds, mpo, "Link", "OpLink")
        id = gen_rand_id()
        for b = 1 : numsites
            terms[b][id] = mpo[b]
        end
    end
    
    return CouplingModel(sites, terms)
end

#################################################################################

