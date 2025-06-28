
"""
SUPER INEFFICIENT CODES FOR CREATING `MPO` FROM `OpStrings`.
"""

#################################################################################

"""
    _chunksum_mpoterms(mpoterms::Vector2{ITensor},
                       chunksize::Int,
                       maxdim::Int,
                       mindim::Int,
                       cutoff::Float64,
                       svd_alg::String)::Vector{ITensor}

Internal helper function to iteratively sum and compress MPO terms.
"""
function _chunksum_mpoterms(mpoterms::Vector2{ITensor},
                            chunksize::Int,
                            maxdim::Int,
                            mindim::Int,
                            cutoff::Float64,
                            svd_alg::String)::Vector{ITensor}
    
    if length(mpoterms) == 1
        return mpoterms[begin]
        
    elseif length(mpoterms) <= chunksize
        
        summed_mpoterms = _directsum(mpoterms)       
        numsites = length(summed_mpoterms)
        
        for dir = ["left", "right"]
            locs = dir == "left" ? (1 : numsites - 1) : reverse(1 : numsites - 1) 
            for b in locs                

                phi = summed_mpoterms[b] * summed_mpoterms[b+1]
                uinds = uniqueinds(summed_mpoterms[b],
                                   summed_mpoterms[b+1])
                U, S, V = svd(phi, uinds;
                              maxdim=maxdim,
                              mindim=mindim,
                              cutoff=cutoff,
                              alg = svd_alg,
                              lefttags = "OpLink,l=$b",
                              righttags = "OpLink,l=$b")

                norm = ITensors.norm(S)
                S /= norm
                if dir == "left"
                    V *= S
                else
                    U *= S
                end                
                summed_mpoterms[b] = U * norm^0.5
                summed_mpoterms[b+1] = V * norm^0.5
            end
        end
        return summed_mpoterms
    else
        chunks = _divide_by_chunksize(collect(1 : length(mpoterms)),
                                      chunksize)
        new_mpoterms = Vector2{ITensor}(undef, length(chunks))

        @threaded_loop for ii = 1 : length(chunks)          
            new_mpoterms[ii] = _chunksum_mpoterms(mpoterms[chunks[ii]],
                                                  chunksize,
                                                  maxdim,
                                                  mindim,
                                                  cutoff,
                                                  svd_alg)
        end
                
        return _chunksum_mpoterms(new_mpoterms,
                                  chunksize,
                                  maxdim,
                                  mindim,
                                  cutoff,
                                  svd_alg)
    end
end

#################################################################################

"""
    function ITensorMPS.MPO(os::OpStrings{T1},
                            sites::Vector{Index{T2}};
                            maxdim::Int = typemax(Int),
                            mindim::Int = 1,
                            cutoff::Float64 = Float64_threshold(),
                            svd_alg::String = "divide_and_conquer",
                            chunksize::Int = 12) where {T1 <: Number, T2}

Creates `MPO` from `os::OpStrings`.
The present version uses recursive SVDs to create the MPO. Very inefficient when number of
Hamiltonian terms is large. Future updates will solve the problem.

#### Named arguments and their default values:
 - `maxdim::Int = typemax(Int)`: Maximum MPO bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum MPO bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `chunksize::Int = 12`. Maximum size of the chunks on which recursive SVDs are performed.
"""
function ITensorMPS.MPO(os::OpStrings{T1},
                      sites::Vector{Index{T2}};
                      maxdim::Int = typemax(Int),
                      mindim::Int = 1,
                      cutoff::Float64 = Float64_threshold(),
                      svd_alg::String = "divide_and_conquer",
                      chunksize::Int = 12
                      )::MPO where {T1 <: Number, T2}
    
    numsites::Int = length(sites)
    os = removeIdsZeros(os)
    os = bosonize(os, sites)
    os = mergeterms(os)
    
    for term in os 
        if minsite(term) < 1 || maxsite(term) > numsites
            error("`MPO()`: Out-of-bounds error in OpStrings !!")
        end
    end
    
    numterms = length(os)
    chunks = _divide_by_chunksize(collect(1:numterms), chunksize)
    summed_mpoterms = Vector2{ITensor}(undef, length(chunks))


    @threaded_loop for dummy = 1 : length(chunks)

        collected_mpoterms::Vector2{ITensor} = []        
        for ii in chunks[dummy]
            term = os[ii]
            coeff = coefficient(term)
            ops = operators(term)
            
            positions = Int[x.second for x in ops]
            tensors = ITensor[op(x.first, sites[x.second]) for x in ops]    
            signcoeff = sign(coeff)
            abscoeff = abs(coeff)^(1.0/length(term))
            tensors .*= abscoeff
            tensors[begin] *= signcoeff
            
            single_mpoterm = Vector{ITensor}(undef, numsites)
            single_mpoterm[positions] = tensors
            
            for b = 1 : numsites
                !isassigned(single_mpoterm, b) && (single_mpoterm[b] = op("Id", sites[b]))
            end

            _add_oplinks!(single_mpoterm)            
            
            push!(collected_mpoterms, single_mpoterm)
        end
        
        summed_mpoterms[dummy] = _chunksum_mpoterms(collected_mpoterms,
                                                    chunksize,
                                                    maxdim,
                                                    mindim,
                                                    cutoff,
                                                    svd_alg)        
    end
    
    final_mpoterms = _chunksum_mpoterms(summed_mpoterms,
                                        chunksize,
                                        maxdim,
                                        mindim,
                                        cutoff,
                                        svd_alg)
    return MPO(final_mpoterms)
end

#################################################################################

