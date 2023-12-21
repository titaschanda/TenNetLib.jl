
#################################################################################

function _parity_sign(P::Vector{Int})::Int
    L = length(P)
    s = +1
    for i in 1:L, j in (i + 1):L
        (P[j] != P[i]) && (s *= sign(P[j] - P[i]))
    end
    return s
end

#################################################################################

"""
    function bosonize(oppair::Vector{Pair{String, Int}},
                      sites::Vector{Index{T}}) where T

Given a string of operator names with positions in the form of `Vector{Pair{String, Int}}`
"bosonizes" the operator string by putting the Jordan-Wigner string "F" at appropriate places.

#### Arguments
 - `oppair::Vector{Pair{String, Int}}`: Input string of operator names with positions.
 - `sites::Vector{Index}`: The entire site `Index`s as per ITensors' convention.

#### Return values
 - `::Int`: Even (+1) or odd (-1) permutation.
 - `::Vector{Pair{String, Int}}`:  modified string of operator names and positions.
"""
function bosonize(oppair::Vector{Pair{String, Int}},
                  sites::Vector{Index{T}}
                  )::Tuple{Int, Vector{Pair{String, Int}}} where T

    oppair = copy(oppair)
    parity_pos = Vector{Int}()
    for pp in oppair
        if has_fermion_string(pp.first, sites[pp.second])
            push!(parity_pos, pp.second)
        end
    end

    #************************************************************

    fsize = length(parity_pos)

    if fsize % 2 != 0
        error("`bosonize()`: Cannot handle odd number of fermionic ops !!")
    end

    numperm = _parity_sign(parity_pos)
    
    #************************************************************

    if fsize != 0
        sort!(oppair, by = x -> x.second)
        minpos = oppair[begin].second
        maxpos = oppair[end].second    
        oppairNew = Vector{Pair{String, Int}}()
        fermion_on_left = false
        for ii = minpos : maxpos
            oppos = findall(x -> x.second == ii, oppair)
            if length(oppos) == 0 && fermion_on_left
                push!(oppairNew, "F" => ii)
            elseif length(oppos) != 0
                for pos in oppos
                    isfermi = has_fermion_string(oppair[pos].first, sites[oppair[pos].second])
                    isfermi && (fermion_on_left = !fermion_on_left)
                    push!(oppairNew, oppair[pos].first => ii)
                end
                fermion_on_left && push!(oppairNew, "F" => ii)
            end
        end
    else
        oppairNew = oppair
    end

    #************************************************************

    merged_dict = Dict{Int, String}()

    for (value, key) in oppairNew
        if haskey(merged_dict, key)
            merged_dict[key] = string(merged_dict[key], " * ", value)
        else
            merged_dict[key] = value
        end
    end

    oppairNew =
        Pair{String, Int}[value => key for (key, value) in merged_dict]
    sort!(oppairNew, by = x -> x.second)
    
    #************************************************************
      
    return numperm, oppairNew
end

#################################################################################
