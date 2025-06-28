
#################################################################################

"""
    struct OpString{T <: Number}
        coeff::T
        ops::Vector{Pair{String, Int}}
    end

Holds operator strings (operator names with corresponding positions) along with
the coefficient.
 - `coeff::T`: Coeffcient of the operator string.
 - `ops::Vector{Pair{String, Int}}`: String of operator names along with the positions.
"""
struct OpString{T <: Number}
    coeff::T
    ops::Vector{Pair{String, Int}}
end

#################################################################################

"""
    function Base.copy(opstr::OpString)

Shallow copy of `OpString`.
"""
Base.copy(opstr::OpString) = OpString(Base.copy(opstr.coeff), Base.copy(opstr.ops))

#################################################################################

"""
    function coefficient(opstr::OpString)

Returns the coefficient of the `OpString`.
"""
coefficient(opstr::OpString)::Number =
    abs(imag(opstr.coeff)) < Float64_threshold() ? real(opstr.coeff) : opstr.coeff

#################################################################################

"""
    function operators(opstr::OpString)

Returns the operator string of the `OpString`.
"""
operators(opstr::OpString) = opstr.ops

#################################################################################

"""
    function minsite(opstr::OpString)

Returns the lowest site position in the operator string of the `OpString`.
"""
minsite(opstr::OpString) = minimum(x -> x.second, operators(opstr))

#################################################################################

"""
    function maxsite(opstr::OpString)

Returns the highest site position in the operator string of the `OpString`.
"""
maxsite(opstr::OpString) = maximum(x -> x.second, operators(opstr))

#################################################################################

"""
    function Base.length(opstr::OpString)

Returns the length of the operator string in the `OpString`.
"""
Base.length(opstr::OpString) = Base.length(opstr.ops)

#################################################################################

"""
    function isless(opstr1::OpString{T1},
                    opstr2::OpString{T2}) where {T1 <: Number, T2 <: Number}

Comparisions between two `OpStrings` required for sorting.
"""
function isless(opstr1::OpString{T1},
                opstr2::OpString{T2}) where {T1 <: Number, T2 <: Number}
    ops1 = operators(opstr1)
    ops2 = operators(opstr2)

    if ops1[begin].second != ops2[begin].second
        return ops1[begin].second < ops2[begin].second
    else
        return length(ops1) < length(ops2)
    end
end

#################################################################################

"""
    function removeIds(opstr::OpString{T}) where {T <: Number}

Returns an `OpString` with all "Id" operators removed from the original.
"""
function removeIds(opstr::OpString{T})::OpString{T} where {T <: Number}
    coeff = opstr.coeff
    ops = opstr.ops

    hasId = any(x -> x.first == "Id", ops)
    if hasId
        newops = filter(x -> x.first != "Id", ops)
        return OpString(coeff, newops)
    else
        return opstr
    end
end

#################################################################################

"""
    function bosonize(opstr::OpString{T1},
                      sites::Vector{Index{T2}}) where {T1 <: Number, T2}

Returns an `OpString` after "bosonizing" the original with Jordan-Wigner strings
as needed. See [`bosonize`](@ref).
"""
function bosonize(opstr::OpString{T1},
                  sites::Vector{Index{T2}}
                  )::OpString{T1} where {T1 <: Number, T2}    
    coeff = opstr.coeff
    numperm, ops = bosonize(opstr.ops, sites)
    return OpString(numperm * coeff, ops)
end
    

#################################################################################

"""
    const OpStrings{T} = Vector{OpString{T}}

Collection of `OpString`s.

#### Syntax:
    os = OpStrings()
    os += 1, "Sx" => i, "Sx" => j, "Sx" => k, ....
    os += "Sx" => i, "Sx" => j, "Sx" => k, ....

#### Example:
    os = OpStrings()    
    for j=1:N-1
        os += 1, "Sz" => j, "Sz" => j+1
        os += 0.5, "S+" => j, "S-" => j+1
        os += 0.5, "S-" => j, "S+" => j+1
    end
"""
const OpStrings{T} = Vector{OpString{T}}

#################################################################################

"""
    OpStrings()

Constructor for empty `OpStrings`.
"""
OpStrings() = Vector{OpString{ComplexF64}}()

#################################################################################

Base.:+(os::OpStrings{T}, 
        term::Tuple{Number, Vararg{Pair{String, Int}}}) where T <: Number = 
            push!(os, OpString(T(term[1]), collect(term[2:end])))

Base.:+(os::OpStrings{T}, 
        term::Tuple{Vararg{Pair{String, Int}}}) where T <: Number = 
            push!(os, OpString(T(1.0), collect(term[1:end])))

Base.:+(os::OpStrings{T}, 
        term::Pair{String, Int}) where T <: Number = 
            push!(os, OpString(T(1.0), typeof(term)[term]))

Base.:+(os::OpStrings{T}, 
        coeff::Number, 
        ops::Vector{Pair{String, Int}}) where T <: Number = 
            push!(os, OpString(T{coeff}, ops))

Base.:+(os::OpStrings{T}, 
        ops::Vector{Pair{String, Int}}) where T <: Number = 
            push!(os, OpString(T{1.0}, ops))

#################################################################################

"""
    function removeIdsZeros(os::OpStrings{T}) where {T <: Number}

Returns an `OpStrings` with all "Id" operators removed from the original,
as well as any `OpString` term that has `coeff` less than `Float64_threshold()`.
"""
function removeIdsZeros(os::OpStrings{T})::OpStrings{T} where T <: Number
    osnew = removeIds.(os)
    return filter(x -> abs(coefficient(x)) > Float64_threshold(), osnew)
end

#################################################################################

"""
    function bosonize(os::OpStrings{T1}, sites::Vector{Index{T2}}) where {T1 <: Number, T2}

Returns an `OpStrings` after "bosonizing" the original with Jordan-Wigner strings
as needed. See [`bosonize`](@ref).
"""
function bosonize(os::OpStrings{T1},
                  sites::Vector{Index{T2}}
                  )::OpStrings{T1} where {T1 <: Number, T2}
    return bosonize.(os, Ref(sites))
end

#################################################################################

"""
    function mergeterms(os::OpStrings{T}) where T <: Number

Returns an `OpStrings` where `OpString` elements with exactly same operator
strings has been merged by adding the coefficients.  
"""
function mergeterms(os::OpStrings{T})::OpStrings{T} where T <: Number   
    merged_dict = Dict{Vector{Pair{String, Int}}, T}()
    for term in os
        if haskey(merged_dict, operators(term))
            merged_dict[operators(term)] += term.coeff
        else
            merged_dict[operators(term)] = term.coeff
        end
    end
    
    osnew = OpString{T}[OpString(coeff, ops) for (ops, coeff) in merged_dict]
    sort!(osnew, lt = isless)    
    return osnew
end

#################################################################################

