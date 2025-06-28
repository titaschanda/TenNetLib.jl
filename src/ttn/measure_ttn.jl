
#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::TTN, opten::ITensor)

Returns (complex / real) local expectation value (`::ComplexF64` / `::Float64`)
of a given TTN `psi::TTN`. The operator `opten::ITensor` must be single-site operator.

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::TTN,
                 opten::ITensor)::T where T <: Union{ComplexF64, Float64}
    pos = findsite(psi, opten)
    pass = hasind(opten, siteind(psi, pos)) &&
        hasind(opten, siteind(psi, pos)') &&
        ndims(opten) == 2
    if !pass
        error("`expectC()`: Error in operator tensors !!")
    end
    removeqn = !hasqns(opten) && hasqns(psi) ? true : false
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A

    node = find_sitenode(psi, pos)
    isometrize!(psi, node)
    ket = conditional_removeqns(psi[node])
    ind = commonind(ket, opten)
    bra = dag(prime(ket, ind))
    val = complex(scalar(bra * opten * ket))
    
    if T <: Complex
        return val
    else
        if abs(imag(val)) > 100 * Float64_threshold()
            @warn "`measure(::Float64)`: Got complex number with |imaginary part| > $(100 * Float64_threshold()) !!"
        end
        return real(val)
    end
    
    return complex(val)
end

measure(psi::TTN, opten::ITensor) = measure(ComplexF64, psi, opten)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::TTN, opstr::String, pos::Int)

Returns (complex / real) local expectation value (`::ComplexF64` / `::Float64`) of a
given TTN `psi::TTN` for a given operator name (`opstr::String`) and a site index (`pos::Int`) 

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
measure(::Type{T}, psi::TTN, opstr::String,
        pos::Int) where T <: Union{ComplexF64, Float64} =
            measure(T, psi, op(opstr, siteind(psi, pos)))

measure(psi::TTN, opstr::String, pos::Int) =
    measure(ComplexF64, psi, opstr, Int)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::TTN, opstr::String; kwargs...)

Returns (complex / real) local expectation values (`::Vector{ComplexF64}` / `::Vector{Float64}`)
of a given TTN `psi::TTN` for a given operator name (`opstr::String`) at all the sites.

Optionally, for specific sites, keyword argument `sites` can be specified, e.g.,
`sites = [1, 2, 3]`.

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::TTN, opstr::String;
                 kwargs...) where T <: Union{ComplexF64, Float64}
    N = numsites(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> measure(T, psi, opstr, pos), sites)
end

measure(psi::TTN, opstr::String; kwargs...) =
    measure(ComplexF64, psi, opstr; kwargs)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::TTN, optens::Vector{ITensor})
    
Returns (complex / real) multi-site expectation value (`::ComplexF64` / `::Float64`)
of a given TTN `psi::TTN` for a given vector of single-site operators
(`optens::Vector{ITensor}`).

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::TTN,
                 optens::Vector{ITensor})::T where T <: Union{ComplexF64, Float64}

    sitenodes = Set{Int2}()
    for o in optens
        push!(sitenodes, findnode(psi, o))
    end

    oc = find_eccentric_central_node(psi.graph, sitenodes)
    isometrize!(psi, oc)

    removeqn = false
    for o in optens
        !hasqns(o) && hasqns(psi) && (removeqn = true)
    end
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A    

    lnt = LinkTensorsTTN(psi, optens)
    phi = conditional_removeqns(psi[oc])
    val = complex(scalar(product(lnt, psi, phi) * dag(phi)))

    if T <: Complex
        return val
    else
        if abs(imag(val)) > 100 * Float64_threshold()
            @warn "`measure(::Float64)`: Got complex number with |imaginary part| > $(100 * Float64_threshold()) !!"
        end
        return real(val)
    end
    
    return complex(val)
end

measure(psi::TTN, optens::Vector{ITensor}) = measure(ComplexF64, psi, optens)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::TTN, oppairs::Vector{Pair{String, Int}};
                     isfermions::Bool = true)

Returns (complex / real) multi-site expectation value (`::ComplexF64` / `::Float64`) of
a given TTN `psi::TTN`. `oppairs::Vector{Pair{String, Int}}` contains pairs of operator
names (`String`) and site locations (`Int`).
E.g., for <ψ|OᵢOⱼOₖ... |ψ>, `oppairs = ["O" => i, "O" => j, "O" => k,...]`.

Fermionic JW strings are added automatically for fermionic operators if
`isfermions::Bool = true` (default).

For `ElT = Float64`, if the expectation value is complex, raises a warning!

#### Example:
    valueC = measure(psi, ["Cdag" => 2, "C" => 6, "Cdag" => 9, "C" => 12])
    valueR = measure(Float64, psi, ["Cdag" => 2, "C" => 6, "Cdag" => 9, "C" => 12])
"""
function measure(::Type{T}, psi::TTN, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::T where T <: Union{ComplexF64, Float64}
    
    sites = siteinds(psi)
    if isfermions
        perm, oppairs = bosonize(oppairs, sites)
    end
    
    optens = ITensor[op(x.first, siteind(psi, x.second)) for x in oppairs]
    return isfermions ? perm * measure(T, psi, optens) : measure(T, psi, optens) 
end

measure(psi::TTN, oppairs::Vector{Pair{String, Int}}; isfermions::Bool = true) =
    measure(ComplexF64, psi, oppairs; isfermions)
#################################################################################

