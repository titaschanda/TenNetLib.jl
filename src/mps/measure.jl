
#################################################################################

"""
    function shanon_entropy(p::Vector{Float64})

Compute Shannon Entropy of a given vector `p`.
"""
function shannon_entropy(p::Vector{Float64})::Float64

    norm = sum(p)
    if abs(norm - 1.0) > 100 * Float64_threshold()
        @warn string("WARNING !! `shannon_entropy()` :: Input vector does not have unit norm !!",
                     "Normalizing the input vector !!")
        p /= norm
    end

    S::Float64 = 0.0
    for a in p
        if a > 0.0
            S += -a * log(a)
        end
    end
    return S
end

#################################################################################

"""
    function entropy(psi::MPS, bond::Int)

Compute von Neumann Entropy of a given MPS `psi` at `bond`.
"""
function entropy(psi::MPS, bond::Int)::Float64

    orthogonalize!(psi, bond)
    inds = uniqueinds(psi[bond], psi[bond+1])
    _, _, _, spec = svd(psi[bond], inds)

    eigs = spec.eigs
    if abs(sum(eigs) - 1.0) > 100 * Float64_threshold()
        @warn string("WARNING !! `entropy()` :: Input state does not have unit norm !!",
                     "Normalizing the eigenvalues !!")
        eigs /= sum(eigs)
    end
    
    S::Float64 = 0.0
    for a in eigs
        if a > 0.0
            S += -a * log(a)
        end
    end
    return S
end

#################################################################################

"""
    function entropy(psi::MPS; kwargs...)

Compute von Neumann Entropies of a given MPS `psi` at all the bonds.

Optionally, for specific bonds, keyword argument `bonds` can be specified, e.g.,
`bonds = [1, 2, 3]`.
"""
function entropy(psi::MPS; kwargs...)::Vector{Float64}
    N = length(psi)
    bonds = get(kwargs, :bonds, 1:N-1)
    return map(bond -> entropy(psi, bond), bonds)
end

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::MPS, opten::ITensor)

Returns (complex / real) local expectation value (`::ComplexF64` / `::Float64`)
of a given MPS `psi::MPS`. The operator `opten::ITensor` must be single-site operator.

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::MPS,
                 opten::ITensor)::T where T <: Union{ComplexF64, Float64}

    pos = findsite(psi, opten)
    pass = hasind(opten, siteind(psi, pos)) &&
        hasind(opten, siteind(psi, pos)') &&
        ndims(opten) == 2
    if !pass
        error("`measure`: Error in operator tensors !!")
    end
    removeqn = !hasqns(opten) && hasqns(psi) ? true : false
    
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A

    orthogonalize!(psi, pos)
    ket = conditional_removeqns(psi[pos])
    bra = dag(prime(ket; tags = "Site"))
    val = complex(scalar(bra * opten * ket))
    
    if T <: Complex
        return val
    else
        if abs(imag(val)) > 100 * Float64_threshold()
            @warn "`measure(::Float64)`: Got complex number with |imaginary part| > $(100 * Float64_threshold()) !!"
        end
        return real(val)
    end
end

measure(psi::MPS, opten::ITensor) = measure(ComplexF64, psi, opten)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::MPS, opstr::String, pos::Int)

Returns (complex / real) local expectation value (`::ComplexF64` / `::Float64`) of a
given MPS `psi::MPS` for a given operator name (`opstr::String`) and a site index (`pos::Int`) 

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
measure(::Type{T}, psi::MPS, opstr::String,
        pos::Int) where T <: Union{ComplexF64, Float64} =
            measure(T, psi, op(opstr, siteind(psi, pos)))

measure(psi::MPS, opstr::String, pos::Int) =
    measure(ComplexF64, psi, opstr, Int)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::MPS, opstr::String; kwargs...)

Returns (complex / real) local expectation values (`::Vector{ComplexF64}` / `::Vector{Float64}`)
of a given MPS `psi::MPS` for a given operator name (`opstr::String`) at all the sites.

Optionally, for specific sites, keyword argument `sites` can be specified, e.g.,
`sites = [1, 2, 3]`.

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::MPS, opstr::String;
                 kwargs...) where T <: Union{ComplexF64, Float64}
    N = length(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> measure(T, psi, opstr, pos), sites)
end

measure(psi::MPS, opstr::String; kwargs...) =
    measure(ComplexF64, psi, opstr; kwargs)

#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::MPS, optens::Vector{ITensor})
    
Returns (complex / real) multi-site expectation value (`::ComplexF64` / `::Float64`)
of a given MPS `psi::MPS` for a given vector of single-site operators
(`optens::Vector{ITensor}`).

For `ElT = Float64`, if the expectation value is complex, raises a warning!
"""
function measure(::Type{T}, psi::MPS,
                 optens::Vector{ITensor})::T where T <: Union{ComplexF64, Float64}

    oppairs = Vector{Pair{ITensor, Int}}()
    removeqn = false
    for o in optens
        pos = findsite(psi, o) 
        pass = hasind(o, siteind(psi, pos)) &&
            hasind(o, siteind(psi, pos)') &&
            ndims(o) == 2
        if !pass
            error("`measure()`: Error in operator tensors !!")
        end
        !hasqns(o) && hasqns(psi) && (removeqn = true)
        push!(oppairs, o => pos)
    end
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A
    
    #************************************************************
    
    opdict = Dict{Int, ITensor}()
    for o in oppairs
        if !haskey(opdict, o.second)
            opdict[o.second] = conditional_removeqns(o.first)
        else
            opdict[o.second] = prime(opdict[o.second])
            opdict[o.second] *= conditional_removeqns(o.first)
            opdict[o.second] = mapprime(opdict[o.second], 2, 1)
        end
    end
    
    #************************************************************
    
    if length(opdict) == 1
        fkey = first(keys(opdict)) 
        return expectC(psi, opdict[fkey])
    end

    #************************************************************

    minpos = minimum(keys(opdict))
    maxpos = maximum(keys(opdict))
    
    #************************************************************

    orthogonalize!(psi, minpos)
    ir = commonind(psi[minpos], psi[minpos+1]; tags = "Link")
    C = conditional_removeqns(psi[minpos])
    C *= conditional_removeqns(opdict[minpos])
    C *= conditional_removeqns(dag(prime(prime(psi[minpos], "Site"), ir)))
    for pos = minpos+1 : maxpos-1
        C *= conditional_removeqns(psi[pos])
        if haskey(opdict, pos)
            C *= conditional_removeqns(opdict[pos])
            C *= conditional_removeqns(dag(prime(psi[pos])))
        else
            C *= conditional_removeqns(dag(prime(psi[pos]; tags = "Link")))
        end
    end
    jl = commonind(psi[maxpos], psi[maxpos-1]; tags="Link")
    C *= conditional_removeqns(psi[maxpos])
    C *= conditional_removeqns(opdict[maxpos])
    C *= conditional_removeqns(dag(prime(prime(psi[maxpos], jl), "Site")))

    val = complex(scalar(C))
    
    if T <: Complex
        return val
    else
        if abs(imag(val)) > 100 * Float64_threshold()
            @warn "`measure(::Float64)`: Got complex number with |imaginary part| > $(100 * Float64_threshold()) !!"
        end
        return real(val)
    end
end

measure(psi::MPS, optens::Vector{ITensor}) =
    measure(ComplexF64, psi, optens)
    
#################################################################################

"""
    function measure(::Type{ElT} = ComplexF64, psi::MPS, oppairs::Vector{Pair{String, Int}};
                     isfermions::Bool = true)

Returns (complex / real) multi-site expectation value (`::ComplexF64` / `::Float64`) of
a given MPS `psi::MPS`. `oppairs::Vector{Pair{String, Int}}` contains pairs of operator
names (`String`) and site locations (`Int`).
E.g., for <ψ|OᵢOⱼOₖ... |ψ>, `oppairs = ["O" => i, "O" => j, "O" => k,...]`.

Fermionic JW strings are added automatically for fermionic operators if `isfermions::Bool = true` (default).

For `ElT = Float64`, if the expectation value is complex, raises a warning!

#### Example:
    valueC = measure(psi, ["Cdag" => 2, "C" => 6, "Cdag" => 9, "C" => 12])
    valueR = measure(Float64, psi, ["Cdag" => 2, "C" => 6, "Cdag" => 9, "C" => 12])
"""
function measure(::Type{T}, psi::MPS, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::T where T <: Union{ComplexF64, Float64}
    
    sites = siteinds(psi)
    if isfermions
        perm, oppairs = bosonize(oppairs, sites)
    end
    
    optens = ITensor[op(x.first, siteind(psi, x.second)) for x in oppairs]
    return isfermions ? perm * measure(T, psi, optens) : measure(T, psi, optens) 
end

measure(psi::MPS, oppairs::Vector{Pair{String, Int}}; isfermions::Bool = true) =
    measure(ComplexF64, psi, oppairs; isfermions)

#################################################################################



