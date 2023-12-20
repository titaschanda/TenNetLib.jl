
#################################################################################

"""
    shanon_entropy(p::Vector{Float64})::Float64

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
    entropy(psi::MPS, bond::Int)::Float64

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
    entropy(psi::MPS; kwargs...)::Vector{Float64}

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
    expectC(psi::MPS, opten::ITensor)::ComplexF64

Returns (complex) local expectation value (`::ComplexF64`) of a given MPS `psi::MPS`.
The operator `opten::ITensor` must be single-site operator.
"""
function expectC(psi::MPS, opten::ITensor)::ComplexF64
    pos = findsite(psi, opten)
    pass = hasind(opten, siteind(psi, pos)) &&
        hasind(opten, siteind(psi, pos)') &&
        ndims(opten) == 2
    if !pass
        error("`expectC()`: Error in operator tensors !!")
    end
    removeqn = !hasqns(opten) && hasqns(psi) ? true : false
    
    conditional_removeqns(A) = removeqn && hasqns(A) ? removeqns(A) : A

    orthogonalize!(psi, pos)
    ket = conditional_removeqns(psi[pos])
    bra = dag(prime(ket; tags = "Site"))
    val = scalar(bra * opten * ket)
    return complex(val)
end

#################################################################################

"""
    expectR(psi::MPS, opten::ITensor)::Float64

Returns (real) local expectation value (`::Float64`) of a given MPS `psi::MPS`.
The operator `opten::ITensor` must be single-site operator.

If the expectation value is complex, raises a warning!
"""
function expectR(psi::MPS, opten::ITensor)::Float64
    val = expectC(psi, opten)
    if abs(imag(val)) > 100 * Float64_threshold()
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################

"""
    expectC(psi::MPS, opstr::String, pos::Int)

Returns (complex) local expectation value (`::ComplexF64`) of a given MPS `psi::MPS` for
a given operator name (`opstr::String`) and a site number (`pos::Int`) 
"""
expectC(psi::MPS, opstr::String, pos::Int) = expectC(psi, op(opstr, siteind(psi, pos)))

#################################################################################

"""
    expectR(psi::MPS, opstr::String, pos::Int)

Returns (real) local expectation value (`::Float64`) of a given MPS `psi::MPS` for
a given operator name (`opstr::String`) and a site number (`pos::Int`).

If the expectation value is complex, raises a warning!
"""
expectR(psi::MPS, opstr::String, pos::Int) = expectR(psi, op(opstr, siteind(psi, pos)))

#################################################################################

"""
    expectC(psi::MPS, opstr::String; kwargs...)

Returns (complex) local expectation values (`::Vector{ComplexF64}`) of a given MPS `psi::MPS`
for a given operator name (`opstr::String`) at all the sites.

Optionally, for specific sites, keyword argument `sites` can be specified, e.g.,
`sites = [1, 2, 3]`.
"""
function expectC(psi::MPS, opstr::String; kwargs...)
    N = length(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> expectC(psi, opstr, pos), sites)
end

#################################################################################

"""
    expectR(psi::MPS, opstr::String; kwargs...)

Returns (real) local expectation values (`::Vector{Float64}`) of a given MPS `psi::MPS`
for a given operator name (`opstr::String`) at all the sites.

Optionally, for specific sites, keyword argument `sites` can be specified, e.g.,
`sites = [1, 2, 3]`.

If the expectation value is complex, raises a warning!
"""
function expectR(psi::MPS, opstr::String; kwargs...)
    N = length(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> expectR(psi, opstr, pos), sites)
end

#################################################################################

"""
    expectC(psi::MPS, optens::Vector{ITensor})::ComplexF64

Returns (complex) multi-site expectation value (`::ComplexF64`) of a given MPS `psi::MPS`
for a given vector of single-site operators (`optens::Vector{ITensor}`).
"""
function expectC(psi::MPS, optens::Vector{ITensor})::ComplexF64

    oppairs = Vector{Pair{ITensor, Int}}()
    removeqn = false
    for o in optens
        pos = findsite(psi, o) 
        pass = hasind(o, siteind(psi, pos)) &&
            hasind(o, siteind(psi, pos)') &&
            ndims(o) == 2
        if !pass
            error("`expectC()`: Error in operator tensors !!")
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
    
    return complex(scalar(C))
end

#################################################################################

"""
    expectR(psi::MPS, optens::Vector{ITensor})::Float64

Returns (real) multi-site expectation value (`::Float64`) of a given MPS `psi::MPS`
for a given vector of single-site operators (`optens::Vector{ITensor}`).

If the expectation value is complex, raises a warning!
"""
function expectR(psi::MPS, optens::Vector{ITensor})::Float64
    val = expectC(psi, optens)
    if abs(imag(val)) > 100 * Float64_threshold()
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################

"""
    expectC(psi::MPS, oppairs::Vector{Pair{String, Int}};
            isfermions::Bool = true)::ComplexF64

Returns (complex) multi-site expectation value (`::ComplexF64`) of a given MPS `psi::MPS`.
`oppairs::Vector{Pair{String, Int}}` contains pairs of operator names (`String`) and
site locations (`Int`). 

E.g., for <ψ|OᵢOⱼOₖ... |ψ>, `oppairs = ["O" => i, "O" => j, "O" => k,...]`.

Fermionic JW strings are added automatically for fermionic operators if `isfermions::Bool = true` (default).
"""
function expectC(psi::MPS, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::ComplexF64
    
    sites = siteinds(psi)
    if isfermions
        perm, oppairs = bosonize(oppairs, sites)
    end
    
    optens = ITensor[op(x.first, siteind(psi, x.second)) for x in oppairs]
    return isfermions ? perm * expectC(psi, optens) : expectC(psi, optens) 
end

#################################################################################

"""
    expectR(psi::MPS, oppairs::Vector{Pair{String, Int}};
            isfermions::Bool = true)::Float64

Returns (real) multi-site expectation value (`::Float64`) of a given MPS `psi::MPS`.
`oppairs::Vector{Pair{String, Int}}` contains pairs of operator names (`String`) and
site locations (`Int`). 

E.g., for <ψ|OᵢOⱼOₖ... |ψ>, `oppairs = ["O" => i, "O" => j, "O" => k,...]`.

Fermionic JW strings are added automatically for fermionic operators if `isfermions::Bool = true` (default).

If the expectation value is complex, raises a warning!
"""
function expectR(psi::MPS, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::Float64
    val = expectC(psi, oppairs; isfermions)
    if abs(imag(val)) > 100 * Float64_threshold()
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################



