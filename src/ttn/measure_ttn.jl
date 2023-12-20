
#################################################################################

function expectC(psi::TTN, opten::ITensor)::ComplexF64
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
    val = scalar(bra * opten * ket)
    return complex(val)
end
    
#################################################################################

function expectR(psi::TTN, opten::ITensor)::Float64
    val = expectC(psi, opten)
    if abs(imag(val)) > 100 * Float64_threshold()
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################

expectC(psi::TTN, opstr::String, pos::Int) = expectC(psi, op(opstr, siteind(psi, pos)))

#################################################################################

expectR(psi::TTN, opstr::String, pos::Int) = expectR(psi, op(opstr, siteind(psi, pos)))

#################################################################################

function expectC(psi::TTN, opstr::String; kwargs...)
    N = numsites(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> expectC(psi, opstr, pos), sites)
end

#################################################################################

function expectR(psi::TTN, opstr::String; kwargs...)
    N = numsites(psi)
    sites = get(kwargs, :sites, 1:N)
    return map(pos -> expectR(psi, opstr, pos), sites)
end

#################################################################################

function expectC(psi::TTN, optens::Vector{ITensor})::ComplexF64

    sitenodes = Set{Int2}()
    for o in optens
        push!(sitenodes, findnode(psi, o))
    end

    oc = find_eccentric_central_node(psi.graph, sitenodes)
    isometrize!(psi, oc)

    lnt = LinkTensorsTTN(psi, optens)
    phi = psi[oc]
    return complex(scalar(product(lnt, psi, phi) * dag(phi)))
end
        
#################################################################################

function expectR(psi::TTN, optens::Vector{ITensor})::Float64
    val = expectC(psi, optens)
    if abs(imag(val)) > 100 * Float64_threshold()
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################

function expectC(psi::TTN, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::ComplexF64
    
    sites = siteinds(psi)
    if isfermions
        perm, oppairs = bosonize(oppairs, sites)
    end
    
    optens = ITensor[op(x.first, siteind(psi, x.second)) for x in oppairs]
    return isfermions ? perm * expectC(psi, optens) : expectC(psi, optens) 
end

#################################################################################

function expectR(psi::TTN, oppairs::Vector{Pair{String, Int}};
                 isfermions::Bool = true)::Float64
    val = expectC(psi, oppairs; isfermions)
    if abs(imag(val)) > 100 * _Float64_Threshold
        @warn "`expectR()`: Got complex number with |imaginary part| > $(100 * Float64_threshold()), use `expectC instead !!"
    end
    return real(val)
end

#################################################################################
