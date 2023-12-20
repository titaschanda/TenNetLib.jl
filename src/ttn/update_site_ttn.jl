
#################################################################################

function update_position!(sysenv::StateEnvsTTN, solver, node::Int2;
                          time_step::Union{Float64, ComplexF64, Nothing},
                          normalize::Bool,
                          maxdim::Int,
                          mindim::Int,
                          cutoff::Float64,
                          svd_alg::String,
                          kwargs...)
    
    position!(sysenv, node;
              normalize=normalize,
              maxdim=maxdim,
              mindim=mindim,
              cutoff=cutoff,
              svd_alg=svd_alg)
    phi = sysenv.psi[node]
    
    energy, phi = solver(sysenv, phi, time_step; kwargs...)
    sysenv.psi[node] = phi
    return energy
end

#################################################################################

function subspace_expand!(psi::TTN, node::Int2, nextnode::Int2,
                          max_expand_dim::Int, noise::Float64)::Nothing

    A = psi.tensors[node]
    B = psi.tensors[nextnode]

    ind_to_update = commonind(B, A)
    indsA = uniqueinds(A, B)
    indsB = uniqueinds(B, A)
    
    ind_padB = indexintersection(indsB, dag.(indsA);
                                 maxdim = max_expand_dim,
                                 dir = dir(ind_to_update))    
    ind_padA = dag(ind_padB)

    ## padding B
    padding_tensor_B = randomITensor(ind_padB, indsB)
    padding_tensor_B *= noise * norm(B) / norm(padding_tensor_B)
    enlargedB, sumindB = directsum(B => ind_to_update,
                                   padding_tensor_B => ind_padB;
                                   tags = tags(ind_to_update))
    psi.tensors[nextnode] = enlargedB
    ##

    ## padding A
    padding_tensor_A = randomITensor(ind_padA, indsA)
    padding_tensor_A *= noise * norm(A) / norm(padding_tensor_A)
    enlargedA, sumindA = directsum(A => dag(ind_to_update),
                                   padding_tensor_A => ind_padA)
    replaceind!(enlargedA, sumindA, dag(sumindB))
    psi.tensors[node] = enlargedA
    ##

    return nothing
end

#################################################################################
