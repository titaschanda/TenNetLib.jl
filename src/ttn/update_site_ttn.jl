
#################################################################################

"""
    function update_position!(sysenv::StateEnvsTTN, solver, node::Int2;
                              time_step::Union{Float64, ComplexF64, Nothing},
                              normalize::Bool,
                              maxdim::Int,
                              mindim::Int,
                              cutoff::Float64,
                              svd_alg::String,
                              kwargs...)

Moves the orthogonality center to `pos` and update StateEnvsTTN at position `pos` by `solver`.

#### Arguments:
 - `sysenv::StateEnvsTTN`
 - `solver`: Solver for update. Currently only `eig_solver` is supported.
 - `pos::Int`: Position of the node to be updated.

#### Named arguments and their default values:
 - `time_step::Union{Float64, ComplexF64, Nothing} = nothing`: Time step for future functionality.
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`.
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 5`.
 - `solver_maxiter::Int = 2`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `::Union{Float64, ComplexF64}`: Energy. It is complex if `ishermitian == false`.
"""
function update_position!(sysenv::StateEnvsTTN, solver, node::Int2;
                          time_step::Union{Float64, ComplexF64, Nothing},
                          normalize::Bool,
                          maxdim::Int,
                          mindim::Int,
                          cutoff::Float64,
                          svd_alg::String,
                          kwargs...)::Union{Float64, ComplexF64}
    
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

"""
    function subspace_expand!(psi::TTN, node::Int2, nextnode::Int2,
                              max_expand_dim::Int, noise::Float64)

Enlarges the bond domension between `node` and `nextnode` by `max_expand_dim` using the
subspace expansion.
The parameter `noise` controls stength of the perturbation.
"""
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
    padding_tensor_B = random_itensor(ind_padB, indsB)
    padding_tensor_B *= noise * norm(B) / norm(padding_tensor_B)
    enlargedB, sumindB = directsum(B => ind_to_update,
                                   padding_tensor_B => ind_padB;
                                   tags = tags(ind_to_update))
    psi.tensors[nextnode] = enlargedB
    ##

    ## padding A
    padding_tensor_A = random_itensor(ind_padA, indsA)
    padding_tensor_A *= noise * norm(A) / norm(padding_tensor_A)
    enlargedA, sumindA = directsum(A => dag(ind_to_update),
                                   padding_tensor_A => ind_padA)
    replaceind!(enlargedA, sumindA, dag(sumindB))
    psi.tensors[node] = enlargedA
    ##

    return nothing
end

#################################################################################
