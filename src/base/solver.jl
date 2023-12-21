
#################################################################################
"""
    function eig_solver(env, phi0::ITensor, time_step::Nothing; kwargs...)

Solver to find smallest eigenvalue corresponding to "matrix" `env` and
input vector `phi0`.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`
 - `solver_tol::Float64 = 1E-14`.
 - `solver_krylovdim::Int = 5`.
 - `solver_maxiter::Int = 2`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = false`.
 - `solver_check_convergence::Bool = false`.

#### Return values:
 - `::Float64`: Smallest eigenvalue.
 - `::ITensor`: Eigenstate corresponding to the smallest eigenvalue.
"""
function eig_solver(env, phi0::ITensor,
                    time_step::Nothing; 
                    kwargs...)::Tuple{Float64, ITensor}
    howmany = 1
    which = get(kwargs, :solver_which_eigenvalue, :SR)
    solver_kwargs = (;
                     ishermitian = get(kwargs, :ishermitian, true),
                     tol = get(kwargs, :solver_tol, 1E-14),
                     krylovdim = get(kwargs, :solver_krylovdim, 5),
                     maxiter = get(kwargs, :solver_maxiter, 2),
                     verbosity = get(kwargs, :solver_outputlevel, 0),
                     eager = get(kwargs, :solver_eager, false)
                     )
    vals, vecs, info = eigsolve(env, phi0, howmany, which; solver_kwargs...)

    check_conv = get(kwargs, :solver_check_convergence, false)
    if check_conv && info.converged < 1
        error("`eig_solver()` not converged !!")
    end

    return vals[1], vecs[1]
end

#################################################################################

"""
    function exp_solver(env, phi0::ITensor, time_step::Union{Float64, ComplexF64}; kwargs...)

Exponentiation solver to find `exp(env * phi0 * time_step)`.

#### Named arguments for `solver` and their default values:
See the documentation of KrylovKit.jl.
 - `ishermitian::Bool = true`
 - `solver_tol::Float64 = 1E-12`.
 - `solver_krylovdim::Int = 30`.
 - `solver_maxiter::Int = 10`.
 - `solver_outputlevel::Int = 0`: See `verbosity` in KrylovKit.jl.
 - `solver_eager::Bool = true`.
 - `solver_check_convergence::Bool = true`.

#### Return values:
 - `::Float64`: `NaN`. 
 - `::ITensor`: Exponentiated `ITensor`.
"""
function exp_solver(env, phi0::ITensor,
                    time_step::Union{Float64, ComplexF64}; 
                    kwargs...)::Tuple{Float64, ITensor}
    if isnothing(time_step)
        error("`exp_solver()` is not defined with `time_step=$time_step` !!")
    end
    solver_kwargs = (;
                     ishermitian = get(kwargs, :ishermitian, true),
                     tol = get(kwargs, :solver_tol, 1E-12),
                     krylovdim = get(kwargs, :solver_krylovdim, 30),
                     maxiter = get(kwargs, :solver_maxiter, 100),
                     verbosity = get(kwargs, :solver_outputlevel, 0),
                     eager = get(kwargs, :solver_eager, true),
                     )
    psi, info = exponentiate(env, time_step, phi0; solver_kwargs...)

    check_conv = get(kwargs, :solver_check_convergence, false)
    if check_conv && info.converged < 1
        error("`eig_solver()` not converged !!")
    end

    return NaN, psi
end  

#################################################################################
