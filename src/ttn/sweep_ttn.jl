
#################################################################################

"""
    mutable struct SweepDataTTN
        sweepcount::Int
        maxchi::Vector{Int}
        energy::Vector{Float64}
    end

Holds historical data after each (full)sweep of the TTN. Requires for convergence check etc.
 - `sweepcount::Int`: Number of fullsweeps performed.
 - `maxchi::Vector{Int}`: Maximum MPS bond/link dimensions after every sweep.
 - `energy::Vector{Float64}`: Energies after every sweep.

#### Default constructor:
 - `SweepDataTTN()`: Initialize an empty `SweepData` object.
"""
mutable struct SweepDataTTN
    sweepcount::Int
    maxchi::Vector{Int}
    energy::Vector{Float64}
end

SweepDataTTN() = SweepDataTTN(0, Int[], Float64[])

#################################################################################

"""
    Base.copy(swdata::SweepDataTTN)

Shallow copy of `SweepDataTTN`.
"""
Base.copy(swdata::SweepDataTTN) = SweepDataTTN(swdata.sweepcount,
                                               Base.copy(swdata.energy),
                                               Base.copy(swdata.entropy)
                                               )
#################################################################################

"""
    function default_sweeppath(psi::TTN)

For a given a TTN having the default hierarchical binary tree structure, teturns the default
path for sweeping through the TTN.

#### Return values:
 - `Vector{Int2}`: The list of nodes as the path to be followed during a (half)sweep.
"""
function default_sweeppath(psi::TTN)
    
    sweeppath = Vector{Int2}()    
    nsites = _minimum_power2_greater_than(numsites(psi))
    nlayers = trailing_zeros(nsites)

    for ll = nlayers-1 : -1 : 1, nn = 1 : nsites >> ll
        nnpos = (nlayers - ll) % 2 == 1 ? nn : (nsites >> ll) - nn + 1
        node = (ll, nnpos)
        if ll == 1 && !(node in psi.graph.nodes)
            continue
        else
            push!(sweeppath, node)
        end
    end
    return sweeppath
end

#################################################################################

"""
    function fullsweep!(sysenv::StateEnvsTTN, sweeppath::Vector{Int2}, solver,
                        swdata::SweepDataTTN; kwargs...)

Perform a fullsweep of the TTN by `solver`.

#### Arguments:
 - `sysenv::StateEnvsTTN`.
 - `sweeppath::Vector{Int2}`: The list of nodes as the path to be followed during a (half)sweep.
   Must have each node atleast once. For the next halfsweep the reverse path is followed.
 - `solver`: Solver for update. Currently only `eig_solver` is supported.
 - `swdata::SweepDataTTN`.

#### Named arguments and their default values:
 - `time_step::Union{Float64, ComplexF64, Nothing} = nothing`: Time step for future
   functionality.
 - `normalize::Bool = true`: Whether to normalize after update.
 - `maxdim::Int = typemax(Int)`: Maximum bond dimension after SVD truncation.
 - `mindim::Int = 1`: Minimum bond dimension after SVD truncation.
 - `cutoff::Float64 = Float64_threshold()`: Cutoff for SVD truncation.
 - `svd_alg::String = "divide_and_conquer"`.
 - `noise::Float64 = 0.0`.
 - `expand_dim::Int = 0` if `noise == 0` else `20`. Dimension to be expanded (on top of `maxdim`)
   during subspace expansion.
 - `max_expand_dim::Int = 2 * expand_dim`. maximum dimension to be expanded during subspace
   expansion.
 - `expand_numiter::Int = 4`. Number of iteration for subspace expansion. Should be greater
   than 1.
 - `linkwise_maxdim::Union{Nothing, Dict{LinkTypeTTN, Int}} = nothing`. Specifies maximum bond
   dimension of a particular link.  
 - `outputlevel::Int = 1`. If `0` prints no information, for `1` outputs after
   every fullsweep, if `2` prints at every update step.

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
 - `::Float64`: Change in Energy ΔE

`swdata::SweepDataTTN` gets updated.
"""
function fullsweep!(sysenv::StateEnvsTTN, sweeppath::Vector{Int2}, solver,
                    swdata::SweepDataTTN; kwargs...)
    
    @assert issetequal(sysenv.psi.graph.nodes, Set(sweeppath))

    outputlevel::Int = get(kwargs, :outputlevel, 1)
    
    # time_step for future
    time_step::Union{Float64, ComplexF64, Nothing} = get(kwargs, :time_step, nothing)    
    
    # SVD/Decomp Algs    
    maxdim::Int = get(kwargs, :maxdim, typemax(Int))
    mindim::Int = get(kwargs, :mindim, 1)
    cutoff::Float64 = get(kwargs, :cutoff, Float64_threshold())
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")    

    normalize::Bool = get(kwargs, :normalize, true)
    noise::Float64 = get(kwargs, :noise, 0.0)
    expand_dim::Int = get(kwargs, :expand_dim,
                          abs(noise) < 100 * Float64_threshold() ?
                              0 : 20)
    max_expand_dim::Int = get(kwargs, :max_expand_dim, 2 * expand_dim)
    expand_numiter::Int = get(kwargs, :expand_numiter, 4)

    linkwise_maxdim::Union{Nothing, Dict{LinkTypeTTN, Int}} =
        get(kwargs, :linkwise_maxdim, nothing)


    if expand_dim != 0 && expand_numiter < 2
        error(string("`fullsweep!()`: `expand_numiter=$expand_numiter` cannot be less than 2",
                     " for `expand_dim=$expand_dim` !!"))
    end

    energy::Float64 = NaN
    swdata.sweepcount += 1
    sw_time::Float64 = NaN

    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end

    
    if abs(noise) < 100 * Float64_threshold()
        sw_time = @elapsed begin
            for node in [sweeppath; reverse(sweeppath)]
            #for node in sweeppath
                energy =  update_position!(sysenv, solver, node;
                                           kwargs...,
                                           time_step = time_step, 
                                           normalize = normalize,
                                           maxdim = maxdim,
                                           mindim = mindim,
                                           cutoff = -1.0,
                                           svd_alg = svd_alg)

                if (outputlevel > 1)
                    @printf("At sweep %d updated node (%d,%d) => Energy %s, MaxLinkDim %d\n",
                            swdata.sweepcount, node[1], node[2], energy,
                            ITensors.maxdim(sysenv.psi[node]))
                    flush(stdout)
                end            
            end
        end
    else
        sw_time = @elapsed begin
            central_node = find_eccentric_central_node(sysenv.psi.graph)
            for ii in [1 : length(sweeppath); length(sweeppath):-1:1]
                            
                node = sweeppath[ii]

                if node == central_node
                    energy =  update_position!(sysenv, solver, node;
                                               kwargs...,
                                               time_step = time_step, 
                                               normalize = normalize,
                                               maxdim = maxdim,
                                               mindim = mindim,
                                               cutoff = -1.0,
                                               svd_alg = svd_alg)
                    if (outputlevel > 1)
                        @printf("At sweep %d updated node (%d,%d) => Energy %s, MaxLinkDim %d\n",
                                swdata.sweepcount, node[1], node[2], energy,
                                ITensors.maxdim(sysenv.psi[node]))
                        flush(stdout)
                    end                           
                else
                    nextnode = nextnode_in_path(sysenv.psi.graph,
                                                node, central_node)
                
                    position!(sysenv, node;
                              normalize = normalize,
                              maxdim = maxdim,
                              mindim = mindim,
                              cutoff = -1.0,
                              svd_alg = svd_alg,
                              node_to_skip = nextnode)

                    subspace_expand!(sysenv.psi, node, nextnode, max_expand_dim, noise)
                    sysenv.psi.orthocenter = nextnode
            
                    link = LinkTypeTTN(node, nextnode)
                    linkmaxdim = !isnothing(linkwise_maxdim) &&
                        haskey(linkwise_maxdim, link) ? linkwise_maxdim[link] : maxdim
                
                
                    for dummy = 1 : expand_numiter
                        newmaxdim = dummy == expand_numiter ?
                            linkmaxdim : (dummy == 1 ?
                            linkmaxdim + max_expand_dim : linkmaxdim + expand_dim)
                        newnode = dummy % 2 == 1 ? node : nextnode

                        energy =  update_position!(sysenv, solver, newnode;
                                                   kwargs...,
                                                   time_step = time_step, 
                                                   normalize = normalize,
                                                   maxdim = newmaxdim,
                                                   mindim = mindim,
                                                   cutoff = cutoff,
                                                   svd_alg = svd_alg)
                    
                        if (outputlevel > 1)
                            @printf("At sweep %d updated node (%d,%d) => Energy %s, MaxLinkDim %d\n",
                                    swdata.sweepcount, newnode[1], newnode[2], energy,
                                    ITensors.maxdim(sysenv.psi[newnode]))
                            flush(stdout)
                        end
                    end
                end
            end
        end
    end
    
    push!(swdata.maxchi, maxlinkdim(sysenv.psi))
    push!(swdata.energy, energy)
    if (swdata.sweepcount > 1)
        enerr = swdata.energy[end] - swdata.energy[end-1]
    else
        enerr = NaN
    end
    
    if (outputlevel > 0)
        @printf("-----------------------------------------------------------------------------------\n")
        @printf("At sweep %d => E=%s, MaxLinkDim=%d, Noise=%0.2E\n",
                swdata.sweepcount, swdata.energy[end], swdata.maxchi[end], noise)
        @printf("At sweep %d => ΔE=%.2E, Time=%0.3f\n",
                swdata.sweepcount, enerr, sw_time)
        @printf("-----------------------------------------------------------------------------------\n")
        flush(stdout)
    end
    
    if (outputlevel > 1)
        @printf("###################################################################################\n")
        flush(stdout)
    end
    
    return enerr
end

#################################################################################
