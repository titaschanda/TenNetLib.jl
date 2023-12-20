
#################################################################################

function _distribute_site_positions(numsites, loc = Int[])::Vector{Int}
    
    if length(loc) == 0
        rem = numsites % 2
        push!(loc, numsites÷2)
        push!(loc, numsites÷2 + rem)
        numsites = _minimum_power2_greater_than(numsites)
    else
        newloc = Int[]
        power = trailing_zeros(length(loc))
        for m in loc
            rem = m % 2
            push!(newloc, power % 2 == 1 ? m÷2 + rem : m÷2)
            push!(newloc, power % 2 == 1 ? m÷2 : m÷2 + rem)
        end
        loc = newloc
    end
    
    numsites ÷= 2
    power = trailing_zeros(numsites)
    if power > 1
        return _distribute_site_positions(numsites, loc)
    else
        pos = Int[]
        for m in loc
            if m == 1
                push!(pos, 1)
                push!(pos, 0)
            elseif m==2
                push!(pos, 1)
                push!(pos, 1)
            else
                error("`_distribute_site_positions()`: SOMETHING IS WRONG !!")
            end
        end
        return pos
    end
end

#################################################################################

function _default_graph_from_numsites(numsites::Int)::Tuple{Graph{Int2}, Dict{Int, Int2}}

    site_positions = _distribute_site_positions(numsites)
    pow2numsites = length(site_positions)
    numlayers = trailing_zeros(pow2numsites)

    sitelocs = Dict{Int, Int2}()
    graph = Graph{Int2}()
    
    sitecount = 1
    for ll = 1 : numlayers-2, nn = 1 : pow2numsites >> ll        
        if ll == 1 && (site_positions[2*nn-1] == 0 || site_positions[2*nn] == 0) 
            sitelocs[sitecount] = (ll+1, (nn+1)÷2)
            sitecount += 1
            continue
        elseif ll == 1
            sitelocs[sitecount] = (ll, nn)
            sitecount += 1
            sitelocs[sitecount] = (ll, nn)
            sitecount += 1
        end                
        addedge!(graph, (ll, nn), (ll+1, (nn+1)÷2))    
    end
    addedge!(graph, (numlayers-1, 1), (numlayers-1, 2))
    return graph, sitelocs
end

#################################################################################

function _ttn_ind_reducedim!(graph::Graph{Int2}, inds::Dict{Int2, Vector{Index{QNBlocks}}},
                             maxdim::Int)

    for node in graph.nodes
        indlocs = findall(x -> dim(x) > maxdim && hastags(x, "Link"), inds[node])

        for indloc in indlocs
            ind = inds[node][indloc]
            
            othernode = only(filter(x -> ind in inds[x], graph[node]))
            otherindloc = findfirst(x -> x == ind, inds[othernode])
        
            qnblocks = space(ind)
            blockdims::Vector{Int} = Int[x.second for x in qnblocks]
        
            if sum(blockdims) > maxdim
                blockdims = Int.(round.(maxdim * (blockdims / sum(blockdims))))
                for ii = 1 : length(qnblocks)           
                    qnblocks[ii] = qnblocks[ii].first => blockdims[ii]
                end
                qnblocks = filter(x -> x.second > 0, qnblocks)
            end
        
            inds[node][indloc] = Index(qnblocks; dir = dir(ind), tags = tags(ind))
            inds[othernode][otherindloc] = dag(inds[node][indloc])

        end
    end
end

#################################################################################

function _ttn_ind_cleanup_one!(graph::Graph{Int2}, inds::Dict{Int2, Vector{Index{QNBlocks}}},
                               node::Int2, fixedind::Index{QNBlocks})

    indlocs = findall(x -> x != fixedind && hastags(x, "Link"), inds[node])
    for loc in indlocs
        indold = inds[node][loc] 
        inds[node][loc] = indexintersection([indold],
                                            dag.(filter(x -> x != indold, inds[node]));
                                            dir = dir(indold),
                                            tags = tags(indold))
        
        othernode = only(filter(x -> indold in inds[x], graph[node]))
        otherloc = findfirst(x -> x == indold, inds[othernode])
        inds[othernode][otherloc] = dag(inds[node][loc])
    end
end

#################################################################################

function _ttn_ind_cleanup!(graph::Graph{Int2}, inds::Dict{Int2, Vector{Index{QNBlocks}}},
                           center_node::Int2, firstind::Index{QNBlocks})

    @assert firstind in inds[center_node]
    bfs_path = nodes_from_bfs(graph, center_node; reverse=false)
    
    _ttn_ind_cleanup_one!(graph, inds, center_node, firstind)
    for node in bfs_path[2:end]        
        prevnode = nextnode_in_path(graph, node, center_node)
        
        fixedind = commonind(inds[node], inds[prevnode])
        _ttn_ind_cleanup_one!(graph, inds, node, fixedind)
    end
end


#################################################################################

# TO DO : ADD TTN from graph and sitelocs

#################################################################################

"""
    function TTN(sites::Vector{Index{T}}, chi::Int, 
                 qn::QN = QN()) where T

Returns a TTN object from `sites` with initial bond dimension `chi` and global QN sector `qn`.
The TTN structure is a default binary graph. Automatically handles situations where the number
of sites is not a power of 2.
"""
function TTN(sites::Vector{Index{T}}, chi::Int, 
             qn::QN = QN()) where T
    
    numsites = length(sites)
    for b = 1 : numsites
        if !hastags(sites[b], "Site")
            error("""`TTN()`: Input `sites` must have `tag="Site" !!""")
        end
        if !hastags(sites[b], "n=$b")
            error("""`TTN()`: Input `sites` must have `tag="n=$b" !!""")
        end 
    end

    graph, sitelocs = _default_graph_from_numsites(numsites)

    inds = Dict{Int2, Vector{Index{T}}}()    
    for node in graph.nodes        
        inds[node] = Index{T}[]
    end
    for b = 1 : numsites
        node = sitelocs[b]
        push!(inds[node], sites[b])
    end

    center_node = find_eccentric_central_node(graph, collect(values(sitelocs)))
    bfs_path = nodes_from_bfs(graph, center_node; reverse=true)

    for node in bfs_path[begin:end-1]
        nextnode = nextnode_in_path(graph, node, center_node)
        # for QN give initial chi 2^30
        nextind = combineinds(inds[node];
                              maxdim = hasqns(sites) ? (1 << 30) : chi,
                              tags = "Link,$(node)")
        push!(inds[node], dag(nextind))
        push!(inds[nextnode], nextind)
    end
    
    if hasqns(sites)
        qnindex = Index(qn => 1; tags = "QN", dir = ITensors.In)
        push!(inds[center_node], qnindex)        
        _ttn_ind_cleanup!(graph, inds, center_node, qnindex)
        _ttn_ind_reducedim!(graph, inds, chi)
    end
    
    tensors = Dict{Int2, ITensor}()
    for node in graph.nodes
        tensors[node] = randomITensor(inds[node])
        normalize!(tensors[node])
    end

    ttn = TTN(sites, graph, tensors, Int2())
    isometrize_full!(ttn, center_node; normalize=true, cutoff=0.0)
    return ttn
end

#################################################################################
