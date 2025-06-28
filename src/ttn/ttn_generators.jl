
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

"""
    function default_graph_sitenodes(N::Int)

Given the total number of sites, `N::Int`, generates the default hierarchical
binary tree graph and a `Dict{Int, Int2}` object that maps each site to the corresponding node.
Automatically handles situations where the number of sites is not a power of 2.

#### Return values:
 - `::Graph{Int2}`: Default hierarchical tree graph to accomodate `N` number of sites
   in a TTN.
 - `::Dict{Int, Int2}`: Maps each site to the corresponding node.

#### Example:
```
graph, sitenodes = default_graph_sitenodes(32)

sitenodes[1] == (1,1) # true
sitenodes[2] == (1,1) # true
sitenodes[3] == (1,2) # true
sitenodes[4] == (1,2) # true
sitenodes[31] == (1,16) # true
sitenodes[32] == (1,16) # true
```
"""
function default_graph_sitenodes(numsites::Int)::Tuple{Graph{Int2}, Dict{Int, Int2}}

    site_positions = _distribute_site_positions(numsites)
    pow2numsites = length(site_positions)
    numlayers = trailing_zeros(pow2numsites)

    sitenodes = Dict{Int, Int2}()
    graph = Graph{Int2}()
    
    sitecount = 1
    for ll = 1 : numlayers-2, nn = 1 : pow2numsites >> ll        
        if ll == 1 && (site_positions[2*nn-1] == 0 || site_positions[2*nn] == 0) 
            sitenodes[sitecount] = (ll+1, (nn+1)÷2)
            sitecount += 1
            continue
        elseif ll == 1
            sitenodes[sitecount] = (ll, nn)
            sitecount += 1
            sitenodes[sitecount] = (ll, nn)
            sitecount += 1
        end                
        addedge!(graph, (ll, nn), (ll+1, (nn+1)÷2))    
    end
    addedge!(graph, (numlayers-1, 1), (numlayers-1, 2))
    return graph, sitenodes
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

"""
    function randomTTN(sites::Vector{Index{T}}, graph::Graph{Int2},
                       sitenodes::Dict{Int, Int2}, chi::Int, qn::QN = QN()) where T

Returns a TTN object, having random elements, from site `Index`s `sites`,
the underlyting `graph`, `sitenodes::Dict{Int, Int2}` that maps each site to the
corresponding node, initial bond dimension `chi`, and (optional) global QN sector `qn`.
The structure is determined by the input `graph` object.

**Note**: This function can be used to generate any loop-free tensor network.

**Note**: For QN conserving TTN, the bond dimension might be off by one or two from `chi`.
"""
function randomTTN(sites::Vector{Index{T}}, graph::Graph{Int2},
                   sitenodes::Dict{Int, Int2}, chi::Int, qn::QN = QN()) where T
    
    numsites = length(sites)
    for b = 1 : numsites
        if !hastags(sites[b], "Site")
            error("""`randomTTN()`: Input site `Index` must have `tag="Site" !!""")
        end
    end

    if length(sitenodes) != numsites
        error("""`randomTTN()`: Size mismatch between `sites` and `sitenodes` !!""")
    end

    for b = 1 : numsites
        if !(b in keys(sitenodes))
            error("""`randomTTN()`: `sitenodes` does not have site "$b" !!""")
        end
        if !(sitenodes[b] in graph.nodes)
            error("""`randomTTN()`: `sitenodes["$b"]` does not exist in `graph` !!""")
        end
    end

    if has_cycle(graph)
        error("""`randomTTN()`: `graph` has loops in it !!""")
    end
    
    inds = Dict{Int2, Vector{Index{T}}}()    
    for node in graph.nodes        
        inds[node] = Index{T}[]
    end
    for b = 1 : numsites
        node = sitenodes[b]
        push!(inds[node], sites[b])
    end

    center_node = find_eccentric_central_node(graph, collect(values(sitenodes)))
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

"""
    function default_randomTTN(sites::Vector{Index{T}}, chi::Int, qn::QN = QN()) where T

Returns a TTN object, having random elements, from site `Index`s `sites`, initial bond dimension
`chi` and (optional) global QN sector `qn`. The structure is a default hierarchical binary
tree graph. Automatically handles situations where the number of sites is not a power of 2.

**Note**: For QN conserving TTN, the bond dimension might be off by one or two from `chi`.
"""
function default_randomTTN(sites::Vector{Index{T}}, chi::Int, qn::QN = QN()) where T
    
    numsites = length(sites)
    graph, sitenodes = default_graph_sitenodes(numsites)
    return randomTTN(sites, graph, sitenodes, chi, qn)
end

#################################################################################

