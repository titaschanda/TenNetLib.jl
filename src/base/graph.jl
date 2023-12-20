
#################################################################################

"""
    mutable struct Graph{T}
        nodes::Set{T}
        edges::Dict{T, Set{T}}
    end

Undirected graph where nodes are of type `T`.
 - `nodes::Set{T}`: Holds the nodes of the graph.
 - `edges::Dict{T, Set{T}}`: Holds the edges of the graph. For a given node as `key`,
    holds all the connected nodes as `value` of the `Dict`.
"""
mutable struct Graph{T}
    nodes::Set{T}
    edges::Dict{T, Set{T}}
end

#################################################################################

"""
    Graph{T}()

Contructor of empty `Graph`.
"""
Graph{T}() where T = Graph{T}(Set{T}(), Dict{T, Set{T}}())

#################################################################################

"""
    Graph{T}(nodes::Set{T}) where T
    Graph{T}(nodes::Vector{T}) where T

Constructor of `Graph`.
 - `nodes`: Intial nodes.
"""
Graph{T}(nodes::Set{T}) where T = 
    Graph{T}(nodes, Dict(key => Set{T}() for key in nodes))

Graph{T}(nodes::Vector{T}) where T = Graph{T}(Set(nodes))

#################################################################################

"""
    Base.copy(graph::Graph{T}) where T

Shallow copy of `Graph`.
"""
function Base.copy(graph::Graph{T}) where T
    newedges = Dict{T, Set{T}}()
    for (key, value) in graph.edges
        newedges[key] = Base.copy(value)
    end    
    return Graph{T}(Base.copy(graph.nodes), newedges)
end

#################################################################################

"""
    Base.getindex(graph::Graph{T}, node::T) where T

Returns the `edges` for a given `node::T` as an index. 
"""
Base.getindex(graph::Graph{T}, node::T) where T =
    Base.getindex(graph.edges, node)

#################################################################################

"""
    Base.setindex!(graph::Graph{T}, neighbors::Set{T}, node::T) where T
    Base.setindex!(graph::Graph{T}, neighbors::Vector{T}, node::T) where T

Sets the `edges=neighbors` for a given `node::T` as an index.
"""
function Base.setindex!(graph::Graph{T}, neighbors::Set{T}, node::T) where T 
    Base.setindex!(graph.edges, neighbors, node)
    push!(graph.nodes, node)
    return graph
end

Base.setindex!(graph::Graph{T}, neighbors::Vector{T}, node::T) where T =
    Base.setindex!(graph, Set(neighbors), node)

#################################################################################

"""
    Base.:(==)(graph1::Graph{T}, graph2::Graph{T}) where T

Equality between two `Graph` objects.
"""
Base.:(==)(graph1::Graph{T}, graph2::Graph{T}) where T =
    graph1.nodes == graph2.nodes &&
    graph1.edges == graph2.edges

#################################################################################

"""
    getnodes(graph::Graph{T}) where T = copy(graph.nodes)

Returns (shallow copy of) nodes in the `Graph`.
"""
getnodes(graph::Graph{T}) where T = copy(graph.nodes)

#################################################################################

"""
    hasnode(graph::Graph{T}, node::T) where T

Checks whether a node is in the `Graph`.
"""
hasnode(graph::Graph{T}, node::T) where T = node in graph.nodes

#################################################################################

"""
    addnode!(graph::Graph{T}, node::T)::Graph{T} where T

Adds a node to the `Graph`.
"""
function addnode!(graph::Graph{T}, node::T)::Nothing where T
    if !hasnode(graph, node)
        push!(graph.nodes, node)
        graph[node] = Set{T}()
    end
    return nothing
end

#################################################################################

"""
    addedge!(graph::Graph{T}, node1::T, node2::T)::Graph{T} where T

Adds an edge between `node1` and `node2`. If `node1` or `node2` are not present in
the `Graph`, they are added.
"""
function addedge!(graph::Graph{T}, node1::T, node2::T)::Nothing where T    
    if hasnode(graph, node1)
        node1 != node2 &&  push!(graph[node1], node2)
    else
        graph[node1] = Set(T[node2])
    end  
    if hasnode(graph, node2)
        node1 != node2 && push!(graph[node2], node1)
    else
        graph[node2] = Set(T[node1])
    end    
    return nothing
end

#################################################################################

"""
    isneighbor(graph::Graph{T}, node1::T, node2::T)::Bool where T

Checks whether `node1` and `node2` are connected by an edge.
"""
function isneighbor(graph::Graph{T}, node1::T, node2::T)::Bool where T
    if !hasnode(graph, node1)
        error("`isneighbor()` : Input `node=$node1 does not exist in the `Graph` !!")
    end
    if !hasnode(graph, node2)
        error("`isneighbor()` : Input `node=$node2 does not exist in the `Graph` !!")
    end    
    return node2 in graph[node1] && node1 in graph[node2]
end

#################################################################################

"""
    bfs(graph::Graph{T},
        source::T,
        destination::Union{Nothing, T} = nothing
        )::Tuple{Dict{T, Int}, Dict{T, T}} where T

Performs a full BFS starting from the node `source` (and optionally, to the `destination`).

#### Return values
 - `::Dict{T, Int}`: Distances of nodes (key) in the BFS path.
 - `::Dict{T, T}`: Parents of nodes (key) in the BFS path.
"""
function bfs(graph::Graph{T},
             source::T,
             destination::Union{Nothing, T} = nothing
             )::Tuple{Dict{T, Int}, Dict{T, T}} where T 

    q = Deque{T}()
    push!(q, source)    
    #visited = Dict{T, Bool}()
    visited = Set{T}()
    parents = Dict{T, T}()
    distances = Dict{T, Int}()
    distances[source] = 0

    #visited[source] = true     
    push!(visited, source)
    
    while !isempty(q)
        node = first(q)
        popfirst!(q)        
        for neighbor in graph[node]
            #if !haskey(visited, neighbor) || !visited[neighbor]
            if !(neighbor in visited)
                push!(q, neighbor)
                #visited[neighbor] = true
                push!(visited, neighbor)
                distances[neighbor] = distances[node] + 1
                parents[neighbor] = node
                
                if !isnothing(destination) && neighbor == destination
                    return distances, parents
                end
            end
        end
    end
    return distances, parents
end

#################################################################################

"""
    nodes_from_bfs(graph::Graph{T}, source::T;
                   reverse::Bool = false)::Vector{T} where T

    nodes_from_bfs(graph::Graph{T}, source::T,
                   destinations::Union{Set{T}, Vector{T}};
                   reverse::Bool = false)::Vector{T} where T

Returns the `Vector` of nodes in the BFS path from the `source` (optionally, towards the
`destinations`). If `reverse=true` returns the reverse order of nodes.
"""
function nodes_from_bfs(graph::Graph{T}, source::T;
                        reverse::Bool = false)::Vector{T} where T
    distances, _ = bfs(graph, source)
    nodes = sort(collect(keys(distances)); by = x -> distances[x], rev = reverse)
    return nodes
end


function nodes_from_bfs(graph::Graph{T}, source::T,
                        destinations::Union{Set{T}, Vector{T}};
                        reverse::Bool = false)::Vector{T} where T

    distances, parents = bfs(graph, source)
    
    for dest_node in destinations
        if  dest_node !=source && !(dest_node in keys(parents))
            error(string("`nodes_from_bfs()`: `destination=$dest_node` is not reachable ",
                         "from `source=$source` !!"))
        end
    end
    
    result_nodes = Set{T}()
    for dest_node in destinations
        if dest_node != source
            current_node = dest_node
            push!(result_nodes, current_node)
            while haskey(parents, current_node)
                push!(result_nodes, parents[current_node])
                current_node = parents[current_node]
            end
        else
            push!(result_nodes, source)
        end
    end
    
    nodevec = collect(result_nodes)
    sort!(nodevec; by = x -> distances[x], rev = reverse)
    return nodevec
end

#################################################################################

"""
    shortest_path(graph::Graph{T}, source::T, destination::T)::Vector{T} where T

Finds the shortest path between `source` and `destination` in a `Graph`.
"""
function shortest_path(graph::Graph{T}, source::T, destination::T)::Vector{T} where T
    if source == destination
        return [source]
    end
    
    _, parents = bfs(graph, source, destination)
    
    if  destination !=source && !(destination in keys(parents))
        error(string("`shortest_path()`: `destination=$destination` is not reachable ",
                     "from `source=$source` !!"))
    end
    
    nodelist = Vector{T}()
    push!(nodelist, destination)
    while haskey(parents, destination)
        push!(nodelist, parents[destination])
        destination = parents[destination]
    end
    return reverse(nodelist)
end


#################################################################################

"""
    nextnode_in_path(graph::Graph{T}, source::T, destination::T, n=1) where T

Finds the next node in the shortest path between `source` and `destination` in a `Graph`.
Optionally, finds `n`th node in the path.
"""
nextnode_in_path(graph::Graph{T}, source::T, destination::T, n=1) where T = 
    shortest_path(graph, source, destination)[n+1]

#################################################################################

"""
    _has_cycle_dfs(graph::Graph{T}, node::T, parents::Union{Nothing, T}, 
                   visited::Dict{T, Bool})::Bool where T

Internal helper function to check whether a `Graph` has cycles/loops in it.
"""
function _has_cycle_dfs(graph::Graph{T}, node::T, parents::Union{Nothing, T}, 
                        visited::Dict{T, Bool})::Bool where T
    visited[node] = true
    for neighbor in graph[node]
        if (!haskey(visited, neighbor) || !visited[neighbor]) && (parents !== neighbor)
            if _has_cycle_dfs(graph, neighbor, node, visited)
                return true
            end
        elseif visited[neighbor] && neighbor !== parents
            return true
        end
    end
    return false
end

#################################################################################

"""
    has_cycle(graph::Graph{T}) where T

Checks whether a `Graph` has cycles/loops in it.
"""
function has_cycle(graph::Graph{T}) where T
    visited = Dict{T, Bool}()
    for node in graph.nodes
        if !haskey(visited, node) || !visited[node]
            if _has_cycle_dfs(graph, node, nothing, visited)
                return true
            end
        end
    end
    return false
end

#################################################################################

"""
    find_sum_central_node(graph::Graph{T}) where T
    find_sum_central_node(graph::Graph{T}, nodes::Set{T}) where T
    find_sum_central_node(graph::Graph{T}, nodes::Vector{T}) where T

Finds the center node of a `Graph`. Optionally, when `nodes::Vector{T}` is specified, finds
the center node with respect to the `nodes`.
"""
function find_sum_central_node(graph::Graph{T},
                               nodes::Union{Nothing, Set{T}} = nothing) where T

    if isempty(graph.nodes) || has_cycle(graph)
        error("`find_sum_central_node()`: Graph is empty or contains cycles !!")
    end
    
    center_node::T = T()
    min_sum_distances::Int = typemax(Int)
    
    for node in graph.nodes
        distances, _ = bfs(graph, node)        
        !isnothing(nodes) && (distances = filter(x -> x.first in nodes, distances))
        sum_distances = sum(values(distances))
        if sum_distances < min_sum_distances
            center_node = node
            min_sum_distances = sum_distances
        end
    end 
    return center_node
end

find_sum_central_node(graph::Graph{T}, nodes::Vector{T}) where T =
    find_sum_central_node(graph, Set(nodes))

#################################################################################

"""
    find_eccentric_central_node(graph::Graph{T}) where T
    find_eccentric_central_node(graph::Graph{T}, nodes::Set{T}) where T
    find_eccentric_central_node(graph::Graph{T}, nodes::Vector{T}) where T

Finds the center node of a `Graph`. Optionally, when `nodes::Vector{T}` is specified, finds
the center node with respect to the `nodes`.
"""
function find_eccentric_central_node(graph::Graph{T},
                                     nodes::Union{Nothing, Set{T}} = nothing) where T

    if isempty(graph.nodes) || has_cycle(graph)
        error("`find_eccentric_central_node()`: Graph is empty or contains cycles !!")
    end
    
    center_node::T = T()
    max_distances::Int = typemax(Int)
    
    for node in graph.nodes
        distances, _ = bfs(graph, node)
        !isnothing(nodes) && (distances = filter(x -> x.first in nodes, distances))
        max_dist = maximum(values(distances))
        if max_dist < max_distances
            center_node = node
            max_distances = max_dist
        end
    end 
    return center_node
end

find_eccentric_central_node(graph::Graph{T}, nodes::Vector{T}) where T =
    find_eccentric_central_node(graph, Set(nodes))

#################################################################################
