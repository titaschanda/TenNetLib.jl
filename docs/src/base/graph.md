# The `Graph` object

TeNLib.jl defines its own undirected graph object, the `Graph`. Currently only used for the Tree
Tensor Network (TTN) codes.

```@docs
Graph
Graph{T}() where T
Graph{T}(nodes::Set{T}) where T
getnodes(graph::Graph{T}) where T
hasnode(graph::Graph{T}, node::T) where T
addnode!(graph::Graph{T}, node::T) where T
addedge!(graph::Graph{T}, node1::T, node2::T) where T
isneighbor(graph::Graph{T}, node1::T, node2::T) where T
bfs(graph::Graph{T}, source::T, destination::Union{Nothing, T} = nothing) where T
nodes_from_bfs(graph::Graph{T}, source::T;reverse::Bool = false) where T
shortest_path(graph::Graph{T}, source::T, destination::T) where T
nextnode_in_path(graph::Graph{T}, source::T, destination::T, n=1) where T
has_cycle(graph::Graph{T}) where T
find_sum_central_node(graph::Graph{T}, nodes::Union{Nothing, Set{T}} = nothing) where T
find_eccentric_central_node(graph::Graph{T}, nodes::Union{Nothing, Set{T}} = nothing) where T
Base.copy(graph::Graph{T}) where T
Base.getindex(graph::Graph{T}, node::T) where T
Base.setindex!(graph::Graph{T}, neighbors::Set{T}, node::T) where T
Base.:(==)(graph1::Graph{T}, graph2::Graph{T}) where T
```