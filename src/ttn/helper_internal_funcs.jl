
#################################################################################

"""
    _minimum_power2_greater_than(n::Int)

Internal function to find the minimum number that is a power of 2 and greater than
the input number.
"""
function _minimum_power2_greater_than(n::Int)
    if n <= 0
        Error("`_minimum_power2_greater_than()`: Input must be a positive integer")
    end    
    n -= 1
    n |= (n >>> 1)    
    while n & (n + 1) != 0
        n |= n >>> 1
    end    
    return n + 1
end

#################################################################################

function _get_links(psi::TTN, node::Int2,
                    nextnode::Int2)::Tuple{LinkTypeTTN,
                                           Vector{LinkTypeTTN}}
    next_link = LinkTypeTTN(node, nextnode)
    prev_nodes = collect(filter(x -> x != nextnode, psi.graph[node]))
    site_poses = findsites(psi, psi[node])
    if length(site_poses) > 0
        append!(prev_nodes, Int2[(0, x) for x in site_poses])
    end
    prev_links = LinkTypeTTN[LinkTypeTTN(node, x) for x in prev_nodes]
    return next_link, prev_links
end

function _get_links(psi::TTN, node::Int2)::Vector{LinkTypeTTN}    

    nodes = collect(psi.graph[node])
    site_poses = findsites(psi, psi[node])
    if length(site_poses) > 0
        append!(nodes, Int2[(0, x) for x in site_poses])
    end
    links = LinkTypeTTN[LinkTypeTTN(node, x) for x in nodes]
    return links
end

#################################################################################
