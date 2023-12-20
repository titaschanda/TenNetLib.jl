#################################################################################

"""
    const Int2 = Tuple{Int, Int}

Tuple of two `Int`. 
"""
const Int2 = Tuple{Int, Int}
Int2() = (typemin(Int), typemin(Int))

#################################################################################

"""
    struct LinkTypeTTN
        first::Int2
        second::Int2
    end

Link type for TTN.
"""
struct LinkTypeTTN
    first::Int2
    second::Int2
end

Base.hash(obj::LinkTypeTTN, h::UInt) = hash(obj.first) ⊻ hash(obj.second, h) 
Base.:(==)(obj1::LinkTypeTTN, obj2::LinkTypeTTN) =
    (obj1.first == obj2.first && obj1.second == obj2.second) ||
    (obj1.first == obj2.second && obj1.second == obj2.first)

#################################################################################
