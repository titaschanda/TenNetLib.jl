
#################################################################################

"""
    const Vector2{T} = Vector{Vector{T}}

`Vector` of `Vector`s.
"""
const Vector2{T} = Vector{Vector{T}}

################################################################################

"""
    const IDType = ITensors.IDType

Type of randomly generated id. It is usually `UInt64`.
"""
const IDType = ITensors.IDType

################################################################################

"""
    const IDTensors = Dict{IDType, ITensor}

Dictionary of `key` = randomly generated id and `value` = `ITensor` objects.
"""
const IDTensors = Dict{IDType, ITensor}

################################################################################

"""
    function Base.eltype(idt::IDTensors)

Returns the element type (e.g., `Float64` or `ComplexF64`) of the `IDTensors`.
"""
function Base.eltype(idt::IDTensors)
    if length(idt) == 0
        return Bool
    else        
        return scalartype(collect(values(idt)))
    end
end

#################################################################################

