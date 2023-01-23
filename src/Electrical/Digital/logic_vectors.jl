import Base: size, axes, getindex

const LogicOrNumber = Union{Logic, Number}

struct StdULogicVector{N} <: AbstractArray{Logic, N}
    logic::AbstractArray{Logic}
    level::AbstractArray{Int}
    function StdULogicVector(l::AbstractArray{<:LogicOrNumber})
        N = ndims(l)
        l = AbstractArray{Logic}(convert.(Logic, l))
        levels = get_logic_level.(l)
        new{N}(l, levels)
    end
end

struct StdLogicVector{N} <: AbstractArray{Logic, N}
    logic::AbstractArray{Logic}
    level::AbstractArray{Int}
    function StdLogicVector(l::AbstractArray{<:LogicOrNumber})
        N = ndims(l)
        l = AbstractArray{Logic}(convert.(Logic, l))
        levels = get_logic_level.(l)
        new{N}(l, levels)
    end
end

Base.size(l::T) where {T <: Union{StdULogicVector, StdLogicVector}} = size(l.logic)

Base.axes(l::T) where {T <: Union{StdULogicVector, StdLogicVector}} = axes(l.logic)

function Base.getindex(s::T,
                       i::Int) where {T <: Union{StdULogicVector{N}, StdLogicVector{N}}
                                      } where {N}
    getindex(s.logic, i)
end
function Base.getindex(s::T, i1::Int, i2::Int,
                       I::Int...) where {T <: Union{StdULogicVector{N}, StdLogicVector{N}}
                                         } where {N}
    getindex(s.logic, i1, i2, I...)
end

# predefined vectors

std_ulogic = StdULogicVector([U, X, F0, F1, Z, W, L, H, DC])
UX01 = StdULogicVector([U, X, F0, F1])
UX01Z = StdULogicVector([U, X, F0, F1, Z])
X01 = StdULogicVector([X, F0, F1])
X01Z = StdULogicVector([X, F0, F1, Z])
