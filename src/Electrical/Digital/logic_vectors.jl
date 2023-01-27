import Base: size, axes, getindex

const LogicOrNumber = Union{Logic, Number}

struct StdULogicVector{N} <: AbstractArray{Logic, N}
    logic::AbstractArray{Logic}
    level::AbstractArray{Int}
    function StdULogicVector(l::AbstractArray)
        N = ndims(l)
        l = AbstractArray{Logic}(convert.(Logic, l))
        new{N}(l, get_logic_level.(l))
    end
end

struct StdLogicVector{N} <: AbstractArray{Logic, N}
    logic::AbstractArray{Logic}
    level::AbstractArray{Int}
    function StdLogicVector(l::AbstractArray)
        N = ndims(l)
        l = AbstractArray{Logic}(convert.(Logic, l))
        new{N}(l, get_logic_level.(l))
    end
end

const LogicVector = Union{StdULogicVector, StdLogicVector}

Base.size(l::L) where {L <: LogicVector} = size(l.logic)

Base.axes(l::L) where {L <: LogicVector} = axes(l.logic)

function Base.getindex(s::L,
                       i::Int) where {L <: LogicVector}
    getindex(s.logic, i)
end
function Base.getindex(s::L, i1::Int, i2::Int,
                       I::Int...) where {L <: LogicVector}
    getindex(s.logic, i1, i2, I...)
end

get_logic_level(s::L) where {L <: LogicVector} = s.level

# predefined vectors
std_ulogic = StdULogicVector([U, X, F0, F1, Z, W, L, H, DC])
UX01 = StdULogicVector([U, X, F0, F1])
UX01Z = StdULogicVector([U, X, F0, F1, Z])
X01 = StdULogicVector([X, F0, F1])
X01Z = StdULogicVector([X, F0, F1, Z])
