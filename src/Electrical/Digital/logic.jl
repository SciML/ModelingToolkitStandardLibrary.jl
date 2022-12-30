struct Logic
    level::Int
end

function Base.show(io::IO, ::MIME"text/plain", l::Logic)
    if l.level == 1
        print(io, "U")
    elseif l.level == 2
        print(io, "X")
    elseif l.level == 3
        print(io, "F0")
    elseif l.level == 4
        print(io, "F1")
    elseif l.level == 5
        print(io, "Z")
    elseif l.level == 6
        print(io, "W")
    elseif l.level == 7
        print(io, "L")
    elseif l.level == 8
        print(io, "H")
    elseif l.level == 9
        print(io, "DC")
    else
        print(io, "Invalid logic level")
    end
end

# Aliases for 9 Digital Logic levels
const Uninitialized = Logic(1)
const ForcingUnknown = Logic(2)
const ForcingZero = Logic(3)
const ForcingOne = Logic(4)
const HighImpedence = Logic(5)
const WeakUnknown = Logic(6)
const WeakZero = Logic(7)
const WeakOne = Logic(8)
const DontCare = Logic(9)

const U = Uninitialized
const X = ForcingUnknown
const F0 = ForcingZero
const F1 = ForcingOne
const Z = HighImpedence
const W = WeakUnknown
const L = WeakZero
const H = WeakOne
const DC = DontCare

Base.zero(l::Logic) = F0
Base.zero(l::Type{Logic}) = F0
Base.one(l::Logic) = F1
Base.one(l::Type{Logic}) = F1

# Helpers to convert 1 and 0 to their `Logic` counterparts
function Base.convert(l::Type{Logic}, i::Number)
    if i == zero(i)
        zero(l)
    elseif i == one(i)
        one(l)
    else
        throw("This isn't a valid `Logic` value")
    end
end
Base.convert(l::Type{Logic}, i::Logic) = i

convert_to_logic(i::Number) = convert(Logic, i)
convert_to_logic(i::Logic) = i
