"""
    Gain(; name, k)

Output the product of a gain value with the input signal.

# Parameters:

  - `k`: Scalar gain

# Connectors:

  - `input`
  - `output`
"""
@component function Gain(; name, k = nothing)
    @named siso = SISO()
    @unpack u, y = siso

    pars = @parameters begin
        k = k, [description = "Gain"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        y ~ k * u
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, siso)
end
Gain(k; name) = Gain(; k, name)

"""
    MatrixGain(; K::AbstractArray, name)

Output the product of a gain matrix with the input signal vector.

# Structural parameters:

  - `K`: Matrix gain

# Connectors:

  - `input`
  - `output`
"""
@component function MatrixGain(; name, K = nothing)
    nout = size(K, 1)
    nin = size(K, 2)

    pars = @parameters begin
    end

    systems = @named begin
        input = RealInput(; nin = nin)
        output = RealOutput(; nout = nout)
    end

    vars = @variables begin
    end

    equations = Equation[
        [output.u[i] ~ sum(K[i, j] * input.u[j] for j in 1:nin)
         for i in 1:nout]...
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Sum(; input__nin::Int, name)

Output the sum of the elements of the input port vector.
Input port dimension can be set with `input__nin`

# Connectors:

  - `input`
  - `output`
"""
@component function Sum(; name, input__nin = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        input = RealInput(; nin = input__nin)
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ sum(input.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Feedback(; name)

Output difference between reference input (input1) and feedback input (input2).

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@component function Feedback(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ input1.u - input2.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Add(; name, k1 = 1.0, k2 = 1.0)

Output the sum of the two scalar inputs.

# Parameters:

  - `k1`: Gain for first input
  - `k2`: Gain for second input

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@component function Add(; name, k1 = 1.0, k2 = 1.0)
    pars = @parameters begin
        k1 = k1, [description = "Gain of Add input1"]
        k2 = k2, [description = "Gain of Add input2"]
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ k1 * input1.u + k2 * input2.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Add(; name, k1 = 1.0, k2 = 1.0, k3 = 1.0)

Output the sum of the three scalar inputs.

# Parameters:

  - `k1`: Gain for first input
  - `k2`: Gain for second input
  - `k3`: Gain for third input

# Connectors:

  - `input1`
  - `input2`
  - `input3`
  - `output`
"""
@component function Add3(; name, k1 = 1.0, k2 = 1.0, k3 = 1.0)
    pars = @parameters begin
        k1 = k1, [description = "Gain of Add input1"]
        k2 = k2, [description = "Gain of Add input2"]
        k3 = k3, [description = "Gain of Add input3"]
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput()
        input3 = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ k1 * input1.u + k2 * input2.u + k3 * input3.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Product(; name)

Output product of the two inputs.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@component function Product(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ input1.u * input2.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Division(; name)

Output first input divided by second input.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@component function Division(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput(guess = 1.0) # denominator can not be zero
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ input1.u / input2.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Power(; name)

Output the exponential with base as the first input and exponent as second input i.e u1^u2

# Connectors:

  - `base`
  - `exponent`
  - `output`
"""
@component function Power(; name)
    pars = @parameters begin
    end

    systems = @named begin
        base = RealInput()
        exponent = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ base.u^exponent.u
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Modulo(; name)

Output the remainder when the first input is divided by second input.

# Connectors:

  - `dividend`
  - `divisor`
  - `remainder`
"""
@component function Modulo(; name)
    pars = @parameters begin
    end

    systems = @named begin
        dividend = RealInput()
        divisor = RealInput(guess = 1.0) # denominator can not be zero
        remainder = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        remainder.u ~ mod(dividend.u, divisor.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    UnaryMinus(; name)

Output the product of -1 and the input.

# Connectors:

  - `input`
  - `output`
"""
@component function UnaryMinus(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ -(input.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Floor(; name)

Output the floor rounding of the input.

# Connectors:

  - `input`
  - `output`
"""
@component function Floor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ floor(input.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Ceil(; name)

Output the ceiling rounding of the input.

# Connectors:

  - `input`
  - `output`
"""
@component function Ceil(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ ceil(input.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    StaticNonLinearity(func; name)

Applies the given function to the input.

If the given function is not composed of simple core methods (e.g. sin, abs, ...), it has to be registered via `@register_symbolic func(u)`

# Connectors:

  - `input`
  - `output`
"""
@component function StaticNonLinearity(; name, func = nothing)
    @named siso = SISO()
    @unpack u, y = siso

    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        y ~ func(u)
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, siso)
end
StaticNonLinearity(func; name) = StaticNonLinearity(; func, name)

"""
    Abs(; name)

Output the absolute value of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Abs(; name) = StaticNonLinearity(abs; name)

"""
    Sign(; name)

Output the sign of the input

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sign(; name) = StaticNonLinearity(sign; name)

"""
    Sqrt(; name)

Output the square root of the input (input >= 0 required).

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sqrt(; name) = StaticNonLinearity(sqrt; name)

"""
    Sin(; name)

Output the sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sin(; name) = StaticNonLinearity(sin; name)

"""
    Cos(; name)

Output the cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Cos(; name) = StaticNonLinearity(cos; name)

"""
    Tan(; name)

Output the tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Tan(; name) = StaticNonLinearity(tan; name)

"""
    Asin(; name)

Output the arc sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Asin(; name) = StaticNonLinearity(asin; name)

"""
    Acos(; name)

Output the arc cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Acos(; name) = StaticNonLinearity(acos; name)

"""
    Atan(; name)

Output the arc tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Atan(; name) = StaticNonLinearity(atan; name)

"""
    Atan2(; name)

Output the arc tangent of the input.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@component function Atan2(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ atan(input1.u, input2.u)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Sinh(; name)

Output the hyperbolic sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sinh(; name) = StaticNonLinearity(sinh; name)

"""
    Cosh(; name)

Output the hyperbolic cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Cosh(; name) = StaticNonLinearity(cosh; name)

"""
    Tanh(; name)

Output the hyperbolic tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Tanh(; name) = StaticNonLinearity(tanh; name)

"""
    Exp(; name)

Output the exponential (base e) of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Exp(; name) = StaticNonLinearity(exp; name)

"""
    Log(; name)

Output the natural (base e) logarithm of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Log(; name) = StaticNonLinearity(log; name)

"""
    Log10(; name)

Output the base 10 logarithm of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Log10(; name) = StaticNonLinearity(log10; name)
