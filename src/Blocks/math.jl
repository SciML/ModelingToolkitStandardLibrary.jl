"""
    Gain(k=1; name)

Output the product of a gain value with the input signal.

# Parameters:
- `k`: Scalar gain
"""
function Gain(k=1; name)
    @named siso = SISO()
    @unpack u, y = siso
    pars = @parameters k=k
    eqs = [
        y ~ k * u
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), siso)
end

"""
    MatrixGain(K::AbstractArray; name)

Output the product of a gain matrix with the input signal vector.

# Parameters:
- `K`: Matrix gain
"""
function MatrixGain(K::AbstractArray; name)
    nout, nin = size(K)
    @named input = RealInput(;nin=nin)
    @named output = RealOutput(;nout=nout)
    eqs = [
        output.u[i] ~ sum(K[i,j] * input.u[j] for j in 1:nin) for i in 1:nout # FIXME: if array equations work
    ]
    compose(ODESystem(eqs, t, [], []; name=name), [input, output])
end

"""
    Sum(n::Int; name)

Output the sum of the elements of the input port vector.

# Parameters:
- `n`: Input port dimension
"""
function Sum(n::Int; name)
    @named input = RealInput(;nin=n)
    @named output = RealOutput()
    eqs = [
        output.u ~ sum(input.u)
    ]
    compose(ODESystem(eqs, t, [], []; name=name), [input, output])
end
    
"""
    Feedback(;name)

Output difference between reference input (input1) and feedback input (input2).
"""
function Feedback(;name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs= [
        output.u ~ input1.u - input2.u
    ]
    return compose(ODESystem(eqs, t, [], []; name=name), input1, input2, output)
end

"""
    Add(;name, k1=1, k2=1)

Output the sum of the two scalar inputs.

# Parameters:
- `k1`: Gain for first input
- `k2`: Gain for second input
"""
function Add(;name, k1=1, k2=1)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    pars = @parameters begin
        k1=k1
        k2=k2
    end
    eqs= [
        output.u ~ k1 * input1.u + k2 * input2.u
    ]
    return compose(ODESystem(eqs, t, [], pars; name=name), input1, input2, output)
end

"""
    Add(;name, k1=1, k2=1,k3=1)

Output the sum of the three scalar inputs.

# Parameters:
- `k1`: Gain for first input
- `k2`: Gain for second input
- `k3`: Gain for third input
"""
function Add3(;name, k1=1, k2=1, k3=1)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named input3 = RealInput()
    @named output = RealOutput()
    pars = @parameters begin
        k1=k1
        k2=k2
        k3=k3
    end
    eqs= [
        output.u ~ k1 * input1.u + k2 * input2.u + k3 * input3.u
    ]
    return compose(ODESystem(eqs, t, [], pars; name=name), input1, input2, input3, output)
end

"""
    Product(;name)

Output product of the two inputs.
"""
function Product(;name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs= [
        output.u ~ input1.u * input2.u
    ]
    return compose(ODESystem(eqs, t, [], []; name=name), input1, input2, output)
end

"""
    Division(;name)

Output first input divided by second input.
"""
function Division(;name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs= [
        output.u ~ input1.u / input2.u
    ]
    return compose(ODESystem(eqs, t, [], []; name=name), input1, input2, output)
end


"""
    StaticNonLinearity(func ;name)

Applies the given function to the input. 

If the given function is not composed of simple core methods (e.g. sin, abs, ...), it has to be registered via `@register_symbolic func(u)`
"""
function StaticNonLinearity(func; name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [y ~ func(u)]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
    Abs(;name)

Output the absolute value of the input.
"""
Abs(;name) = StaticNonLinearity(abs; name)

"""
    Sign(;name)

Output the sign of the input
"""
Sign(;name) = StaticNonLinearity(sign; name)

"""
    Sqrt(;name)

Output the square root of the input (input >= 0 required).
"""
Sqrt(;name) = StaticNonLinearity(sqrt; name)

"""
    Sin(;name)

Output the sine of the input.
"""
Sin(;name) = StaticNonLinearity(sin; name)

"""
    Cos(;name)

Output the cosine of the input.
"""
Cos(;name) = StaticNonLinearity(cos; name)

"""
    Tan(;name)

Output the tangent of the input.
"""
Tan(;name) = StaticNonLinearity(tan; name)

"""
    Asin(;name)

Output the arc sine of the input.
"""
Asin(;name) = StaticNonLinearity(asin; name)

"""
    Acos(;name)

Output the arc cosine of the input.
"""
Acos(;name) = StaticNonLinearity(acos; name)

"""
    Atan(;name)

Output the arc tangent of the input.
"""
Atan(;name) = StaticNonLinearity(atan; name)

"""
    Atan2(;name)

Output the arc tangent of the input.
"""
function Atan2(;name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs = [
        output.u ~ atan(input1.u, input2.u)
    ]
    compose(ODESystem(eqs, t, [], []; name=name), [input1, input2, output])
end

"""
    Sinh(;name)

Output the hyperbolic sine of the input.
"""
Sinh(;name) = StaticNonLinearity(sinh; name)

"""
    Cosh(;name)

Output the hyperbolic cosine of the input.
"""
Cosh(;name) = StaticNonLinearity(cosh; name)

"""
    Tanh(;name)

Output the hyperbolic tangent of the input.
"""
Tanh(;name) = StaticNonLinearity(tanh; name)

"""
    Exp(;name)

Output the exponential (base e) of the input.
"""
Exp(;name) = StaticNonLinearity(exp; name)

"""
    Log(;name)

Output the natural (base e) logarithm of the input.
"""
Log(;name) = StaticNonLinearity(log; name)

"""
    Log10(;name)

Output the base 10 logarithm of the input.
"""
Log10(;name) = StaticNonLinearity(log10; name)