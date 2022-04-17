"""
Output the product of a gain value with the input signal.
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
Output the product of a gain matrix with the input signal vector.
"""
function MatrixGain(K::AbstractArray; name)
    ny, nu = size(K)
    @named miso = MIMO(;nin=nu, nout=ny)
    @unpack u, y = miso
    eqs = [
        y[i] ~ sum(K[i,j] * u[j] for j in 1:nin) for i in 1:nout # FIXME: if array equations work
    ]
    extend(ODESystem(eqs, t, [], []; name=name), miso)
end

"""
Output the sum of the elements of the input vector.
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
Output difference between commanded and feedback input.
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
Output the sum of the two inputs.
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
Output the sum of the three inputs.
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
Output the absolute value of the input.
"""
function Abs(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ abs(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the sign of the input
"""
function Sign(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ sign(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the square root of the input (input >= 0 required).
"""
function Sqrt(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ sqrt(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the sine of the input.
"""
function Sin(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ sin(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the cosine of the input.
"""
function Cos(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ cos(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the tangent of the input.
"""
function Tan(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ tan(uP)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the arc sine of the input.
"""
function Asin(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ asin(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the arc cosine of the input.
"""
function Acos(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ acos(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the arc tangent of the input.
"""
function Atan(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ atan(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
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
Output the hyperbolic sine of the input.
"""
function Sinh(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ sinh(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the hyperbolic cosine of the input.
"""
function Cosh(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ cosh(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the hyperbolic tangent of the input.
"""
function Tanh(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ tanh(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the exponential (base e) of the input.
"""
function Exp(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ exp(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the natural (base e) logarithm of the input.
"""
function Log(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ log(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end

"""
Output the base 10 logarithm of the input.
"""
function Log10(;name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [
        y ~ log10(u)
    ]
    extend(ODESystem(eqs, t, [], []; name=name), siso)
end