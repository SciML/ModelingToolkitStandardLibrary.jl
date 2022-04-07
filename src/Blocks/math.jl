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
    nout, nin = size(K)
    @named miso = MIMO(;nin=nin, nout=nout)
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
        output.u ~ sum(input.u[i] for i in 1:n) # FIXME: if array equations work
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

# TODO:
# Cos	Output the cosine of the input
# Tan	Output the tangent of the input
# Asin	Output the arc sine of the input
# Acos	Output the arc cosine of the input
# Atan	Output the arc tangent of the input
# Atan2	Output atan(u1/u2) of the inputs u1 and u2
# Sinh	Output the hyperbolic sine of the input
# Cosh	Output the hyperbolic cosine of the input
# Tanh	Output the hyperbolic tangent of the input
# Exp	Output the exponential (base e) of the input
# Log	Output the natural (base e) logarithm of the input (input > 0 required)
# Log10 Output the base 10 logarithm of the input (input > 0 required)