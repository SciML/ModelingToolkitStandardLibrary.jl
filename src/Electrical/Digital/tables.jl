# AND gate
AndTable = [
            # U X F0 F1 Z W L H DC
            U U F0 U U U F0 U U        # U
            U X F0 X X X F0 X X        # X
            F0 F0 F0 F0 F0 F0 F0 F0 F0 # F0
            U X F0 F1 X X F0 F1 X      # F1
            U X F0 X X X F0 X X        # Z
            U X F0 X X X F0 X X        # W
            F0 F0 F0 F0 F0 F0 F0 F0 F0 # L
            U X F0 F1 X X F0 F1 X      # H
            U X F0 X X X F0 X X]       # DC

# NOT gate
NotTable = [U, X, F1, F0, X, X, F1, F0, X]

# OR gate
OrTable = [
           # U X F0 F1 Z W L H DC
           U U U F1 U U U F1 U        # U
           U X X F1 X X X F1 X        # X
           U X F0 F1 X X F0 F1 X      # F0
           F1 F1 F1 F1 F1 F1 F1 F1 F1 # F1
           U X X F1 X X X F1 X        # Z
           U X X F1 X X X F1 X        # W
           U X F0 F1 X X F0 F1 X      # L
           F1 F1 F1 F1 F1 F1 F1 F1 F1 # H
           U X X F1 X X X F1 X]       # DC

# XOR gate
XorTable = [
            # U X F0 F1 Z W L H DC
            U U U U U U U U U     # U
            U X X X X X X X X     # X
            U X F0 F1 X X F0 F1 X # F0
            U X F1 F0 X X F1 F0 X # F1
            U X X X X X X X X     # Z
            U X X X X X X X X     # W
            U X F0 F1 X X F0 F1 X # L
            U X F1 F0 X X F1 F0 X # H
            U X X X X X X X X]    # DC

function _xor2(a, b)
    XorTable[logic_level_dict[a], logic_level_dict[b]]
end

function _xor(x...)
    y = [x[1]]
    for i in 2:lastindex(x)
        push!(y, _xor2(x[i], y[i - 1]))
    end
    return y[end]
end

@register_symbolic _xor(x...)
@register_symbolic _xor(a, b)
@register_symbolic _xor(a, b, c)
@register_symbolic _xor(a, b, c, d)
@register_symbolic _xor(a, b, c, d, e)
