# 9-level logic
@parameters X

logicmap = Dict(0 => 0,
                1 => 1,
                :U => 2, # Uninitialized
                :X => 3, # Unknown
                :Z => 4, # High Impedence
                :W => 5, # Weak Unknown
                :L => 6, # Weak Low                
                :H => 7, # Weak High
                :D => 8) # Don't Care; standard representation is `-`
                
 
logic(l) = try logicmap[l] 
  catch(ex)
    X
  end

# SISO NOT gate
NotTable = OffsetArray([
              # 0  1  :U  :X  :Z  :W  :L  :H  :D
                1; 0; :U; :X; :X; :X; 1; 0; :X
            ], 0:8)

function _not(x)
  i = logic(x)
  typeof(i) != Int && return X
  NotTable[i]
end
@register_symbolic _not(x)

# MISO AND gate
AndTable = OffsetArray([
              # 0  1 :U :X :Z :W :L :H :D
                0  0  0  0  0  0  0  0  0 # 0
                0  1 :U :X :X :X  0  1 :X # 1
                0 :U :U :U :U :U  0 :U :U # :U
                0 :X :U :X :X :X  0 :X :X # :X
                0 :X :U :X :X :X  0 :X :X # :Z
                0  1 :U :X :X :X  0 :X :X # :W
                0  0  0  0  0  0  0  0  0 # :L
                0  1 :U :X :X :X  0  1 :X # :H
                0 :X :U :X :X :X  0 :X :X # :D
            ], 0:8, 0:8)

function _and2(a, b)
  i, j = logic(a), logic(b)
  (typeof(i) != Int || typeof(j) != Int) && return X
  AndTable[i, j]
end
function _and(x...)
  y = [_and2(x[1], x[1])]
  for i in 2:length(x)
    push!(y, _and2(x[i], y[i-1]))
  end
  return y[end]
end
@register_symbolic _and(x...)
@register_symbolic _and(a, b)
@register_symbolic _and(a, b, c)
@register_symbolic _and(a, b, c, d)
@register_symbolic _and(a, b, c, d, e)

# MISO OR gate
OrTable = OffsetArray([
              # 0  1 :U :X :Z :W :L :H :D
                0  1 :U :X :X :X  0  1 :X # 0
                1  1  1  1  1  1  1  1  1 # 1
                :U 1 :U :U :U :U :U  1 :U # :U
                :X 1 :U :X :X :X :X  1 :X # :X
                :X 1 :U :X :X :X :X  1 :X # :Z
                :X 1 :U :X :X :X :X  1 :X # :W
                0  1 :U :X :X :X  0  1 :X # :L
                1  1  1  1  1  1  1  1  1 # :H
                :X 1 :U :X :X :X :X  1 :X # :D
                ], 0:8, 0:8)

function _or2(a, b)
  i, j = logic(a), logic(b)
  (typeof(i) != Int || typeof(j) != Int) && return X
  OrTable[i, j]
end

function _or(x...)
  y = [_or2(x[1], x[1])]
  for i in 2:length(x)
    push!(y, _or2(x[i], y[i-1]))
  end
  return y[end]
end
@register_symbolic _or(x...)
@register_symbolic _or(a, b)
@register_symbolic _or(a, b, c)
@register_symbolic _or(a, b, c, d)
@register_symbolic _or(a, b, c, d, e)
@register_symbolic _or(a, b, c, d, e, f, g, h)

# MISO :XOR gate

XorTable = OffsetArray([
              #  0  1 :U :X :Z :W :L :H :D
                 0  1 :U :X :X :X  0  1 :X # 0
                 1  0 :U :X :X :X  1  0 :X # 1
                :U :U :U :U :U :U :U :U :U # :U
                :X :X :U :X :X :X :X :X :X # :X
                :X :X :U :X :X :X :X :X :X # :Z
                :X :X :U :X :X :X :X :X :X # :W
                 0  1 :U :X :X :X  0  1 :X # :L
                 1  0 :U :X :X :X  1  0 :X # :H
                :X :X :U :X :X :X :X :X :X # :D
                ], 0:8, 0:8)

function _xor2(a, b)
  i, j = logic(a), logic(b)
  (typeof(i) != Int || typeof(j) != Int) && return X
  XorTable[i, j]
end

function _xor(x...)
  y = [_xor2(x[1], 0)]
  for i in 2:length(x)
    push!(y, _xor2(x[i], y[i-1]))
  end
  return y[end]
end
@register_symbolic _xor(x...)
@register_symbolic _xor(a, b)
@register_symbolic _xor(a, b, c)
@register_symbolic _xor(a, b, c, d)
@register_symbolic _xor(a, b, c, d, e)

# TODO: revisit y[1] for all miso gates for 9-level logic