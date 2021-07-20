using Base: NamedTuple
################################################################################
#
#    This file contains utility functions
#           1. Functions for parameter containers (redefined NamedTuples)
#               1.1 Adding parameter containers
#               1.2 Converting parameter containers into vectors.
#           2. Positions
################################################################################
import Base.+
import Base.-
####################### Parameters #############################################
# Abstract type
import Base.NamedTuple as Parameter

# Operators

+(p::Parameter, q::Parameter) = merge(p, q)


## Convert parameter (named tuple) to vector for use in fitting
function parameter_to_vec( p :: Parameter ) 
    l = length(p)
    v = Vector{Float64}()
    for (k, vals) in zip(keys(p), p)
        if typeof(vals) <: AbstractFloat
            append!(v, vals)
        else
            append!(v, vec(vals))
        end
    end
    return v
end

## Convert vector back to parameter for use after fitting
function vec_to_parameter!(p::Parameter, v::Vector{Float64})  
    i = 1
    j = 1
    t = Vector{Any}(undef, length(p))
    for val in p
        l = length(val)
        if l == 1
            t[j] =  v[i]
            i += 1
        else
            s = size(val)
            t[j] = reshape( v[i:(i+l-1)], s)
            i += l
        end
        j += 1
    end
    p = Parameter{keys(p)}(t)
    return p
end



####################### Positions #############################################
struct Position 
    x :: Float64
    y :: Float64 
    z :: Float64
    type ::Symbol
    Position(x, y, z) = new(x, y, z)
end
function norm(r::Position)
    return sqrt(r.x^2 + r.y^2 + r.z^2)
end

+(r1::Position, r2::Position) = Position(r1.x+r2.x, r1.y + r2.y, r1.z + r2.z)
-(r1::Position, r2::Position) = Position(r1.x-r2.x, r1.y-r2.y, r1.z-r2.z)
