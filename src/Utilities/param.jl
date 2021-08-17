####################### Parameters #############################################
# Abstract type
import Base.NamedTuple as Parameter

# Operators
⊕(p::Parameter, q::Parameter) = merge(p, q)
≅(p::Tuple, q::Tuple) = length(p) == length(q) && all([size(pp) == size(qq) for (pp, qq) in zip(p, q)])
≅(x::NamedTuple{N,T}, y::NamedTuple{N2,T2}) where {N,T,N2,T2} =
  all( [key1 == key2 for (key1, key2) in zip(keys(x), keys(y))] ) && (values(x) ≅ values(y))

+(p::Tuple, q::Tuple) = (length(p) == length(q)) ? [v1 + v2 for (v1, v2) in zip(values(p), values(q))] : AssertionError("p and q must have the same length and sizes.")
+(p::Parameter, q::Parameter) = ( p ≅ q ) ? Parameter{Tuple(keys(p))}( values(p) + values(q) ) : AssertionError("p and q must have the same keys and sizes.")

-(p::Tuple, q::Tuple) = (length(p) == length(q)) ? [v1 - v2 for (v1, v2) in zip(values(p), values(q))] : AssertionError("p and q must have the same length and sizes.")
-(p::Parameter, q::Parameter) = ( p ≅ q ) ? Parameter{Tuple(keys(p))}( values(p) - values(q) ) : AssertionError("p and q must have the same keys and sizes.")
   
*(a::Real, p::Parameter) = Parameter{Tuple(keys(p))}( a .* values(p) )
        
        


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