# ####################### Parameters #############################################
# import Base.NamedTuple as Parameter
# import Base.+
# function +(p::Parameter, q::Parameter)
#     names_p = keys(p)
#     names_q = keys(q)
#     new_names_q = []
#     for (i, key) in enumerate(names_q)
#         j = 1
#         while !isempty( intersect(names_p, [key] ) )
#             key = Symbol( string(key) * "_$j") 
#             j+=1
#         end
#         println(new_names_q)
#         println(key)
#         new_names_q = (new_names_q..., (key, )...)
#     end

#     merge(Parameter{names_p}(values(p)), Parameter{new_names_q}(values(q)))
# end
        





## Convert parameter (named tuple) to vector for use in fitting
function parameter_to_vec(p::Parameter)
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
            t[j] = v[i]
            i += 1
        else
            s = size(val)
            t[j] = reshape(v[i:(i+l-1)], s)
            i += l
        end
        j += 1
    end
    p = Parameter{keys(p)}(t)
    return p
end