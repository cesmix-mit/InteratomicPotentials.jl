################################################################################
#
#    This file contains virial methods for a variety of empirical atomic potentials  
#       Because virial depends on the force, a lot of the work done for typing force
#       methods can be used to simplify virial implementations.
################################################################################

############################## Virial ################################
function virial(r::Position, p::Potential)
    f = force(r, p)
    return f[1] * r.x + f[2] * r.y + f[3] * r.z
end

function virial(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    f = force(r, p, type1, type2)
    return f[1] * r.x + f[2] * r.y + f[3] * r.z
end

############################ Vectorize ###############################
function virial(r::Vector{Position}, p::Potential)
    n = length(r)
    v = 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            v +=  virial(rtemp, p)
        end
    end
    return v
end

function virial(r::Vector{Position}, p::MixedPotential)
    n = length(r)
    v = 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            v +=  virial(rtemp, p, r[i].type, r[j].type)
        end
    end
    return v
end