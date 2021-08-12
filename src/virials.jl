################################################################################
#
#    This file contains virial and virial stress tensor methods for a variety of empirical atomic potentials  
#       Because virial depends on the force, a lot of the work done for typing force
#       methods can be used to simplify virial implementations.
#       Right now the methods are divided across simple potentials (type::Potential),
#       the gan potential (type::Gan <: MixedPotential) and the snap potential (type::SNAP <: FittedPotential).
################################################################################

#############################################################################
############################## Simple Virial ################################
#############################################################################
function virial(r::Position, p::EmpiricalPotential)
    f = force(r, p)
    return f[1] * r.x + f[2] * r.y + f[3] * r.z
end

function virial(r::Vector{Position}, p::EmpiricalPotential)
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

function virial(c::Configuration, p::EmpiricalPotential)
    return virial(c.Positions, p)
end

function virial(r::Vector{Configuration}, p::EmpiricalPotential)
    n = length(r)
    v = Vector{Real}(undef, n)
    for i = 1:n 
        v[i] = virial(r[i], p)
    end
    return v
end

function virial_stress(r::Position, p::EmpiricalPotential)
    f = force(r, p)
    vi = [r.x, r.y, r.z] * f'
    return [vi[1, 1], vi[2, 2], vi[3, 3], vi[2, 3], vi[1, 3], vi[1,2]]
end


function virial_stress(r::Vector{Position}, p::EmpiricalPotential)
    n = length(r)
    v = Vector{Real}(undef, 6)
    v[:] .= 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            v += virial_stress(rtemp, p)
        end
    end
    return v
end

function virial_stress(c::Configuration, p::EmpiricalPotential)
    return virial_stress(c.Positions, p)
end

function virial_stress(r::Vector{Configuration}, p::EmpiricalPotential)
    n = length(r)
    v = [Vector{Real}(undef, 6) for i = 1:n]
    for i = 1:n
        v[i] = virial_stress(r[i], p)
    end
    return v
end

#####################################################################
####################### Mixed Potentials (GaN) ######################
#####################################################################
function virial(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    f = force(r, p, type1, type2)
    return f[1] * r.x + f[2] * r.y + f[3] * r.z
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

function virial(c::Configuration, p::MixedPotential)
    return virial(c.Positions, p)
end

function virial(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    v = zeros(n)
    for i = 1:n
        v[i] = virial(r[i].Positions, p)
    end
    return v
end

function virial_stress(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    f = force(r, p, type1, type2)
    vi = [r.x, r.y, r.z] * f'
    return [vi[1, 1], vi[2, 2], vi[3, 3], vi[2, 3], vi[1, 3], vi[1,2]]
end

function virial_stress(r::Vector{Position}, p::MixedPotential)
    n = length(r)
    v = zeros(6)
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            v +=  virial_stress(rtemp, p, r[i].type, r[j].type)
        end
    end
    return v
end

function virial_stress(c::Configuration, p::MixedPotential)
    return virial_stress(c.Positions, p)
end

function virial_stress(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    v = [zeros(6) for i = 1:n]
    for i = 1:n
        v[i] = virial_stress(r[i].Positions, p)
    end
    return v
end



#######################################################################
########################### SNAP ######################################
#######################################################################

function virial(c:: Configuration, p::SNAP)
    A = get_snap(c, p)
    virial = A[end-5:end, :] * p.β
    return sum(virial[1:3])
end

function virial(r:: Vector{Configuration}, p::SNAP)
    n = length(r)
    v = Vector{Real}(undef, n)
    for j = 1:n
        v[j] = virial(r[j], p)
    end
    return v
end

function virial_stress(c:: Configuration, p::SNAP)
    A = get_snap(c, p)
    virial = A[end-5:end, :] * p.β
    return virial
end

function virial_stress(r:: Vector{Configuration}, p::SNAP)
    n = length(r)
    v = [Vector{Real}(undef, 6) for i = 1:n]
    for i = 1:n
        v[i] = virial_stress(r[i], p)
    end
    return v
end



