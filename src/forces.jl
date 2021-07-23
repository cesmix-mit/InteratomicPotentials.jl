################################################################################
#
#    This file contains force methods for a variety of empirical atomic potentials  
#
################################################################################


############################## Lennard Jones ###################################

function force(r::Position, p::LennardJones)
    d = norm(r)
    return 12.0 * 4.0 * ( p.σ^12 / d^13 - 0.5*p.σ^6 / d^7 ) .* [r.x, r.y, r.z] ./ d
end

##############################   Born-Mayer  ###################################

function force(r::Position, p::BornMayer)
    d = norm(r) 
    return p.A /p.ρ * exp(-d / p.ρ ) .* [r.x, r.y, r.z] ./ d
end

##############################   Coulomb  ###################################

function force(r:: Position, p::Coulomb)
    d = norm(r)
    return p.q_1 * p.q_2 / (4.0 * π * p.ε0 * d) .* [r.x, r.y, r.z] ./ d
end

##############################   GaN  ###################################


function force(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return force(r, p.c) + force(r, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return force(r, p.c) + force(r, p.lj_N_N)
    else 
        return force(r, p.c) + force(r, p.bm_Ga_N)
    end
end

##############################  SNAP  ###################################

function force(c::Configuration, p::SNAP; return_positions = false)
    A = get_snap(c, p)
    force = A[2:end-6, :] * p.β
    if return_positions
        l = Int(length(force) / 3)
        force = reshape(force, 3, l)'
        force_position = Vector{Position}(undef, l)
        for i = 1:l
            force_position[i] = Position(force[i, 1], force[i, 2], force[i,3])
        end
        return force_position
    else
        return force
    end
end

############################## Vectorize ################################

function force(r::Vector{Position}, p::Potential)
    n = length(r)
    f = [zeros(3) for j = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            f[i] +=  force(rtemp, p) 
            f[j] -= f[i]
        end
    end
    return f
end

function force(r::Vector{Position}, p::MixedPotential)
    n = length(r)
    f = [zeros(3) for j = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            f[i] +=  force(rtemp, p, r[i].type, r[j].type) 
            f[j] += f[i]
        end
    end
    return f
end