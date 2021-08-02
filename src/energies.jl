################################################################################
#
#    This file contains potential energy methods 
#    for a variety of empirical atomic potentials
#
################################################################################


############################## Lennard Jones ###################################

function potential_energy(r::Position, p::LennardJones)
    d = p.σ / norm(r)
    return 4.0 * p.ϵ * ( d^12 - d^6 )
end

##############################   Born-Mayer  ###################################

function potential_energy(r::Position, p::BornMayer)
    return p.A * exp(-norm(r) / p.ρ)
end

##############################   Coulomb  ###################################

function potential_energy(r::Position, p::Coulomb)
    return p.q_1 * p.q_2 / (4.0 * π * p.ϵ0 * norm(r))
end

############################## Vectorize ################################

function potential_energy(r::Vector{Position}, p::Potential)
    n = length(r)
    pe = 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            pe +=  potential_energy(rtemp, p)
        end
    end
    return pe
end

function potential_energy(c::Configuration, p::Potential)
    return potential_energy(c.Positions, p)
end

function potential_energy(r::Vector{Configuration}, p::Potential)
    n = length(r)
    p_e = zeros(n)
    for i = 1:n
        p_e[i] = potential_energy(r[i], p)
    end
end

#########################################################################
##############################   GaN  ###################################
#########################################################################
function potential_energy(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return potential_energy(r, p.c) + potential_energy(r, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return potential_energy(r, p.c) + potential_energy(r, p.lj_N_N)
    else 
        return potential_energy(r, p.c) + potential_energy(r, p.bm_Ga_N)
    end
end

function potential_energy(r::Vector{Position}, p::MixedPotential)
    n = length(r)
    pe = 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            pe +=  potential_energy(rtemp, p, r[i].type, r[j].type)
        end
    end
    return pe
end


function potential_energy(c::Configuration, p::MixedPotential)
    return potential_energy(c.Positions, p)
end

function potential_energy(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    pe = zeros(n)
    for i = 1:n
        pe[i] = potential_energy(r[i], p)
    end
    return pe
end

#########################################################################
##############################  SNAP  ###################################
#########################################################################
function potential_energy(c::Configuration, p::SNAP)
    A = get_snap(c, p)
    p_e = dot(A[1, :]  ,  p.β )
    return p_e
end

function potential_energy(r::Vector{Configuration}, p::SNAP)
    n = length(r)
    pe = zeros(n)
    for i = 1:n
        A = get_snap(r[i], p)
        pe[i] = dot(A[1, :] , p.β)
    end
    return pe
end







