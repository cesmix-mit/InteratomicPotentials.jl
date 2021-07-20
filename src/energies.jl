################################################################################
#
#    This file contains types potential energies
#    for a variety of empirical atomic potentials
#
################################################################################


############################## Lennard Jones ###################################

function potential_energy(r::Position, p::LennardJones)
    d = p.σ / norm(r)
    return 4.0 * p.ϵ * ( d^12 - d^6 )
end

function force(r::Position, p::LennardJones)
    d = norm(r)
    return 12.0 * 4.0 * ( (p.σ / d)^13 - 0.5*(p.σ / d)^7 ) .* [r.x, r.y, r.z] ./ d
end
##############################   Born-Mayer  ###################################

function potential_energy(r::Position, p::BornMayer)
    return p.A * exp(-norm(r) / p.ρ)
end

function force(r::Position, p::BornMayer)
    d = norm(r) 
    return p.A * d * exp(-d / p.ρ ) .* [r.x, r.y, r.z] ./ d
end

##############################   Coulomb  ###################################

function potential_energy(r::Position, p::Coulomb)
    return p.q_1 * p.q_2 / (4.0 * π * p.ε0 * norm(r))
end

function force(r:: Position, p::Coulomb)
    d = norm(r)
    return p.q_1 * p.q_2 / (4.0 * π * p.ε0 * d) .* [r.x, r.y, r.z] ./ d
end

##############################   GaN  ###################################

function potential_energy(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return potential_energy(r, p.c) + potential_energy(r, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return potential_energy(r, p.c) + potential_energy(r, p.lj_N_N)
    else 
        return potential_energy(r, p.c) + potential_energy(r, p.bm_Ga_N)
    end
end

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
function potential_energy(r::Vector{Position}, p::SNAP)
    println("Not yet implmented.")
    return 0.0
end

function force(r::Vector{Position}, p::SNAP)
    println("Not yet implemented.")
    return 0.0
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