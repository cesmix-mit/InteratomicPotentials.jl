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








