#########################################################################
##############################   GaN  ###################################
#########################################################################

function force(r::Atom, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return force(r.Position, p.c) + force(r.Position, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return force(r.Position, p.c) + force(r.Position, p.lj_N_N)
    else 
        return force(r.Position, p.c) + force(r.Position, p.bm_Ga_N)
    end
end


function force(r::Vector{Atom}, p::MixedPotential)
    n = length(r)
    # f = Array{Float64}(undef, 3, n)
    f = [zeros(3) for j = 1:n]
    for i = 1:(n-1)
        ri = r[i].Position
        for j = (i+1):n
            rj = r[j].Position
            rtemp = ri - rj
            f[i] +=  force(rtemp, p, r[i].type, r[j].type) 
            f[j] -= f[i]
        end
    end
    return f
end

function force(c::Configuration, p::MixedPotential)
    return force(c.Atoms, p)
end

function force(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    f = [[zeros(3) for j = 1:r[i].num_atoms] for i = 1:n]
    
    for i = 1:n
        f[i] = force(r[i], p)
    end
    return f
end







