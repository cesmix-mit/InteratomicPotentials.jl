#########################################################################
##############################   GaN  ###################################
#########################################################################

function force(r::Vector{<:Real}, p::GaN, Type1::Symbol, Type2::Symbol)
    
    if (Type1 == :Ga) && (Type2 == :Ga) # Ga-Ga interaction
        return force(r, p.c) + force(r, p.lj_Ga_Ga)
    elseif (Type1 == :N) && (Type2 == :N) # N-N interaction
        return force(r, p.c) + force(r, p.lj_N_N)
    else 
        return force(r, p.c) + force(r, p.bm_Ga_N)
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
            f[i] +=  force(rtemp, p, r[i].Type, r[j].Type) 
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
    
    f = [[zeros(3) for j = 1:length(r[i].Atoms)] for i = 1:n]
    
    for i = 1:n
        f[i] = force(r[i], p)
    end
    return f
end







