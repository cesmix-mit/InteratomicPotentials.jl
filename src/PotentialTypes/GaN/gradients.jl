##############################  Gradients  ###################################

function grad_potential_energy(r::Atom, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return grad_potential_energy(r.Position, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return grad_potential_energy(r.Position, p.lj_N_N)
    else 
        return grad_potential_energy(r.Position, p.bm_Ga_N)
    end
end

function grad_potential_energy(r::Vector{Atom}, p::GaN)
    n = length(r)
    # Initialize 
    rtemp = r[1].Position - r[2].Position
    pe = [0. * grad_potential_energy(rtemp, p, :Ga, :Ga), 0. * grad_potential_energy(rtemp, p, :N, :N), 0. * grad_potential_energy(rtemp, p, :Ga, :N)]
    for i = 1:(n-1)
        ri = r[i].Position
        ri_type = r[i].type
        for j = (i+1):n
            rj = r[j].Position
            rj_type = r[j].type
            rtemp =  ri - rj
            if ri_type == r[j].type 
                if r[i].type == :Ga
                    pe[1] += grad_potential_energy(rtemp, p, ri_type, rj_type)
                else
                    pe[2] += grad_potential_energy(rtemp, p, ri_type, rj_type)
                end
            else
                pe[3] += grad_potential_energy(rtemp, p, ri_type, rj_type)
            end
        end
    end
    return pe
end

function grad_force(r::Atom, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return grad_force(r.Position, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return grad_force(r.Position, p.lj_N_N)
    else 
        return grad_force(r.Position, p.bm_Ga_N)
    end
end

function grad_force(r::Vector{Atom}, p::GaN)
    # Initialize 
    rtemp = r[1].Position - r[2].Position
    f = [ [0. * grad_force(rtemp, p, :Ga, :Ga), 0. * grad_force(rtemp, p, :N, :N), 0. * grad_force(rtemp, p, :Ga, :N)] for j = 1:n]
    for i = 1:(n-1)
        ri = r[i].Position
        for j = (i+1):n
            rj = r[j].Position
            rtemp = ri - rj
            if r[i].type == r[j].type 
                if r[i].type == :Ga
                    f[i][1] += grad_force(rtemp, p, r[i].type, r[j].type)
                    f[j][1] -= f[i][1]
                else
                    f[i][2] += grad_force(rtemp, p, r[i].type, r[j].type)
                    f[j][2] -= f[i][2]
                end
            else
                f[i][3] += grad_force(rtemp, p, r[i].type, r[j].type)
                f[j][3] -= f[i][3]
            end
        end
    end
    return f
end

function grad_virial(r::Atom, p::GaN, type1::Symbol, type2::Symbol)
    df = grad_force(r, p, type1, type2)
    if type1 != type2
        dfdA = df[:dfdA]
        dfdrho = df[:dfdρ]
        return (dvdA = dfdA[1]*r.Position[1] + dfdA[2]*r.Position[2] + dfdA[3]*r.Position[3], 
            dvdρ = dfdrho[1]*r.Position[1] + dfdrho[2]*r.Position[2] + dfdrho[3]*r.Position[3])
    else
        dfde = df[:dfdϵ]
        dfdsig = df[:dfdσ]
    return (dvdϵ = dfde[1]*r.Position[1] + dfde[2]*r.Position[2] + dfde[3]*r.Position[3], 
            dvdσ = dfdsig[1]*r.Position[1] + dfdsig[2]*r.Position[2] + dfdsig[3]*r.Position[3])
    end
end

function grad_virial(r::Vector{Atom}, p::GaN)
    n = length(r)
    # Initialize 
    rtemp = r[1].Position - r[2].Position
    pe = [0. * grad_virial(rtemp, p, :Ga, :Ga), 0. * grad_virial(rtemp, p, :N, :N), 0. * grad_virial(rtemp, p, :Ga, :N)]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            if r[i].type == r[j].type 
                if r[i].type == :Ga
                    pe[1] += grad_virial(rtemp, p, r[i].type, r[j].type)
                else
                    pe[2] += grad_virial(rtemp, p, r[i].type, r[j].type)
                end
            else
                pe[3] += grad_virial(rtemp, p, r[i].type, r[j].type)
            end
        end
    end
    return pe
end