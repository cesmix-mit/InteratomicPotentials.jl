################################################################################
#
#    This file contains methods for gradients of energies, forces, and stresses 
#    (virial) with respect to the parameters of each method.
#
################################################################################

############################## Lennard Jones ###################################

function grad_potential_energy(r::Position, p::LennardJones)
    d = p.σ / norm(r)
    g = (dpdϵ = 4.0 * (d^12 - d^6), dpdσ = 24.0 * p.ϵ / p.σ * (2*d^12 - d^6) )

    return g
end

function grad_force(r::Position, p::LennardJones)
    d = norm(r)
    g = (dfdϵ = 48.0 * ( (p.σ / d)^13 - 0.5*(p.σ / d)^7 ) .* [r.x, r.y, r.z] ./ d, 
    dfdσ = 144.0 * p.ϵ / p.σ * ( 4.0*p.σ^12 / d^13 - p.σ^6 / d^7 ) .* [r.x, r.y, r.z] ./ d )
    return g 
end

function grad_virial(r::Position, p::LennardJones)
    df = grad_force(r, p)
    dfde = df[1]
    dfdsig = df[2]
    
    v = (dvdϵ = dfde[1]*r.x + dfde[2]*r.y + dfde[3]*r.z, 
    dvdσ = dfdsig[1]*r.x + dfdsig[2]*r.y + dfdsig[3]*r.z)
    return v

end
##############################   Born-Mayer  ###################################

function grad_potential_energy(r::Position, p::BornMayer)
    d = norm(r)
    return (dpdA = exp(-d / p.ρ),
            dpdρ = p.A * d * exp(-d/p.ρ) / p.ρ^2 )
end

function grad_force(r::Position, p::BornMayer)
    d = norm(r)
    return (dfdA = 1.0 /p.ρ * exp(-d / p.ρ ) .* [r.x, r.y, r.z] ./ d, 
            dfdρ = p.A / p.ρ^3 * exp(-d / p.ρ )*(d - p.ρ)  .* [r.x, r.y, r.z] ./ d)
end

function grad_virial(r::Position, p::BornMayer)
    df = grad_force(r, p)
    dfdA = df[:dfdA]
    dfdrho = df[:dfdρ]
    return (dvdA = dfdA[1]*r.x + dfdA[2]*r.y + dfdA[3]*r.z, 
            dvdρ = dfdrho[1]*r.x + dfdrho[2]*r.y + dfdrho[3]*r.z)
end
##############################   Coulomb  ###################################

function grad_potential_energy(r::Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dpdnull = 0.0)
end
function grad_force(r:: Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dfdnull =  0 .* [r.x, r.y, r.z] ./ d )
end

function grad_virial(r::Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dvdnull =  0 .* [r.x, r.y, r.z] ./ d )
end
##############################   GaN  ###################################

function grad_potential_energy(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return grad_potential_energy(r, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return grad_potential_energy(r, p.lj_N_N)
    else 
        return grad_potential_energy(r, p.bm_Ga_N)
    end
end

function grad_force(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    
    if (type1 == :Ga) && (type2 == :Ga) # Ga-Ga interaction
        return grad_force(r, p.lj_Ga_Ga)
    elseif (type1 == :N) && (type2 == :N) # N-N interaction
        return grad_force(r, p.lj_N_N)
    else 
        return grad_force(r, p.bm_Ga_N)
    end
end

function grad_virial(r::Position, p::GaN, type1::Symbol, type2::Symbol)
    df = grad_force(r, p, type1, type2)
    if type1 != type2
        dfdA = df[:dfdA]
        dfdrho = df[:dfdρ]
        return (dvdA = dfdA[1]*r.x + dfdA[2]*r.y + dfdA[3]*r.z, 
            dvdρ = dfdrho[1]*r.x + dfdrho[2]*r.y + dfdrho[3]*r.z)
    else
        dfde = df[:dfdϵ]
        dfdsig = df[:dfdσ]
    return (dvdϵ = dfde[1]*r.x + dfde[2]*r.y + dfde[3]*r.z, 
            dvdσ = dfdsig[1]*r.x + dfdsig[2]*r.y + dfdsig[3]*r.z)
    end
end


##############################  SNAP  ###################################
function grad_potential_energy(r::Position, p::SNAP)
    println("Not yet implmented.")
    return 0.0
end

function grad_force(r::Position, p::SNAP)
    println("Not yet implemented.")
    return 0.0
end

function grad_virial(r::Position, p::SNAP)
    println("Not yet implemented")
    return 0.0
end
############################## Vectorize ################################
function grad_potential_energy(r::Vector{Position}, p::EmpiricalPotential)
    n = length(r)
    pe = 0. * grad_potential_energy(r[1] - r[2], p)
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            pe +=  grad_potential_energy(rtemp, p)
        end
    end
    return pe
end

function grad_potential_energy(r::Vector{Position}, p::GaN)
    n = length(r)
    # Initialize 
    rtemp = r[1] - r[2]
    pe = [0. * grad_potential_energy(rtemp, p, :Ga, :Ga), 0. * grad_potential_energy(rtemp, p, :N, :N), 0. * grad_potential_energy(rtemp, p, :Ga, :N)]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            if r[i].type == r[j].type 
                if r[i].type == :Ga
                    pe[1] += grad_potential_energy(rtemp, p, r[i].type, r[j].type)
                else
                    pe[2] += grad_potential_energy(rtemp, p, r[i].type, r[j].type)
                end
            else
                pe[3] += grad_potential_energy(rtemp, p, r[i].type, r[j].type)
            end
        end
    end
    return pe
end

function grad_force(r::Vector{Position}, p::EmpiricalPotential)
    n = length(r)
    f = [0. * grad_force(r[1] - r[2], p) for j = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            f[i] +=  grad_force(rtemp, p) 
            f[j] -= f[i]
        end
    end
    return f
end

function grad_force(r::Vector{Position}, p::GaN)
    # Initialize 
    rtemp = r[1] - r[2]
    f = [ [0. * grad_force(rtemp, p, :Ga, :Ga), 0. * grad_force(rtemp, p, :N, :N), 0. * grad_force(rtemp, p, :Ga, :N)] for j = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
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

function grad_virial(r::Vector{Position}, p::EmpiricalPotential)
    n = length(r)
    v = 0. * grad_virial(r[1] - r[2], p)
    for i = 1:(n-1)
        for j = (i+1):n
            rtemp = r[i] - r[j]
            v +=  grad_virial(rtemp, p)
        end
    end
    return v
end

function grad_virial(r::Vector{Position}, p::GaN)
    n = length(r)
    # Initialize 
    rtemp = r[1] - r[2]
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