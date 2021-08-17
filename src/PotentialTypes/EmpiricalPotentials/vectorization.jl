############################## Energy ################################

function potential_energy(r::Vector{Position}, p::EmpiricalPotential)
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

function potential_energy(c::Configuration, p::EmpiricalPotential)
    return potential_energy(c.Positions, p)
end

function potential_energy(r::Vector{Configuration}, p::EmpiricalPotential)
    n = length(r)
    p_e = Vector{Real}(undef, n)
    for i = 1:n
        p_e[i] = potential_energy(r[i], p)
    end
    return p_e
end

############################## Force ################################

function force(r::Vector{Position}, p::EmpiricalPotential; rcut = 2.25)
    n = length(r)
    # f = Array{Float64}(undef, n, 3)
    f = [zeros(Real, 3) for j = 1:n]
    for i = 1:n
        for j = (i+1):n
            rtemp = r[i] - r[j]
            if (norm(rtemp) < 1e-8) 
                continue
            else
                f[i] += force(rtemp, p) 
                f[j] -= f[i]
            end
        end
    end
    return f
end

function force(c::Configuration, p::EmpiricalPotential)
    return force(c.Positions, p; rcut = maximum(c.radii))
end

function force(r::Vector{Configuration}, p::EmpiricalPotential)
    n = length(r)
    f = [[Vector{Real}(undef, 3) for j = 1:r[i].num_atoms] for i = 1:n]
    # println("f ", f)
    for i = 1:n
        # println("force ", force(r[i], p))
        f[i] = force(r[i], p)
    end
    return f
end

############################ Virial (Stress) ##############################

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