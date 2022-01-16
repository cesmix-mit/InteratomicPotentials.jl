# ############################## Energy ################################

function potential_energy(a::AtomsBase.Atom, p::EmpiricalPotential)
    return potential_energy(ustrip.(a.position), p)
end

function potential_energy(a::AbstractSystem, p::EmpiricalPotential)
    pe = 0.0
    r = a.particles
    N = length(r)
    for i = 1:(N-1)
        ai = ustrip.(r[i].position)
        for j = (i+1):N
            aj = ustrip.(r[j].position)
            rtemp = ai - aj
            pe +=  potential_energy(rtemp, p)
        end
    end
    return pe
end

# ############################## Force ################################
function force(a::AtomsBase.Atom, p::EmpiricalPotential)
    return force(ustrip.(a.position), p)
end

function force(s::AbstractSystem, p::EmpiricalPotential)
    r = s.particles
    n = length(r)
    f = [ zeros(3) for j = 1:n]
    for i = 1:n
        ai = ustrip.(r[i].position)
        for j = (i+1):n
            aj = ustrip.(r[j].position)
            rtemp = ai - aj
            if (norm(rtemp) < 1e-8) 
                continue
            else
                f[i] += force(rtemp, p)
                f[j] -= f[i]
            end
        end
    end
    return SVector{n}([SVector{3}(fi) for fi in f])
end
# ############################## Virial ################################
function virial(a::AtomsBase.Atom, p::EmpiricalPotential)
    return virial(ustrip.(a.position), p)
end

function virial(S::AbstractSystem, p::EmpiricalPotential)
    r = S.particles
    n = length(r)
    v = 0.0
    for i = 1:(n-1)
        ri = ustrip.(r[i].position)
        for j = (i+1):n
            rj = ustrip.(r[j].position)
            rtemp = ri - rj
            v +=  virial(rtemp, p)
        end
    end
    return v
end

function virial_stress(a::AtomsBase.Atom, p::EmpiricalPotential)
    return virial_stress(ustrip.(a.position), p)
end

function virial_stress(S::AbstractSystem, p::EmpiricalPotential)
    v = zeros(Real, 6)
    r = S.particles
    n = length(r)
    for i = 1:(n-1)
        ri = ustrip.(r[i].position)
        for j = (i+1):n
            rj = ustrip.(r[j].position)
            rtemp = ri - rj
            v += virial_stress(rtemp, p)
        end
    end
    return SVector{6}(v)
end
# function potential_energy(a::Atom, p::EmpiricalPotential)
#     return potential_energy(a.Position, p)
# end

# function potential_energy(r::SVector{N, Atom} where N <: Int, p::EmpiricalPotential)
#     pe = 0.0
#     # N = length(r)
#     for i = 1:(N-1)
#         ai = r[i].Position
#         for j = (i+1):N
#             aj = r[j].Position
#             rtemp = ai - aj
#             pe +=  potential_energy(rtemp, p)
#         end
#     end
#     return pe
# end

# function potential_energy(c::Configuration, p::EmpiricalPotential)
#     return potential_energy(c.Atoms, p)
# end

# function potential_energy(r::SVector{Configuration}, p::EmpiricalPotential)
#     n = length(r)
#     p_e = Vector{Real}(undef, n)
#     for i = 1:n
#         p_e[i] = potential_energy(r[i], p)
#     end
#     return SA[p_e[i] for i = 1:n]
# end

# ############################## Force ################################

# function force(a::Atom, p::EmpiricalPotential)
#     return force(a.Position, p)
# end

# function force(r::SVector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     f = [zeros(Real, 3) for j = 1:n]
#     for i = 1:n
#         ai = r[i].Position
#         for j = (i+1):n
#             aj = r[j].Position
#             rtemp = ai - aj
#             if (norm(rtemp) < 1e-8) 
#                 continue
#             else
#                 f[i] += force(rtemp, p) 
#                 f[j] -= f[i]
#             end
#         end
#     end
#     return f
# end

# function force(c::Configuration, p::EmpiricalPotential)
#     return force(c.Atoms, p)
# end

# function force(r::Vector{Configuration}, p::EmpiricalPotential)
#     n = length(r)
#     f = [[Vector{Real}(undef, 3) for j = 1:length(r[i].Atoms)] for i = 1:n]
#     for i = 1:n
#         f[i] = force(r[i], p)
#     end
#     return f
# end

# ############################ Virial (Stress) ##############################

# function virial(a::Atom, p::EmpiricalPotential)
#     return virial(a.Position, p)
# end

# function virial(r::Vector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     v = 0.0
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             v +=  virial(rtemp, p)
#         end
#     end
#     return v
# end

# function virial(c::Configuration, p::EmpiricalPotential)
#     return virial(c.Atoms, p)
# end

# function virial(r::Vector{Configuration}, p::EmpiricalPotential)
#     n = length(r)
#     v = Vector{Real}(undef, n)
#     for i = 1:n 
#         v[i] = virial(r[i], p)
#     end
#     return v
# end

# function virial_stress(a::Atom, p::EmpiricalPotential)
#     return virial_stress(a.Position, p)
# end


# function virial_stress(r::Vector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     v = Vector{Real}(undef, 6)
#     v[:] .= 0.0
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             v += virial_stress(rtemp, p)
#         end
#     end
#     return v
# end

# function virial_stress(c::Configuration, p::EmpiricalPotential)
#     return virial_stress(c.Atoms, p)
# end

# function virial_stress(r::Vector{Configuration}, p::EmpiricalPotential)
#     n = length(r)
#     v = [Vector{Real}(undef, 6) for i = 1:n]
#     for i = 1:n
#         v[i] = virial_stress(r[i], p)
#     end
#     return v
# end

# ############################## Gradients ################################
# function grad_potential_energy(r::Vector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     pe = 0. * grad_potential_energy(r[1].Position - r[2].Position, p)
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             pe +=  grad_potential_energy(rtemp, p)
#         end
#     end
#     return pe
# end



# function grad_force(r::Vector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     f = [0. * grad_force(r[1].Position - r[2].Position, p) for j = 1:n]
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             f[i] += grad_force(rtemp, p) 
#             f[j] -= f[i]
#         end
#     end
#     return f
# end



# function grad_virial(r::Vector{Atom}, p::EmpiricalPotential)
#     n = length(r)
#     v = 0. * grad_virial(r[1].Position - r[2].Position, p)
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             v +=  grad_virial(rtemp, p)
#         end
#     end
#     return v
# end
