# ##############################  Gradients  ###################################

# function grad_potential_energy(r::Atom, p::GaN, Type1::Symbol, Type2::Symbol)
    
#     if (Type1 == :Ga) && (Type2 == :Ga) # Ga-Ga interaction
#         return grad_potential_energy(r.Position, p.lj_Ga_Ga)
#     elseif (Type1 == :N) && (Type2 == :N) # N-N interaction
#         return grad_potential_energy(r.Position, p.lj_N_N)
#     else 
#         return grad_potential_energy(r.Position, p.bm_Ga_N)
#     end
# end

# function grad_potential_energy(r::Vector{Atom}, p::GaN)
#     n = length(r)
#     # Initialize 
#     rtemp = r[1].Position - r[2].Position
#     pe = [0. * grad_potential_energy(rtemp, p, :Ga, :Ga), 0. * grad_potential_energy(rtemp, p, :N, :N), 0. * grad_potential_energy(rtemp, p, :Ga, :N)]
#     for i = 1:(n-1)
#         ri = r[i].Position
#         ri_Type = r[i].Type
#         for j = (i+1):n
#             rj = r[j].Position
#             rj_Type = r[j].Type
#             rtemp =  ri - rj
#             if ri_Type == r[j].Type 
#                 if r[i].Type == :Ga
#                     pe[1] += grad_potential_energy(rtemp, p, ri_Type, rj_Type)
#                 else
#                     pe[2] += grad_potential_energy(rtemp, p, ri_Type, rj_Type)
#                 end
#             else
#                 pe[3] += grad_potential_energy(rtemp, p, ri_Type, rj_Type)
#             end
#         end
#     end
#     return pe
# end

# function grad_force(r::Atom, p::GaN, Type1::Symbol, Type2::Symbol)
    
#     if (Type1 == :Ga) && (Type2 == :Ga) # Ga-Ga interaction
#         return grad_force(r.Position, p.lj_Ga_Ga)
#     elseif (Type1 == :N) && (Type2 == :N) # N-N interaction
#         return grad_force(r.Position, p.lj_N_N)
#     else 
#         return grad_force(r.Position, p.bm_Ga_N)
#     end
# end

# function grad_force(r::Vector{Atom}, p::GaN)
#     # Initialize 
#     rtemp = r[1].Position - r[2].Position
#     f = [ [0. * grad_force(rtemp, p, :Ga, :Ga), 0. * grad_force(rtemp, p, :N, :N), 0. * grad_force(rtemp, p, :Ga, :N)] for j = 1:n]
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             if r[i].Type == r[j].Type 
#                 if r[i].Type == :Ga
#                     f[i][1] += grad_force(rtemp, p, r[i].Type, r[j].Type)
#                     f[j][1] -= f[i][1]
#                 else
#                     f[i][2] += grad_force(rtemp, p, r[i].Type, r[j].Type)
#                     f[j][2] -= f[i][2]
#                 end
#             else
#                 f[i][3] += grad_force(rtemp, p, r[i].Type, r[j].Type)
#                 f[j][3] -= f[i][3]
#             end
#         end
#     end
#     return f
# end

# function grad_virial(r::Atom, p::GaN, Type1::Symbol, Type2::Symbol)
#     df = grad_force(r, p, Type1, Type2)
#     if Type1 != Type2
#         dfdA = df[:dfdA]
#         dfdrho = df[:dfdρ]
#         return (dvdA = dfdA[1]*r.Position[1] + dfdA[2]*r.Position[2] + dfdA[3]*r.Position[3], 
#             dvdρ = dfdrho[1]*r.Position[1] + dfdrho[2]*r.Position[2] + dfdrho[3]*r.Position[3])
#     else
#         dfde = df[:dfdϵ]
#         dfdsig = df[:dfdσ]
#     return (dvdϵ = dfde[1]*r.Position[1] + dfde[2]*r.Position[2] + dfde[3]*r.Position[3], 
#             dvdσ = dfdsig[1]*r.Position[1] + dfdsig[2]*r.Position[2] + dfdsig[3]*r.Position[3])
#     end
# end

# function grad_virial(r::Vector{Atom}, p::GaN)
#     n = length(r)
#     # Initialize 
#     rtemp = r[1].Position - r[2].Position
#     pe = [0. * grad_virial(rtemp, p, :Ga, :Ga), 0. * grad_virial(rtemp, p, :N, :N), 0. * grad_virial(rtemp, p, :Ga, :N)]
#     for i = 1:(n-1)
#         for j = (i+1):n
#             rtemp = r[i] - r[j]
#             if r[i].Type == r[j].Type 
#                 if r[i].Type == :Ga
#                     pe[1] += grad_virial(rtemp, p, r[i].Type, r[j].Type)
#                 else
#                     pe[2] += grad_virial(rtemp, p, r[i].Type, r[j].Type)
#                 end
#             else
#                 pe[3] += grad_virial(rtemp, p, r[i].Type, r[j].Type)
#             end
#         end
#     end
#     return pe
# end
