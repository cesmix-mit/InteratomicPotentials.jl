# #########################################################################
# ##############################   GaN  ###################################
# #########################################################################

# function potential_energy(r::Vector{<:Real}, p::GaN, Type1::Symbol, Type2::Symbol)
#     if (Type1 == :Ga) && (Type2 == :Ga) # Ga-Ga interaction
#         return potential_energy(r, p.c) + potential_energy(r, p.lj_Ga_Ga)
#     elseif (Type1 == :N) && (Type2 == :N) # N-N interaction
#         return potential_energy(r, p.c) + potential_energy(r, p.lj_N_N)
#     else 
#         return potential_energy(r, p.c) + potential_energy(r, p.bm_Ga_N)
#     end
# end

# function potential_energy(r::Vector{Atom}, p::MixedPotential)
#     n = length(r)
#     pe = 0.0
#     for i = 1:(n-1)
#         ri = r[i].Position
#         for j = (i+1):n
#             rj = r[j].Position
#             rtemp = ri - rj
#             pe +=  potential_energy(rtemp, p, r[i].Type, r[j].Type)
#         end
#     end
#     return pe
# end

# function potential_energy(c::Configuration, p::MixedPotential)
#     return potential_energy(c.Atoms, p)
# end

# function potential_energy(r::Vector{Configuration}, p::MixedPotential)
#     n = length(r)
#     pe = zeros(n)
#     for i = 1:n
#         pe[i] = potential_energy(r[i], p)
#     end
#     return pe
# end
