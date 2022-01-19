##############################   Born-Mayer  ###################################
# mutable struct BornMayer <: EmpiricalPotential
#     A::Float64
#     ρ::Float64
# end

# function BornMayer()
#     #ToDO
#     return BornMayer(1.0, 1.0)
# end

# function get_trainable_params(bm::BornMayer)
#     return (A = bm.A, ρ = bm.ρ)
# end

# function get_nontrainable_params(lj::BornMayer)
#     p = Parameter{}(())
#     return p
# end

# ##############################   Energy  ###################################

# function potential_energy(r::Vector{<:Real}, p::BornMayer)
#     return p.A * exp(-norm(r) / p.ρ)
# end


# ##############################   Force   ###################################

# function force(r::Vector{<:Real}, p::BornMayer)
#     d = norm(r) 
#     return p.A /p.ρ * exp(-d / p.ρ ) .* r ./ d
# end

# ##############################   Gradients  ###################################

# function grad_potential_energy(r::Vector{<:Real}, p::BornMayer)
#     d = norm(r)
#     return (dpdA = exp(-d / p.ρ),
#             dpdρ = p.A * d * exp(-d/p.ρ) / p.ρ^2 )
# end

# function grad_force(r::Vector{<:Real}, p::BornMayer)
#     d = norm(r)
#     return (dfdA = 1.0 /p.ρ * exp(-d / p.ρ ) .* r ./ d, 
#             dfdρ = p.A / p.ρ^3 * exp(-d / p.ρ )*(d - p.ρ)  .* r ./ d)
# end

# function grad_virial(r::Vector{<:Real}, p::BornMayer)
#     df = grad_force(r, p)
#     dfdA = df[:dfdA]
#     dfdrho = df[:dfdρ]
#     return (dvdA = dfdA[1]*r[1] + dfdA[2]*r[2] + dfdA[3]*r[3], 
#             dvdρ = dfdrho[1]*r[1] + dfdrho[2]*r[2] + dfdrho[3]*r[3])
# end
