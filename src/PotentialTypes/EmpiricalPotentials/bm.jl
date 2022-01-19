##############################   Born-Mayer  ###################################
struct BornMayer <: EmpiricalPotential
    A::Float64
    ρ::Float64
    rcutoff :: AbstractFloat
end

get_trainable_params(bm::BornMayer) = Parameter{:A,:ρ}((bm.A, bm.ρ))

get_nontrainable_params(bm::BornMayer) = Parameter{:rcutoff}((bm.rcutoff))


# ##############################   Energy  ###################################

function potential_energy(R::AbstractFloat, bm::BornMayer)
    bm.A * exp(-R / bm.ρ)
end


# ##############################   Force   ###################################

function force(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, bm::BornMayer)
    SVector(bm.A /bm.ρ * exp(-R / bm.ρ ) .* r ./ R)
end

# ##############################   Gradients  ###################################

function grad_potential_energy(r::Vector{<:Real}, bm::BornMayer)
    d = norm(r)
    return (dpdA = exp(-d / bm.ρ),
            dpdρ = bm.A * d * exp(-d/bm.ρ) / bm.ρ^2 )
end

function grad_force(r::Vector{<:Real}, bm::BornMayer)
    d = norm(r)
    return (dfdA = 1.0 /bm.ρ * exp(-d / bm.ρ ) .* r ./ d, 
            dfdρ = bm.A / bm.ρ^3 * exp(-d / bm.ρ )*(d - bm.ρ)  .* r ./ d)
end

function grad_virial(r::Vector{<:Real}, bm::BornMayer)
    df = grad_force(r, bm)
    dfdA = df[:dfdA]
    dfdrho = df[:dfdρ]
    return (dvdA = dfdA[1]*r[1] + dfdA[2]*r[2] + dfdA[3]*r[3], 
            dvdρ = dfdrho[1]*r[1] + dfdrho[2]*r[2] + dfdrho[3]*r[3])
end
