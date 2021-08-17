##############################   Born-Mayer  ###################################
mutable struct BornMayer <: EmpiricalPotential
    A::Float64
    ρ::Float64
end

function BornMayer()
    #ToDO
    return BornMayer(1.0, 1.0)
end

function get_trainable_params(bm::BornMayer)
    return (A = bm.A, ρ = bm.ρ)
end

function get_nontrainable_params(lj::BornMayer)
    p = Parameter{}(())
    return p
end

##############################   Energy  ###################################

function potential_energy(r::Position, p::BornMayer)
    return p.A * exp(-norm(r) / p.ρ)
end


##############################   Force   ###################################

function force(r::Position, p::BornMayer)
    d = norm(r) 
    return p.A /p.ρ * exp(-d / p.ρ ) .* [r.x, r.y, r.z] ./ d
end


