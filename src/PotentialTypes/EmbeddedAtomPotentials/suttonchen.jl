# Contributions by: Stefan Bringuier
# Reference: A. P. Sutton and J. Chen, Phil. Mag. Lett. 61, 139 (1990)

##############################   Sutton-Chen EAM  ###################################
struct SuttonChen <:EmbeddedAtomPotential
    ϵ
    a
    c
    n
    m
    rcutoff
    species :: AbstractVector
end

get_parameters(sc::SuttonChen) = Parameter{ (:ϵ,:a,:c,:n,:m) }((sc.ϵ,sc.a,sc.c,sc.n,sc.m))
set_parameters(p::Parameter{(:ϵ,:a,:c,:n,:m)}, sc::SuttonChen) = SuttonChen(p.ϵ,p.a,p.c,p.n,p.m,sc.rcutoff,sc.species)

deserialize_parameters(p::Parameter{(:ϵ, :a, :c, :n, :m)}, sc::SuttonChen) = [p.ϵ, p.a, p.c, p.n, p.m]
serialize_parameters(p::Vector, sc::SuttonChen) = Parameter{(:ϵ, :a, :c, :n, :m)}( (p[1], p[2], p[3], p[4], p[5]) )

get_hyperparameters(sc::SuttonChen) = Parameter{(:rcutoff,)}( (sc.rcutoff,) )
set_hyperparameters(p::Parameter{(:rcutoff,)}, sc::SuttonChen) = SuttonChen(sc.ϵ, sc.a, sc.c, sc.n, sc.m, p.rcutoff, sc.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, sc::SuttonChen) = [p.rcutoff]
serialize_hyperparameters(p::Vector, sc::SuttonChen) = Parameter{(:rcutoff, )}( (p[1],) )

# ##############################   Cutoff  ###################################
function fc(R::AbstractFloat, sc::SuttonChen)
    return R-sc.rcutoff < 0.0 ? exp(1/(R-sc.rcutoff)) : 0.0
end

function dfc_dR(R::AbstractFloat, sc::SuttonChen)
    return R-sc.rcutoff < 0.0 ? (-1/(R-sc.rcutoff))*exp(1/(R-sc.rcutoff)) : 0.0
end


# ##############################   Density  ###################################
function rho(R::AbstractFloat, sc::SuttonChen)
    ρ = (sc.a / R)^(sc.m) * fc(R,sc)
    return ρ
end

function drho_dR(R::AbstractFloat,r::SVector{3,<:AbstractFloat},sc::SuttonChen)
    fo = (-sc.m/R) * (sc.a/R)^(sc.m) * fc(R,sc) + (sc.a/R)^(sc.m) * dfc_dR(R,sc)
    return SVector(fo .* r ./ R)
end

# ##############################   Energy  ###################################
function potential_energy_repulsive(R::AbstractFloat, sc::SuttonChen)
    return sc.ϵ * (sc.a/R)^(sc.n)  * fc(R,sc)
end

function potential_energy_embedding(ρ::AbstractFloat,sc::SuttonChen)
    return sc.ϵ * sc.c * sqrt(ρ)
end
# ##############################   Force   ###################################

function force_repulsive(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, sc::SuttonChen)
    fo = sc.ϵ * ((-sc.n/R) * (sc.a/R)^(sc.n) * fc(R,sc) + (sc.a/R)^(sc.n) * dfc_dR(R,sc))
    return SVector(fo .* r ./ R)
end

function force_embedding(ρ::AbstractFloat, dρ::SVector{3,<:AbstractFloat}, sc::SuttonChen)
    return sc.ϵ * (sc.c/(2*sqrt(ρ)) * dρ)
end

# ##############################   Gradients  ###################################
# TO DO
