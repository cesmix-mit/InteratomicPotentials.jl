################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

abstract type ArbitraryPotential end
abstract type Potential <:ArbitraryPotential end
abstract type FittedPotential <:ArbitraryPotential end
abstract type MixedPotential <:ArbitraryPotential end

############################## Lennard Jones ###################################
mutable struct LennardJones <: Potential
    ϵ::Float64
    σ::Float64
end

function LennardJones()
    #ToDO
    return LennardJones(1.0, 1.0)
end

function get_trainable_params(lj::LennardJones)
    p = Parameter{:ϵ, :σ}((lj.ϵ, lj.σ))
    return p
end

function get_nontrainable_params(lj::LennardJones)
    p = Parameter{}(())
    return p
end

##############################   Born-Mayer  ###################################
mutable struct BornMayer <: Potential
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

##############################   Coulomb  ###################################
mutable struct Coulomb <: Potential
    q_1::Float64
    q_2::Float64
    ϵ0::Float64
end

function Coulomb()
    #ToDO
    return Coulomb(1.0, 1.0, 1.0)
end

function get_trainable_params(c::Coulomb)
    p = Parameter{}(())
    return p
end

function get_nontrainable_params(c::Coulomb)
    return (q_1 = c.q_1, q_2 = c.q_2, ϵ0 = c.ϵ0)
end

##############################   GaN  ###################################
mutable struct GaN <: MixedPotential
    lj_Ga_Ga::LennardJones
    lj_N_N::LennardJones
    bm_Ga_N::BornMayer
    c::Coulomb
end

function GaN()
    return GaN(LennardJones(), LennardJones(), BornMayer(), Coulomb())
end

function get_trainable_params(g::GaN)
    p = Parameter{(:Ga_ϵ, :Ga_σ)}(values(get_trainable_params(g.lj_Ga_Ga)))
    p = p ⊕ Parameter{(:N_ϵ, :N_σ)}(values(get_trainable_params(g.lj_Ga_Ga)))
    p = p ⊕ get_trainable_params(g.bm_Ga_N)
    return p
end

function get_nontrainable_params(g::GaN)
    p = get_nontrainable_params(g.c)
    return p
end

##############################  SNAP  ###################################
mutable struct SNAP <: FittedPotential
    β :: Vector{Float64}
    r_cutoff_factor :: Float64 
    twojmax         :: Int
    num_atom_types  :: Int
end

function SNAP(r_cutoff_factor::Float64, twojmax::Int, num_atom_types::Int)
    J = twojmax / num_atom_types
    num_coeffs = Int(floor((J + 1.) * (J + 2.) * (( J + (1.5)) / 3. )))
    return SNAP(ones(num_atom_types*num_coeffs+1), r_cutoff_factor, twojmax, num_atom_types)
end

function get_trainable_params(snap::SNAP)
    return (β = snap.β)
end

function get_nontrainable_params(snap::SNAP)
    return (r_cutoff_factor = snap.r_cutoff_factor, twojmax = snap.twojmax)
end

