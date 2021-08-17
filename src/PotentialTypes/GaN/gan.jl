#########################################################################
##############################   GaN  ###################################
#########################################################################
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

########################################################################
include("energies.jl")
include("forces.jl")
include("virials.jl")