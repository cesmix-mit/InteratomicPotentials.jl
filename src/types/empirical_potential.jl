################################################################################
# Types of Empirical Potentials
################################################################################

include("../empirical/lj.jl")
include("../empirical/bm.jl")
include("../empirical/coulomb.jl")
include("../empirical/zbl.jl")
include("../empirical/morse.jl")

export LennardJones, BornMayer, Coulomb, ZBL, Morse

################################################################################
# InteratomicPotentials API implementations for empirical potentials
################################################################################

potential_energy(r::SVector{3}, p::EmpiricalPotential) = potential_energy(norm(r), p)
force(r::SVector{3}, p::EmpiricalPotential) = force(norm(r), r, p)

function _perform_pairwise(f::Function, s::AbstractSystem, p::EmpiricalPotential, data::T) where {T<:Tuple}
    nnlist = neighborlist(s, get_rcutoff(p))

    symbols = atomic_symbol(s)::Vector{Symbol}
    all_data = [Tuple(copy.(data)) for _ ∈ 1:nthreads()]

    if (ismissing(get_species(p)) || unique(symbols) ⊆ get_species(p))
        # Don't need to check each pair for species
        @inbounds @threads for i in 1:length(nnlist)
            for (j, r) in zip(nnlist.j[i], nnlist.r[i])
                f(i, j, r, all_data[threadid()])
            end
        end
    else
        # Need to check each pair for species
        @inbounds @threads for i in 1:length(nnlist)
            if (symbols[i] ∈ get_species(p))
                for (j, r) in zip(nnlist.j[i], nnlist.r[i])
                    if (symbols[j] ∈ get_species(p))
                        f(i, j, r, all_data[threadid()])
                    end
                end
            end
        end
    end

    reduce((x, y) -> x .+ y, all_data)
end

function energy_and_force(s::AbstractSystem, p::EmpiricalPotential)
    data = (zeros(MVector{1,Float64}), zeros(MVector{3,Float64}, length(s)))
    e, f = _perform_pairwise(s, p, data) do i::Int, j::Int, r::SVector{3,Float64}, (e, f)::Tuple{MVector{1,Float64},Vector{MVector{3,Float64}}}
        e .+= potential_energy(r, p)::Float64
        fo = force(r, p)::SVector{3,Float64}
        f[i] += fo
        f[j] -= fo
    end
    (; e=e[] * ENERGY_UNIT, f=SVector{3}.(f) * FORCE_UNIT)
end

function virial_stress(s::AbstractSystem, p::EmpiricalPotential)
    data = (zero(MVector{6,Float64}),)
    v, = _perform_pairwise(s, p, data) do i::Int, j::Int, r::SVector{3,Float64}, (v,)::Tuple{MVector{6,Float64}}
        vi = r * (force(r, p)::SVector{3,Float64})'
        v .+= @SVector [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]]
    end
    SVector{6}(v) * ENERGY_UNIT
end
