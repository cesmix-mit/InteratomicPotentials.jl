module InteratomicPotentialsDFTKExt
using AtomsBase
using DFTK
using InteratomicPotentials
using Unitful
using UnitfulAtomic

function InteratomicPotentials.energy_and_force(system::AbstractSystem,
                                                potential::DFTKPotential)
    model = model_DFT(system, potential.functionals; potential.model_kwargs...)
    basis = PlaneWaveBasis(model; potential.basis_kwargs...)

    scfres = self_consistent_field(basis; potential.scf_kwargs...)
    # cache ψ and ρ as starting point for next calculation
    potential.scf_kwargs[:ψ] = scfres.ψ
    potential.scf_kwargs[:ρ] = scfres.ρ

    (; e=scfres.energies.total * u"hartree", f=compute_forces_cart(scfres) * u"hartree/bohr")
end

end
