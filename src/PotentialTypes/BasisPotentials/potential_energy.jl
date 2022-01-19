function potential_energy(A::AbstractSystem, p::BasisPotential)
    dot(evaluate_basis(A, p.basis_params), p.coefficients)
end