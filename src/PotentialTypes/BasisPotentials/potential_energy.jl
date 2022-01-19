function potential_energy(A::AbstractSystem, p::BasisPotential)
    d = evaluate_basis(A, p.basis_params)
    e = 0.0
    for di in d
        e += dot(di, p.coefficients)
    end
end