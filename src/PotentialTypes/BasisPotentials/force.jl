function force(A::AbstractSystem, p::BasisPotential)
    d = evaluate_basis_d(A, p.basis_params)
    f = fill(zeros(3), length(A))
    for (i, di) in enumerate(d)
        # show(stdout, "text/plain", di)
        # println("\n")

        # show(stdout, "text/plain", p.coefficients)
        # println("\n")
        # Note that we only take f[i] += ... 
        # This is because the descriptors should already account for f[j]
        f[i] += di' * p.coefficients
    end
    SVector{3}.(f)
end