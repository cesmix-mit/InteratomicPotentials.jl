function force(A::AbstractSystem, p::BasisPotential)
    number_of_particles = length(A)
    d = evaluate_basis_d(A, p.basis_params)
    f = [zeros(3) for i = 1:number_of_particles]
    for (i, di) in enumerate(d)
        # Note that we only take f[i] += ... 
        # This is because the descriptors should already account for f[j]
        f[i] += dot(di, p.coefficients)
    end
    return SVector{number_of_particles}([SVector{3}(fi) for fi in f])
end