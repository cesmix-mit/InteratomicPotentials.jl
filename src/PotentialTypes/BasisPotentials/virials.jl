#############################################################################
############################## Simple Virial ################################
#############################################################################
function virial_stress(A::AbstractSystem, p::BasisPotential)
    v = zeros(6)
    d = evaluate_basis_v(A, p.basis_params)
    for di in d 
        v += dot(di, p.coefficients)
    end
    return SVector{6}(v)
end

function virial(A::AbstractSystem, p::BasisPotential)
    v = virial_stress(A, p)
    return v[1] + v[2] + v[3]
end



