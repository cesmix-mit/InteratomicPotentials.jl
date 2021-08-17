#############################################################################
############################## Simple Virial ################################
#############################################################################
function virial(r::Position, p::EmpiricalPotential)
    f = force(r, p)
    return f[1] * r.x + f[2] * r.y + f[3] * r.z
end


function virial_stress(r::Position, p::EmpiricalPotential)
    f = force(r, p)
    vi = [r.x, r.y, r.z] * f'
    return [vi[1, 1], vi[2, 2], vi[3, 3], vi[2, 3], vi[1, 3], vi[1,2]]
end