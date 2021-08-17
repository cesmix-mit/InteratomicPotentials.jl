#######################################################################
########################### SNAP ######################################
#######################################################################

function virial(c:: Configuration, p::SNAP)
    A = get_snap(c, p)
    virial = A[end-5:end, :] * p.β
    return sum(virial[1:3])
end

function virial(r:: Vector{Configuration}, p::SNAP)
    n = length(r)
    v = Vector{Real}(undef, n)
    for j = 1:n
        v[j] = virial(r[j], p)
    end
    return v
end

function virial_stress(c:: Configuration, p::SNAP)
    A = get_snap(c, p)
    virial = A[end-5:end, :] * p.β
    return virial
end

function virial_stress(r:: Vector{Configuration}, p::SNAP)
    n = length(r)
    v = [Vector{Real}(undef, 6) for i = 1:n]
    for i = 1:n
        v[i] = virial_stress(r[i], p)
    end
    return v
end
