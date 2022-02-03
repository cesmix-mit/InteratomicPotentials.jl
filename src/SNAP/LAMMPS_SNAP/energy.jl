#########################################################################
##############################  SNAP  ###################################
#########################################################################
function potential_energy(c::Configuration, p::SNAP)
    A = get_snap(c, p)
    p_e = dot(A[1, :]  ,  p.β )
    return p_e
end

function potential_energy(r::Vector{Configuration}, p::SNAP)
    n = length(r)
    pe = zeros(n)
    for i = 1:n
        A = get_snap(r[i], p)
        pe[i] = dot(A[1, :] , p.β)
    end
    return pe
end
