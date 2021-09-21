#########################################################################
##############################  SNAP  ###################################
#########################################################################
function force(c::Configuration, p::SNAP)
    A = get_snap(c, p)
    force = A[2:end-6, 1:(end-1)] * p.β[1:(end-1)]
    n = length(force)
    num_tuples = Int(n/3)
    force = [[force[3*(j-1)+k] for k = 1:3] for j = 1:num_tuples]
    return force
end

function force(r::Vector{Configuration}, p::SNAP)
    n = length(r)
    f = [[zeros(3) for j = 1:length(r[i].Atoms)] for i = 1:n]
    for i = 1:n 
        f[i] = force(r[i], p)
    end
    return f
end