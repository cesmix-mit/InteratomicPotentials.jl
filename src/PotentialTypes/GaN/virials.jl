#####################################################################
####################### Mixed Potentials (GaN) ######################
#####################################################################
function virial(r::Atom, p::GaN, Type1::Symbol, Type2::Symbol)
    f = force(r, p, Type1, Type2)
    return f[1] * r.Position[1] + f[2] * r.Position[2] + f[3] * r.Position[3]
end


function virial(r::Vector{Atom}, p::MixedPotential)
    n = length(r)
    v = 0.0
    for i = 1:(n-1)
        ri = r[i].Position
        for j = (i+1):n
            rj = r[j].Position
            rtemp = ri - rj
            v +=  virial(rtemp, p, r[i].Type, r[j].Type)
        end
    end
    return v
end

function virial(c::Configuration, p::MixedPotential)
    return virial(c.Atoms, p)
end

function virial(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    v = zeros(n)
    for i = 1:n
        v[i] = virial(r[i].Atoms, p)
    end
    return v
end

function virial_stress(r::Atom, p::GaN, Type1::Symbol, Type2::Symbol)
    f = force(r, p, Type1, Type2)
    vi = r.Position * f'
    return [vi[1, 1], vi[2, 2], vi[3, 3], vi[2, 3], vi[1, 3], vi[1,2]]
end

function virial_stress(r::Vector{Atom}, p::MixedPotential)
    n = length(r)
    v = zeros(6)
    for i = 1:(n-1)
        ri = r[i].Position
        for j = (i+1):n
            rj = r[j].Position
            rtemp = ri - rj
            v +=  virial_stress(rtemp, p, r[i].Type, r[j].Type)
        end
    end
    return v
end

function virial_stress(c::Configuration, p::MixedPotential)
    return virial_stress(c.Atoms, p)
end

function virial_stress(r::Vector{Configuration}, p::MixedPotential)
    n = length(r)
    v = [zeros(6) for i = 1:n]
    for i = 1:n
        v[i] = virial_stress(r[i].Atoms, p)
    end
    return v
end






