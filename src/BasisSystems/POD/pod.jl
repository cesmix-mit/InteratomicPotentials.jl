struct POD <: BasisSystem 
    species::Vector{Symbol}
    rcutoff::Real 
end 

include("lammps_pod.jl")