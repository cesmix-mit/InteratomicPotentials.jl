#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct PotentialStruct

    style::String # single, pair, triplet, quadlet, snap 
    bondtype::String # nonbonded, bonded 
    bondedatomtypes::Vector{Int64} # a list of atom types for the bond 
    potentialfunction::String # a string for potential function    
    mu::Vector{Float64} # linear potential parameters
    rcut::Float64 # cut-off radius 

end

function addpotential(style::String, bondtype::String, bondedatomtypes::Vector{Int64}, potentialfunction::String, mu::Vector{Float64}, rcut::Float64)

    potential = PotentialStruct(style, bondtype, bondedatomtypes, potentialfunction, mu, rcut)

    return potential
end


function getpotential(potential::PotentialStruct)

    style = potential.style
    bondtype = potential.bondtype
    bondedatomtypes = potential.bondedatomtypes
    potentialfunction = potential.potentialfunction   
    mu = potential.mu
    rcut = potential.rcut        
    
    return style, bondtype, bondedatomtypes, potentialfunction, mu, rcut
end



