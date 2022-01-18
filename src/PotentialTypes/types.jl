################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

abstract type ArbitraryPotential end
abstract type EmpiricalPotential <: ArbitraryPotential end
abstract type FittedPotential <: ArbitraryPotential end
abstract type MixedPotential <: ArbitraryPotential end

############################### Empirical Potentials ################################################
include("EmpiricalPotentials/lj.jl")
# include("EmpiricalPotentials/bm.jl")
# include("EmpiricalPotentials/coulomb.jl")
include("EmpiricalPotentials/virials.jl")
include("EmpiricalPotentials/potential_energy.jl")
include("EmpiricalPotentials/force.jl")

################################ SNAP ##############################################################
include("SNAP/snap.jl")

################################ GaN ###############################################################
# include("GaN/gan.jl")
