################################################################################
#
#    This file contains utility functions
#           1. Functions for parameter containers (redefined NamedTuples)
#               1.1 Adding parameter containers
#               1.2 Converting parameter containers into vectors.
#           2. Positions
#           3. Configurations
#               - A Configuration is a type that holds information typically passed to LAMMPS.
#                   - A Configuration is really just a julia wrapper for information contained in DATA files.
#                   - Configuration contains Positions (depending on the number of atoms).
#               - Configurations are particularly important for SNAP implementation, since
#                   SNAP potential, forces, virial, and learning makes use of LAMMPS (and hence should 
#                   use configurations, instead of Positions).
#           4. System (not yet implemented)
#               - A system is a collection of Configurations. (Primarily used for training SNAP potentials)

################################################################################
import Base.+
import Base.-
import Base.*
import Base.copy
import Base.vec

include("param.jl")
include("position.jl")
include("config.jl")





