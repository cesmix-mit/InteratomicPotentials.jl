# The unit convention throughout the package and other packages in the CESMIX Julia Suite
# is to assume all unspecified units to be atomic units as defined in UnitfulAtomic.jl.

# Here we provide constants for the atomic units we use to make the code more readable.

const ENERGY_UNIT = u"hartree"
const FORCE_UNIT = u"hartree / bohr"

const ENERGY_TYPE = typeof(zeros() * ENERGY_UNIT)
const FORCE_TYPE = typeof(zeros() * FORCE_UNIT)
