# InteratomicBasisPotentials

Package that provides add-ons to the [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl) package for interatomic potentials that are defined using basis functions, in particular we support the Spectral Neighbor Analysis Potential (SNAP) [Thompson2015](@cite) and the Atomic Cluster Expansion (ACE) [Drautz2019](@cite). Defines an API for abstract interatomic potentials (currently supporting 2-body interactions) and the atomic configuration interface defined by AtomsBase.jl.

Developed as part of the CESMIX Julia package suite. See also ComposableWorkflows, InteratomicPotentials.jl, and PotentialLearning.jl.


## Conventions

The unit convention throughout the package and other packages in the CESMIX Julia Suite is to assume all unspecified units to be atomic units as defined in the [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl) package. All exposed interfaces should allow for numeric or unitful input. For clarity's sake, it is _strongly recommended_ that user code utilize Unitful wherever possible. Internally, InteratomicPotentials.jl will automatically convert these quantities to be compatible without a significant performance penalty.


## Next Steps

If you would like to use the InteratomicPotentials suite in a molecular dynamics simulator, see [Atomistic.jl](https://github.com/cesmix-mit/Atomistic.jl). There, you will learn more about how the abstract classes provided in the present package can be used in conjuction with the Atomistic API and a variety of MD simulators. 

If you would like to fit your potential parameters to data, see our project at [PotentialLearning.jl](https://github.com/cesmix-mit/PotentialLearning.jl), a work in progress, that aims to provide support for a variety of learning and inference tasks.
