# InteratomicPotentials

<!--
TODO (Dallas): put README type details here
example: https://cesmix-mit.github.io/Atomistic.jl/dev/
-->

## Conventions

The unit convention throughout the package and other packages in the CESMIX Julia Suite is to assume all unspecified units to be atomic units as defined in the [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl) package. All exposed interfaces should allow for numeric or unitful input. For clarity's sake, it is _strongly recommended_ that user code utilize Unitful wherever possible. Internally, InteratomicPotentials.jl will automatically convert these quantities to be compatible without a significant performance penalty.
