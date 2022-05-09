# InteratomicBasisPotentials.jl

[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square")](https://mit-license.org)
[![Build Status](https://github.com/cesmix-mit/InteratomicPotentials.jl/workflows/CI/badge.svg)](https://github.com/cesmix-mit/InteratomicBasisPotentials.jl/actions)
[![codecov](https://codecov.io/gh/cesmix-mit/InteratomicBasisPotentials.jl/branch/main/graph/badge.svg?token=IF6zvl50j9)](https://codecov.io/gh/cesmix-mit/InteratomicBasisPotentials.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesmix-mit.github.io/InteratomicBasisPotentials.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesmix-mit.github.io/InteratomicBasisPotentials.jl/dev)


Extension of [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl) that contains machine learning potentials. Right now we support the Spectral Neighbor Analysis Potential (SNAP) potential of Thompson et al. 2015 (see documentation for bibliography) and a naive hacking of the Atomic Cluster Expansion of Drautz 2019 through the [ACE1.jl](https://github.com/ACEsuit/ACE1.jl/) julia package. We are working toward more complete implementation of these machine learning or data-driven potentials in the context of the CESMIX julia suite that seeks to fit and run these potentials for molecular dynamics. For additional details, see the [CESMIX](https://github.com/cesmix-mit) ecosystem.
