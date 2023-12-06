export DFTKPotential

"""
    DFTKPotential

Integration to use DFTK as an interatomic potential. Make sure to issue
a `using DFTK` to ensure the `energy_and_force` method using this type
is available.
"""
Base.@kwdef struct DFTKPotential <: NonTrainablePotential
    functionals::Vector{Symbol}    = [:gga_x_pbe, :gga_c_pbe]  # default to model_PBE
    model_kwargs::Dict{Symbol,Any} = Dict{Symbol,Any}()
    basis_kwargs::Dict{Symbol,Any} = Dict{Symbol,Any}()
    scf_kwargs::Dict{Symbol,Any}   = Dict{Symbol,Any}()
end

"""
    DFTKPotential(Ecut, kgrid; kwargs...)

Construct a DFTK potential. `Ecut` is the kinetic energy cutoff,
`kgrid` the k-point grid and other `kwargs` are directly passed
to the `DFTKPotential` constructor. See
[https://docs.dftk.org](https://docs.dftk.org) for more details
on using DFTK.
"""
function DFTKPotential(Ecut, kgrid; kwargs...)
    p = DFTKPotential(; kwargs...)
    p.basis_kwargs[:Ecut]  = Ecut
    p.basis_kwargs[:kgrid] = kgrid
    p
end
