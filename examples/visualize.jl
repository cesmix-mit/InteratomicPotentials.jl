import InteratomicPotentials as Potentials
using GLMakie

file_path = joinpath(dirname(pathof(Potentials)), "..", "test", "DATA", "1", "DATA")
radii = 1.5
configuration = Potentials.load_lammps(
    file_path;
    atom_names = [:Ga, :N],
    radii = [radii, radii],
    weights = [1.0, 0.5],
    boundary_type = ["p", "p", "p"],
    units = "metal"
)

import Makie: Plot

# The recipe! - will get called for plot(!)(x::Configuration)
function Makie.plot!(p::Plot(Potentials.Configuration))
    configuration = to_value(p[1]) # first argument is the Configuration

    radii = configuration.masses ./ maximum(configuration.masses)

    for atom_kind in 1:configuration.num_atom_types
        name = configuration.atom_names[atom_kind]
        radius = radii[atom_kind]

        atoms = filter(Atom->Atom.Type == name, configuration.Atoms)
        positions = map(Atom->Point(Atom.Position...), atoms)

        meshscatter!(p, positions, markersize = radius, label=string(name))
    end
end

fig = Figure()
axes = Axis3(fig[1,1], aspect=:data)
plot!(axes, configuration)
