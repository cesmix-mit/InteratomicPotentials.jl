# # Visualizing configurations with Makie.jl
import InteratomicPotentials as Potentials

# Let's load as test data a `GaN` configuration
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

# In order to use Makie we define a recipe
import Makie

# The recipe! - will get called for plot(!)(x::Configuration)
function Makie.plot!(p::Makie.Plot(Potentials.Configuration))
    configuration = to_value(p[1]) # first argument is the Configuration

    colors = [:red, :blue]
    radii = configuration.masses ./ maximum(configuration.masses)

    for atom_kind in 1:configuration.num_atom_types
        name = configuration.atom_names[atom_kind]
        radius = radii[atom_kind]

        atoms = filter(Atom->Atom.Type == name, configuration.Atoms)
        positions = map(Atom->Point(Atom.Position...), atoms)

        meshscatter!(p, positions, markersize = radius, label=string(name), color=colors[atom_kind])
    end
end

# Now we can plot our configuration

using GLMakie
fig = Figure()
axes = Axis3(fig[1,1], aspect=:data)
plot!(axes, configuration)
fig

# Makie allows us to quickly create recordings

function load_configuration(i)
    file_path = joinpath(dirname(pathof(Potentials)), "..", "test", "DATA", string(i), "DATA")
    radii = 1.5
    return Potentials.load_lammps(
        file_path;
        atom_names = [:Ga, :N],
        radii = [radii, radii],
        weights = [1.0, 0.5],
        boundary_type = ["p", "p", "p"],
        units = "metal"
    )
end

# Create a idx `Node`

idx = Node(1)
configuration = @lift load_configuration($idx)

fig = Figure()
axes = Axis3(fig[1,1], aspect=:data)
plot!(axes, configuration)

record(fig, "GaN.mp4", 1:61, framerate=1) do frame
    idx[] = frame
end
