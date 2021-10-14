import InteratomicPotentials as Potentials
using Makie

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
