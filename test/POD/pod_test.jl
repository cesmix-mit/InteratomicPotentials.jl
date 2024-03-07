using AtomsIO

sample_hfo2_sys = load_system("./POD/sample_monoclinic_HfO2.xyz")

lmp_pod = LAMMPS_POD("./POD/sample_6body_hfo2_param.pod", [:Hf,:O])

ld = compute_local_descriptors(sample_hfo2_sys,lmp_pod)
@test typeof(ld) <: Vector{Vector{Float64}}
@test size(ld)[1] == 12
@test size(ld[1])[1] == 1122

#= Reference Test
Reference value obtained with LAMMPS compiled w/ Apple clang version 12.0.5 (clang-1205.0.22.11)
targetting macOS (arm64-apple-darwin22.4.0)

using the eapod branch of the cesmix-mit/lammps @ the following commit:
https://github.com/cesmix-mit/lammps/commit/adf9fc11f17efaba7b62ba332b5c5f5c9f579dee
=#

hfo2_lbp = LBasisPotential(lmp_pod, "./POD/sample_6body_2elem_coeffs.pod")
pe = potential_energy(sample_hfo2_sys,hfo2_lbp)
pe ≈ -111.05104842466731441
