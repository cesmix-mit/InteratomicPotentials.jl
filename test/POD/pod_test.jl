using AtomsIO
#using PotentialLearning 

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
pe_ref = -111.05104842466731441
fref= [ [-0.02746378868026378,   -0.13394934659545488,    0.08657660722828064],
        [ 0.02746406274624429,   -0.13395092253488361,   -0.08652751636541850],
        [ 0.02745734695881662,    0.13395043848863539,   -0.08657510305372776],
        [-0.02746261612019874,    0.13395036752584755,    0.08652643153718376],
        [ 0.09241905860511784,    0.28929482830562214,    0.27397806493641080],
        [-0.09240998240311563,    0.28930056117127434,   -0.27396630756291218],
        [-0.09241824026227247,   -0.28929473367162378,   -0.27397790785541454],
        [ 0.09241063937722754,   -0.28930102483050857,    0.27396588530562915],
        [-0.28422330256274186,    0.18503621409597582,   -0.01150920040065329],
        [ 0.28421014150064283,    0.18502385240845623,    0.01148476431584359],
        [ 0.28422541598839063,   -0.18503778824998007,    0.01150954295384406],
        [-0.28420873514784484,   -0.18502244611335428,   -0.01148526103906378]]

hfo2_lbp = LBasisPotential(lmp_pod, "./POD/sample_6body_2elem_coeffs.pod")
pe = potential_energy(sample_hfo2_sys,hfo2_lbp)
@test pe ≈ pe_ref

f = force(sample_hfo2_sys, hfo2_lbp)
@test f ≈ fref


#= Reference Stress Test
During preliminary testing, many isues arised where the first configuration would produce correct
forces, but later configurations would result in junk forces. Likewise, there were issues that
only appeared with large configurations, or orthogonal systems, etc. 
TODO: wrapped atom example, gas example beyond cutoff

Reference value obtained with LAMMPS compiled w/ Apple clang version 12.0.5 (clang-1205.0.22.11)
targetting macOS (arm64-apple-darwin22.4.0)

using the eapod branch of the cesmix-mit/lammps @ the following commit:
https://github.com/cesmix-mit/lammps/commit/adf9fc11f17efaba7b62ba332b5c5f5c9f579dee
=#

ref_configs = load_trajectory("./POD/stress_test_hfo2_configs.xyz")
@testset "LAMMPS POD Stress Test, First Order" begin 
    lmp_pod2 = LAMMPS_POD("./POD/sample_6body_hfo2_param.pod", [:Hf,:O])
    hfo2_lbp2 = LBasisPotential(lmp_pod2, "./POD/sample_6body_2elem_coeffs.pod")

    @testset for config_num in eachindex(ref_configs)
        config = ref_configs[config_num]
        check_energy = potential_energy(config,hfo2_lbp2)
        @test check_energy ≈ config.system_data.energy
        check_forces = force(config, hfo2_lbp2)
        @test check_forces ≈ config.atom_data.forces
    end
end 

@testset "LAMMPS POD Stress Test, Reverse Order" begin 
    lmp_pod3 = LAMMPS_POD("./POD/sample_6body_hfo2_param.pod", [:Hf,:O])
    hfo2_lbp3 = LBasisPotential(lmp_pod3, "./POD/sample_6body_2elem_coeffs.pod")

    @testset for config_num in reverse(eachindex(ref_configs))
        config = ref_configs[config_num]
        check_energy = potential_energy(config,hfo2_lbp3)
        @test check_energy ≈ config.system_data.energy
        check_forces = force(config, hfo2_lbp3)
        @test check_forces ≈ config.atom_data.forces
    end
end 