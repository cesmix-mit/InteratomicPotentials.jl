
using LAMMPS
using Printf: @sprintf 

# species map can conceivably be separate from species in pod_spec...
struct LAMMPS_POD <: BasisSystem 
    lmp::LMP
    param_file::String
    species_map::Vector{Symbol} # index corresponds to lammps type
    pod_spec::Union{POD,Nothing}
end

function LAMMPS_POD(param_file::String, lammps_species::Vector{Symbol}; parse_param_file=false)
    lmp = initialize_pod_lammps(param_file,lammps_species)
    if parse_param_file 
        println("Sorry, I haven't implementing parsing POD parameter files yet")
    else
        lmp_pod = LAMMPS_POD(lmp,param_file,lammps_species,nothing)
        @warn "Until POD param parser implemented, assuming species_map is the same order as species in POD parameter file"
    end
    lmp_pod
end

function LBasisPotential(lmp_pod::LAMMPS_POD, coeff_fname::String)
    coeffs = vec(readdlm(coeff_fname, ' '; skipstart=1)) 
    LBasisPotential(coeffs,zeros(1),lmp_pod)
end

function initialize_pod_lammps(param_file::String, lammps_species::Vector{Symbol})
    num_types = length(lammps_species)
    atomtype_str = ""
    for elem_symbol in lammps_species
        atomtype_str = atomtype_str * " " * string(elem_symbol)
    end

    lmp = LMP(["-screen", "none"])
    command(lmp, "log none")

    command(lmp, "units         metal")
    command(lmp, "boundary      p p p")
    command(lmp, "atom_style    atomic")
    command(lmp, "neighbor      0.5 bin")
    command(lmp, "neigh_modify  delay 0 every 1 check no")    

    command(lmp, "region main prism 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0")
    command(lmp, "create_box $(num_types) main")

    # Does the cutoff matter here? If it does, need to pass in cutoff, which means I'd have to parse from param file
    command(lmp, "pair_style    zero 10.0")
    command(lmp, "pair_coeff    * * ")

    command(lmp, """compute ld all pod/atom $(param_file) "" "" $(atomtype_str)""")
    command(lmp, """compute dd all podd/atom $(param_file) "" "" $(atomtype_str)""")

    lmp
end

function setup_lammps_system!(A::AbstractSystem, pod::LAMMPS_POD)
     # TODO: Will eventually need to ensure unit consistency
     # For now, just trusting that it's correct.

     lmp = pod.lmp

    atom_syms = atomic_symbol(A)
    unq_syms = unique(atom_syms)
    @assert all(in.(unq_syms, (pod.species_map,)))
    atom_types = map(sym->findfirst(isequal(sym),pod.species_map), atom_syms)

    # the way POD is constructed, the number of types initialized should equal the number of species modelled by POD
    # However, some systems will have only a subset of these species, so need to account for that here
    for unq_sym in unq_syms
        atype = findfirst(isequal(unq_sym),pod.species_map)
        mass  = ustrip(atomic_mass(A)[findfirst(isequal(unq_sym),atom_syms)])

        command(lmp, "mass $(atype) $(mass)")
    end

    bbox = ustrip.(bounding_box(A))
    @assert iszero([bbox[1][2],bbox[1][3],bbox[2][3]]) #needs to conform to lammps triclinic requirements
    bbound_str = ""
    cart_name = ["x", "y", "z"]
    for i in 1:3
        bbound_str = bbound_str * " $(cart_name[i]) final" * (@sprintf " %.27f" 0.0) * (@sprintf " %.27f" bbox[i][i])
    end

    if !iszero(bbox[2][1])
        bbound_str = bbound_str * " xy final" * (@sprintf " %.27f" bbox[2][1])
    end

    if !iszero(bbox[3][1])
        bbound_str = bbound_str * " xz final" * (@sprintf " %.27f" bbox[3][1])
    end

    if !iszero(bbox[3][2])
        bbound_str = bbound_str * " yz final" * (@sprintf " %.27f" bbox[3][2])
    end

    bbound_str = bbound_str * " boundary p p p" 
    command(lmp, "change_box all" * bbound_str)

    atom_pos = ustrip.(position(A))
    for i in 1:length(atom_pos)
        xyz_str = ""
        for j in 1:3
            xyz_str = xyz_str * @sprintf " %.27f" atom_pos[i][j]
        end
        command(lmp, "create_atoms $(atom_types[i]) single" * xyz_str)
    end
end

function compute_local_descriptors(A::AbstractSystem, pod::LAMMPS_POD)
    lmp = pod.lmp
    setup_lammps_system!(A,pod)
    command(lmp, "run 0")

    atomids = extract_atom(lmp, "id")
    sort_idxs = sortperm(atomids)

    raw_ld = extract_compute(lmp,"ld", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)'
    raw_types = extract_atom(lmp,"type")
    sorted_ld = raw_ld[sort_idxs,:]
    sorted_types = raw_types[sort_idxs,:]

    # This block is where the assumption that the species order of pod.species_map should match the POD species order matters
    # Once a param parser is implemented, I won't need this assumption (but will need a bit of extra logic)
    num_pod_types = length(pod.species_map)
    num_atoms = size(sorted_ld)[1]
    num_perelem_ld = size(sorted_ld)[2] + 1 # including 1-body terms
    total_num_ld = num_pod_types*(num_perelem_ld)
    final_ld = [zeros(total_num_ld) for _ in 1:num_atoms] #for consistency, vec{vec} not matrix
    for i in 1:num_atoms
        atype = sorted_types[i] 
        final_ld[i][(atype-1)*num_perelem_ld + 1] = 1.0 # one-body terms
        start = (atype-1)*num_perelem_ld + 2 # +2 accounts for both skipping 1-body term and 1-indexing
        stop  = (atype-1)*num_perelem_ld + num_perelem_ld
        final_ld[i][start:stop] = sorted_ld[i,:]
    end

    command(lmp,"delete_atoms group all")

    final_ld
end
