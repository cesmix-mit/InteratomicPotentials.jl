
using LAMMPS
using Printf: @sprintf 
using LinearAlgebra: norm

# species map can differ in order from species in pod_spec 
# tradeoff between RefValue for the cache vs. just using a mutable struct?
mutable struct LAMMPS_POD <: BasisSystem 
    lmp::LMP
    param_file::String
    species_map::Vector{Symbol} # index corresponds to lammps type
    pod_spec::Union{POD,Nothing}
    num_perelem_ld::Int64
    type_cache::Matrix{Int32} # extra logic for compute dd, this is the output type we receive 
    c_dd_cache::Int64 # extra logic for compute dd
    
    LAMMPS_POD(lmp,param_file,species_map,pod_spec,num_perelem_ld) = new(lmp,param_file,species_map,pod_spec,num_perelem_ld,Int32[-1;;],-1)
end

function LAMMPS_POD(param_file::String, lammps_species::Vector{Symbol}; parse_param_file=false)
    lmp = initialize_pod_lammps(param_file,lammps_species)
    num_perelem_ld = get_num_perelem_ld(lmp,lammps_species)
    if parse_param_file 
        println("Sorry, I haven't implementing parsing POD parameter files yet")
    else
        lmp_pod = LAMMPS_POD(lmp,param_file,lammps_species,nothing,num_perelem_ld)
        #@warn "Until POD param parser implemented, assuming species_map is the same order as species in POD parameter file"
    end
    lmp_pod
end

Base.length(lmp_pod::LAMMPS_POD) = length(lmp_pod.species_map)*lmp_pod.num_perelem_ld

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

    command(lmp, "region main prism 0.0 2.0 0.0 2.0 0.0 2.0 0.0 0.0 0.0") # dummy triclinic system
    command(lmp, "create_box $(num_types) main")
 
    # Does the cutoff matter here? If it does, need to pass in cutoff, which means I'd have to parse from param file
    command(lmp, "pair_style    zero 10.0")
    command(lmp, "pair_coeff    * * ")

    command(lmp, """compute ld all pod/atom $(param_file) "" $(atomtype_str)""")

    lmp
end

function get_num_perelem_ld(lmp::LMP, lammps_species::Vector{Symbol})
    # These masses will get overwritten 
    for i in 1:length(lammps_species)
        command(lmp, "mass $(i) 1.0")
    end
    # need to have one dummy atom for this to work
    command(lmp, "create_atoms 1 single 0.25 0.25 0.25")
    command(lmp, "run 0")

    raw_ld = extract_compute(lmp,"ld", STYLE_ATOM,TYPE_ARRAY)'
    num_perelem_ld = size(raw_ld)[2] + 1 # including 1-body terms

    command(lmp,"delete_atoms group all")

    command(lmp, "pair_style none")
    command(lmp, "pair_style    zero 10.0")
    command(lmp, "pair_coeff    * * ")

    num_perelem_ld
end

function setup_lammps_system!(A::AbstractSystem, pod::LAMMPS_POD)
    # TODO: Will eventually need to ensure unit consistency
    # For now, just trusting that it's correct.

    lmp = pod.lmp
    command(lmp,"delete_atoms group all")

    atom_syms = atomic_symbol(A)
    unq_syms = unique(atom_syms)
    @assert all(in.(unq_syms, (pod.species_map,)))
    atom_types = map(sym->findfirst(isequal(sym),pod.species_map), atom_syms) # atomic sybol --> lammps type based on order of symbol in pod.species_map

    # the way POD is constructed, the number of types initialized should equal the number of species modelled by POD
    # However, some systems will have only a subset of these species, so need to account for that here
    for unq_sym in unq_syms
        atype = findfirst(isequal(unq_sym),pod.species_map)
        mass  = ustrip(atomic_mass(A)[findfirst(isequal(unq_sym),atom_syms)])

        command(lmp, "mass $(atype) $(mass)")
    end

    # Doing it this way avoids the ReinterpretArray type, which is annoying for dispatch on the longest_diag fn
    bbox = [ustrip.(box_arr) for box_arr in bounding_box(A)]
    @assert iszero([bbox[1][2],bbox[1][3],bbox[2][3]]) #needs to conform to lammps triclinic requirements
    bbound_str = ""
    cart_name = ["x", "y", "z"]
    for i in 1:3
        bbound_str = bbound_str * " $(cart_name[i]) final" * (@sprintf " %.27f" 0.0) * (@sprintf " %.27f" bbox[i][i])
    end

    if !iszero(bbox[2][1])
        bbound_str = bbound_str * " xy final" * (@sprintf " %.27f" bbox[2][1])
    else
        bbound_str = bbound_str * " xy final" * (@sprintf " %.27f" 0.0)
    end

    if !iszero(bbox[3][1])
        bbound_str = bbound_str * " xz final" * (@sprintf " %.27f" bbox[3][1])
    else
        bbound_str = bbound_str * " xz final" * (@sprintf " %.27f" 0.0)
    end

    if !iszero(bbox[3][2])
        bbound_str = bbound_str * " yz final" * (@sprintf " %.27f" bbox[3][2])
    else
        bbound_str = bbound_str * " yz final" * (@sprintf " %.27f" 0.0)
    end

    bbound_str = bbound_str * " boundary p p p" 
    command(lmp, "change_box all" * bbound_str)

    max_diag = longest_diagonal_length(bbox[1],bbox[2],bbox[3])
    sort_binsize = 0.5*max_diag + 1.0 # half the max diag length, plus a buffer
    command(lmp, "atom_modify sort 1 $(sort_binsize)")

    # LAMMPS needs wrapped coordinates to create_atoms
    # If this incurs too much a perf penalty, can check to see if needed
    wrapped_pos = wrap_positions(A)
    for i in axes(wrapped_pos,2)
        xyz_str = ""
        for j in 1:3
            xyz_str = xyz_str * @sprintf " %.27f" wrapped_pos[j,i]
        end
        command(lmp, "create_atoms $(atom_types[i]) single" * xyz_str)
    end
end

function compute_local_descriptors(A::AbstractSystem, pod::LAMMPS_POD)
    lmp = pod.lmp
    setup_lammps_system!(A,pod)
    command(lmp, "run 0")

    atomids = extract_atom(lmp, "id", LAMMPS_INT)
    sort_idxs = sortperm(atomids)
    @assert length(A) == length(atomids)

    raw_ld = extract_compute(lmp,"ld", STYLE_ATOM,TYPE_ARRAY)'
    raw_types = extract_atom(lmp,"type", LAMMPS_INT)
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

    command(lmp, "pair_style none")
    command(lmp, "pair_style    zero 10.0")
    command(lmp, "pair_coeff    * * ")

    final_ld
end

function compute_force_descriptors(A::AbstractSystem, pod::LAMMPS_POD)
    lmp = pod.lmp

    atomtype_str = ""
    for elem_symbol in pod.species_map 
        atomtype_str = atomtype_str * " " * string(elem_symbol)
    end

    setup_lammps_system!(A,pod)
    num_atoms = length(A)::Int64
    atomids = extract_atom(lmp, "id", LAMMPS_INT)
    @assert num_atoms == length(atomids)

    sort_idxs = sortperm(atomids)::Vector{Int64}
    @assert sort_idxs == [Int64(i) for i in 1:num_atoms] #so we can use raw_dd
    raw_types = extract_atom(lmp,"type", LAMMPS_INT)::Vector{Int32}
    sorted_types = raw_types[sort_idxs,:]::Array{Int32} # is it Array because I'm slicing it?

    if sorted_types != pod.type_cache
        pod.type_cache = sorted_types # this should be OK because sorted_type is a copy of the lammps types array
        if pod.c_dd_cache == -1
            pod.c_dd_cache = 0
            command(lmp, """compute dd$(pod.c_dd_cache) all pod/global $(pod.param_file) "" $(atomtype_str)""")
        else 
            command(lmp, "uncompute dd$(pod.c_dd_cache)")
            pod.c_dd_cache += 1
            command(lmp, """compute dd$(pod.c_dd_cache) all pod/global $(pod.param_file) "" $(atomtype_str)""")
        end
    end

    command(lmp, "run 0")

    num_pod_types = length(pod.species_map)::Int64
    num_perelem_ld = pod.num_perelem_ld::Int64
    total_num_ld = num_pod_types*(num_perelem_ld)

    raw_dd = extract_compute(lmp, "dd$(pod.c_dd_cache)", STYLE_GLOBAL, TYPE_ARRAY)::Array{Float64,2}
    #raw_dd = extract_compute(lmp, "dd$(pod.c_dd_cache)", STYLE_GLOBAL, TYPE_ARRAY;
    #                         size_2d=(total_num_ld, 3*num_atoms+1))::Array{Float64,2}
    raw_dd = raw_dd'

    final_dd = [[raw_dd[3*i+k+1,1:end] for k in 1:3] for i in 0:num_atoms-1]

    command(lmp, "pair_style none")
    command(lmp, "pair_style    zero 10.0")
    command(lmp, "pair_coeff    * * ")

    final_dd
end


#This should live somewhere else eventually 
function wrap_positions(A::AbstractSystem)
    cart_pos = reduce(hcat,ustrip.(position(A)))

    cry_to_cart = reduce(hcat, ustrip.(bounding_box(A))) # effectively performs transpose of cell matrix
    cart_to_cry = inv(cry_to_cart)

    cry_pos = cart_to_cry*cart_pos
    new_cry_pos =zeros(size(cry_pos))
    for i in axes(new_cry_pos,2)
        new_cry_pos[:,i] = cry_pos[:,i]
        for j in 1:3
            if new_cry_pos[j,i] > 1.0
                new_cry_pos[j,i] -= 1.0
            elseif new_cry_pos[j,i] < 0.0
                new_cry_pos[j,i] += 1.0
            end
        end
    end
    new_cart_pos = cry_to_cart * new_cry_pos
    new_cart_pos
end

#function longest_diagonal_length(a::AbstractArray{T}, b::Vector{T},c::Vector{T}) where T <: Real
function longest_diagonal_length(a::AbstractArray{T}, b::AbstractArray{T},c::AbstractArray{T}) where T <: Real

    # the four possible diagonals in the parallelepiped
    d1 = a + b + c
    d2 = a + b - c
    d3 = a - b + c
    d4 = a - b - c

    max_length = maximum(map(x->norm(x), [d1,d2,d3,d4]))

    max_length
end
