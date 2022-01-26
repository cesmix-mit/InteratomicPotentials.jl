using AtomsBase
using Unitful
using UnitfulAtomic
using StaticArrays

function load_data(;num_entries = 2000, file = "LJCluster/curated_lj_cluster.xyz" )
    systems  = Vector{AbstractSystem}(undef, num_entries)
    energies = Vector{Float64}(undef, num_entries)
    forces    = Vector{Vector{T} where T<:SVector{3, <:Float64}}(undef, num_entries)
    bias = @SVector [4.0, 4.0, 4.0]
    box = [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 9.0]] * 1u"Å"
    bcs = [DirichletZero(), DirichletZero(), DirichletZero()]
    open(file, "r") do io
        count = 1
        while !eof(io) && (count <= num_entries)
            line = readline(io)
            num_atoms = parse(Int, line)
            info_line = split(readline(io))
            energies[count] = parse(Float64, info_line[2][8:end])

            atoms = Vector{Atom}(undef, 13)
            force = Vector{SVector{3, Float64}}(undef, 13) 
            for i = 1:num_atoms
                line = split(readline(io))
                element = Symbol(line[1])
                position = @SVector [ parse(Float64, line[2]), parse(Float64, line[3]), parse(Float64, line[4]) ]
                position += bias
                atoms[i] = Atom(element, position * 1u"Å") 

                force[i] = @SVector [ parse(Float64, line[5]), parse(Float64, line[6]), parse(Float64, line[7])]
            end

            forces[count] = force
            systems[count] = FlexibleSystem(atoms, box, bcs)
            count += 1
        end
        println(count)
    end
    return systems, energies, forces 
end

