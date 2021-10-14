function get_rdf(file::String)
    # Return rdf 
    lines = readlines(file)
    n = length(lines) - 4
    bins = zeros(n)
    rdf = zeros(n)
    coord = zeros(n)
    for j = 1:n
        v = split(lines[4+j], " ")
        bins[j] = parse(Float64, v[2])
        rdf[j] = parse(Float64, v[3])
        coord[j] = parse(Float64, v[4])
    end

    return (bins = bins, rdf = rdf, coord = coord)
end

function get_pe(file::String, dt)
    lines = readlines(file)
    n = length(lines) - 2
    t = zeros(n)
    pe = zeros(n)
    for j = 1:n
        v = split(lines[2+j], " ")
        t[j] = dt*parse(Float64, v[1])
        pe[j] = parse(Float64, v[2])
    end
    return (t = t, pe = pe)
end


function get_msd(file::String, dim::Int)
    lines = readlines(file)
    n = length(lines) - 2
    t_msd = zeros(n)
    msd = zeros(n)
    for j = 1:n
        v = split(lines[2+j], " ")
        t_msd[j] = parse(Float64, v[1])
        msd[j] = parse(Float64, v[2])
    end
    check_n = Int( n*0.9 )
    diffusivity = sum(msd[check_n:end] ./ t_msd[check_n:end]) / (n-check_n)
    if dim == 2
        diffusivity *= 0.25
    else
        diffusivity *= 1.0 / 6.0
    end
    return (t = t_msd, msd = msd, diffusivity = diffusivity)
end

function get_positions(c0::Configuration, path::String, T)
    num_configs = length(T)
    r = Vector{Configuration}(undef, num_configs)
    for (i,t) in enumerate(T)
        file = "DATA.$t"
        c_temp = load_lammps(path*file; atom_names = c0.atom_names, 
                        radii = c0.radii, weights = c0.weights, boundary_type = c0.domain.bound_type)
        r[i] = c_temp
    end

    return r
end

# function animate_atoms(r::Vector{Configuration}, T, path::String; dim = 2, fps = 60)
#     c0 = r[1]
#     n = length(r)
#     if dim == 3
#         anim = @animate for i = 1:n
#             pos = hcat(vec(r[i].Positions)...)
#             scatter(pos[1, :], pos[2, :], pos[3, :], xaxis = ("x", (c0.x_bounds[1], c0.x_bounds[2]), -c0.x_bounds[1]:0.5:c0.x_bounds[2]), 
#                             yaxis = ("y", (c0.y_bounds[1], c0.y_bounds[2]), -c0.y_bounds[1]:0.5:c0.y_bounds[2]), 
#                             zaxis = ("z", (c0.z_bounds[1], c0.z_bounds[2]), -c0.z_bounds[1]:0.5:c0.z_bounds[2]), 
#                             markersize = 20, 
#                             c=colormap("Blues", c0.num_atoms), 
#                             legend = false,
#                             title = string("τ = ", floor(T[i]) ))
#         end
    
#     elseif dim==2
#         anim = @animate for i = 1:n
#             pos = hcat(vec(r[i].Positions)...)
#             scatter(pos[1, :], pos[2, :], xaxis = ("x", (c0.x_bounds[1], c0.x_bounds[2]), -c0.x_bounds[1]:0.5:c0.x_bounds[2]), 
#                             yaxis = ("y", (c0.y_bounds[1], c0.y_bounds[2]), -c0.y_bounds[1]:0.5:c0.y_bounds[2]), 
#                             c=colormap("Blues", c0.num_atoms),
#                             legend = false, 
#                             title = string("τ = ", floor(T[i]) ) )
            
#         end
#     end
    
#     gif(anim, path*"movie.gif", fps = fps)
# end

function run_md(c0::Configuration, lj::LennardJones, Tend::Int, save_dir::String; dim = 3, seed = 1, Temp = 0.5, dt = 0.005, dT = 10)
    d = LMP(["-screen", "none"]) do lmp
        for i = 1
            command(lmp, "log none")
            command(lmp, "units lj")
            command(lmp, "dimension $dim")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            if dim == 2
                command(lmp, "boundary $(c0.domain.bound_type[1]) $(c0.domain.bound_type[2]) p")
            elseif dim == 3
                command(lmp, "boundary $(c0.domain.bound_type[1]) $(c0.domain.bound_type[2]) $(c0.domain.bound_type[3])")
            end

            # Setup box
            command(lmp, "region mybox block $(c0.domain.bounds[1][1]) $(c0.domain.bounds[1][2]) $(c0.domain.bounds[2][1]) $(c0.domain.bounds[2][2]) $(c0.domain.bounds[3][1]) $(c0.domain.bounds[3][2])")
            command(lmp, "create_box $(c0.num_atom_types) mybox")

            # Create atoms
            for j = 1:length(c0.Atoms)
                atom_id = findall(c0.atom_names .== c0.Atoms[j].Type)[1]
                command(lmp, "create_atoms $atom_id single $(c0.Atoms[j].Position[1]) $(c0.Atoms[j].Position[2]) $(c0.Atoms[j].Position[3])")
            end

            command(lmp, "mass 1 1")
            command(lmp, "velocity all create $(Temp) $seed mom yes rot yes dist gaussian")

            # Setup Forcefield
            command(lmp, "pair_style hybrid lj/cut $(maximum(c0.radii))")
            command(lmp, "pair_coeff 1 1 lj/cut $(lj.ϵ) $(lj.σ)")
            
            # Apply fixes
            if dim == 2
                command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed enforce2d")
            elseif dim == 3
                command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed")
            end
            
            # computes
            command(lmp, "compute rdf all rdf 100")      # radial distribution function 
            command(lmp, "fix frdf all ave/time $dT $(Int(Tend/dT)-1) $Tend c_rdf[*] file $(save_dir)tmp.rdf mode vector")
            
            command(lmp, "compute pe all pe")                       # potential energy
            command(lmp, "fix fpe all ave/time 1 100 100 c_pe[*] file $(save_dir)tmp.pe")
            
            command(lmp, "compute msd all msd com yes")
            command(lmp, "fix fmsd all ave/time 1 100 100 c_msd[4] file $(save_dir)tmp.msd")
            
            # Outputs
            command(lmp, "timestep $dt")
            command(lmp, """run $Tend every $dT "write_data $(save_dir)DATA.*" """)
            
            command(lmp, "clear")
        end
    end

    # Return rdf 
    rdf = get_rdf("$(save_dir)tmp.rdf")

    # Return pe 
    pe = get_pe("$(save_dir)tmp.pe", dt)

    # Return msd, diffusivity
    msd = get_msd("$(save_dir)tmp.msd", dim)

    T = dT:dT:Tend
    r = get_positions(c0, save_dir, T)

    return r, rdf, pe, msd

end

function run_md(c0::Configuration, snap::SNAP, Tend::Int, save_dir::String; dim = 3, seed = 1, Temp = 0.5, dt = 0.005, dT = 10)
    create_snap_files(c0, snap, save_dir*"tmp")
    d = LMP(["-screen", "none"]) do lmp
        for i = 1
            command(lmp, "log none")
            command(lmp, "units" * " " * c0.units)
            command(lmp, "dimension $dim")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            if dim == 2
                command(lmp, "boundary $(c0.domain.bound_type[1]) $(c0.domain.bound_type[2]) p")
            elseif dim == 3
                command(lmp, "boundary $(c0.domain.bound_type[1]) $(c0.domain.bound_type[2]) $(c0.domain.bound_type[3])")
            end

            # Setup box
            command(lmp, "region mybox block $(c0.domain.bounds[1][1]) $(c0.domain.bounds[1][2]) $(c0.domain.bounds[2][1]) $(c0.domain.bounds[2][2]) $(c0.domain.bounds[3][1]) $(c0.domain.bounds[3][2])")
            command(lmp, "create_box $(c0.num_atom_types) mybox")


            # Create atoms
            for j = 1:length(c0.Atoms)
                atom_id = findall(c0.atom_names .== c0.Atoms[j].Type)[1]
                command(lmp, "create_atoms $atom_id single $(c0.Atoms[j].Position[1]) $(c0.Atoms[j].Position[2]) $(c0.Atoms[j].Position[3])")
            end

            command(lmp, "mass 1 1")
            command(lmp, "velocity all create $(Temp) $seed mom yes rot yes dist gaussian")

            # Setup Forcefield
            
            cutoff = snap.rcutfac * maximum(c0.radii)
            max_bounds = max(2*cutoff, min( (c0.domain.bounds[1][2] - c0.domain.bounds[1][1]), (c0.domain.bounds[2][2] - c0.domain.bounds[2][1]), (c0.domain.bounds[3][2] - c0.domain.bounds[3][1]) ) )
            # command(lmp, "pair_style hybrid/overlay zero $max_bounds snap")
            # command(lmp, "pair_coeff * * zero")
            command(lmp, "pair_style hybrid/overlay zero $max_bounds zbl 2.0 2.5 snap")
            command(lmp, "pair_coeff * * zero")
            command(lmp, "pair_coeff * * zbl 1 1")
            snap_coeff_path = save_dir*"tmp.snapcoeff"
            snap_param_path = save_dir*"tmp.snapparam"
            command(lmp, "pair_coeff * * snap $(snap_coeff_path) $(snap_param_path) $(c0.atom_names[1])")
            command(lmp, "comm_modify cutoff 4.0")
            # Apply fixes
            # command(lmp, "fix l0 all langevin $Temp $Temp 0.001")
            if dim == 2
                command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed enforce2d")
            elseif dim == 3
                command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed")
            end
            
            # computes
            command(lmp, "compute rdf all rdf 100 cutoff 2.5")      # radial distribution function 
            command(lmp, "fix frdf all ave/time $dT $(Int(Tend/dT)-1) $Tend c_rdf[*] file $(save_dir)tmp.rdf mode vector")
            
            command(lmp, "compute pe all pe")                       # potential energy
            command(lmp, "fix fpe all ave/time 1 100 100 c_pe[*] file $(save_dir)tmp.pe")
            
            command(lmp, "compute msd all msd com yes")
            command(lmp, "fix fmsd all ave/time 1 100 100 c_msd[4] file $(save_dir)tmp.msd")
            
            # Outputs
            command(lmp, "timestep $dt")
            command(lmp, """run $Tend every $dT "write_data $(save_dir)DATA.*" """)
            
            command(lmp, "clear")
        end
    end

    # Return rdf 
    rdf = get_rdf("$(save_dir)tmp.rdf")

    # Return pe 
    pe = get_pe("$(save_dir)tmp.pe", dt)

    # Return msd, diffusivity
    msd = get_msd("$(save_dir)tmp.msd", dim)

    T = dT:dT:Tend
    r = get_positions(c0, save_dir, T)

    return r, rdf, pe, msd

end