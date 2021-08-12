function run_md(c0::Configuration, lj::LennardJones, Tend::Int, save_dir::String; dim = 3, seed = 1, Temp = 0.5, dt = 0.005, dT = 10)
    d = LMP(["-screen", "none"]) do lmp
        for i = 1
            command(lmp, "log none")
            command(lmp, "units lj")
            command(lmp, "dimension $dim")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            if dim == 2
                command(lmp, "boundary $(c0.boundaries[1]) $(c0.boundaries[2]) p")
            elseif dim == 3
                command(lmp, "boundary $(c0.boundaries[1]) $(c0.boundaries[2]) $(c0.boundaries[3])")
            end

            # Setup box
            command(lmp, "region mybox block $(c0.x_bounds[1]) $(c0.x_bounds[2]) $(c0.y_bounds[1]) $(c0.y_bounds[2]) $(c0.z_bounds[1]) $(c0.z_bounds[2])")
            command(lmp, "create_box $(c0.num_atom_types) mybox")

            # Create atoms
            for j = 1:c0.num_atoms
                atom_id = findall(c0.atom_names .== c0.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c0.Positions[j].x) $(c0.Positions[j].y) $(c0.Positions[j].z)")
            end

            command(lmp, "mass 1 1")
            command(lmp, "velocity all create $(Temp) $seed mom yes rot yes dist gaussian")

            # Setup Forcefield
            command(lmp, "pair_style hybrid lj/cut $(maximum(c0.radii))")
            command(lmp, "pair_coeff 1 1 lj/cut $(lj.ϵ) $(lj.σ)")
            
            # Apply fixes
            # command(lmp, "fix l0 all langevin $Temp $Temp 0.001")
            command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed")
            
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
    lines = readlines("$(save_dir)tmp.rdf")
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

    # Return pe 
    lines = readlines("$(save_dir)tmp.pe")
    n = length(lines) - 2
    t = zeros(n)
    pe = zeros(n)
    for j = 1:n
        v = split(lines[2+j], " ")
        t[j] = dt*parse(Float64, v[1])
        pe[j] = parse(Float64, v[2])
    end

    # Return vacf, diffusivity
    lines = readlines("$(save_dir)tmp.msd")
    n = length(lines) - 2
    t_msd = zeros(n)
    msd = zeros(n)
    # diffusivity = 0.0
    for j = 1:n
        v = split(lines[2+j], " ")
        t_msd[j] = parse(Float64, v[1])
        msd[j] = parse(Float64, v[2])
    end
    diffusivity = sum(msd[end-100:end] ./ t_msd[end-100:end]) / 100.0
    if dim == 2
        diffusivity *= 0.25
    else
        diffusivity *= 1.0 / 6.0
    end

    return (bins = bins, rdf = rdf, coord = coord), (t = t, pe = pe), (t = t_msd, msd = msd, diffusivity = diffusivity)

end

function run_md(c0::Configuration, snap::SNAP, Tend::Int, save_dir::String; dim = 3, seed = 1, Temp = 0.5, dt = 0.005, dT = 10)
    create_snap_files(c0, snap, save_dir*"tmp")
    d = LMP(["-screen", "none"]) do lmp
        for i = 1
            command(lmp, "log none")
            command(lmp, "units lj")
            command(lmp, "dimension $dim")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            if dim == 2
                command(lmp, "boundary $(c0.boundaries[1]) $(c0.boundaries[2]) p")
            elseif dim == 3
                command(lmp, "boundary $(c0.boundaries[1]) $(c0.boundaries[2]) $(c0.boundaries[3])")
            end

            # Setup box
            command(lmp, "region mybox block $(c0.x_bounds[1]) $(c0.x_bounds[2]) $(c0.y_bounds[1]) $(c0.y_bounds[2]) $(c0.z_bounds[1]) $(c0.z_bounds[2])")
            command(lmp, "create_box $(c0.num_atom_types) mybox")

            # Create atoms
            for j = 1:c0.num_atoms
                atom_id = findall(c0.atom_names .== c0.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c0.Positions[j].x) $(c0.Positions[j].y) $(c0.Positions[j].z)")
            end

            command(lmp, "mass 1 1")
            command(lmp, "velocity all create $(0.25*Temp) $seed mom yes rot yes dist gaussian")

            # Setup Forcefield
            command(lmp, "pair_style snap")
            snap_coeff_path = save_dir*"tmp.snapcoeff"
            snap_param_path = save_dir*"tmp.snapparam"
            command(lmp, "pair_coeff * * $(snap_coeff_path) $(snap_param_path) $(c0.atom_names[1])")

            # Apply fixes
            # command(lmp, "fix l0 all langevin $Temp $Temp 0.001")
            command(lmp, "fix l0 all nve langevin $Temp $Temp 1.0 $seed")
            
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
    lines = readlines("$(save_dir)tmp.rdf")
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

    # Return pe 
    lines = readlines("$(save_dir)tmp.pe")
    n = length(lines) - 2
    t = zeros(n)
    pe = zeros(n)
    for j = 1:n
        v = split(lines[2+j], " ")
        t[j] = dt*parse(Float64, v[1])
        pe[j] = parse(Float64, v[2])
    end

    # Return vacf, diffusivity
    lines = readlines("$(save_dir)tmp.msd")
    n = length(lines) - 2
    t_msd = zeros(n)
    msd = zeros(n)
    # diffusivity = 0.0
    for j = 1:n
        v = split(lines[2+j], " ")
        t_msd[j] = parse(Float64, v[1])
        msd[j] = parse(Float64, v[2])
    end
    diffusivity = sum(msd[end-100:end] ./ t_msd[end-100:end]) / 100.0
    if dim == 2
        diffusivity *= 0.25
    else
        diffusivity *= 1.0 / 6.0
    end

    return (bins = bins, rdf = rdf, coord = coord), (t = t, pe = pe), (t = t_msd, msd = msd, diffusivity = diffusivity)

end