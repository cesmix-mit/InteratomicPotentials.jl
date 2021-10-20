using InteratomicPotentials 
using StaticArrays

# Test Non-implemented types
@test isa( Angle(), Angle )
@test isa( Bond(), Bond )
@test isa( Dihedral(), Dihedral)
@test isa( Improper(), Improper)

# Test Atom
mass     = 18.0
position = SA_F64[1.0, 1.0, 1.0]
velocity = SA_F64[0.0, 0.0, 0.0]
symbol = :Ar 
@test isa( Atom(mass, position, velocity, symbol ), Atom)
@test isa( Atom(mass, position, velocity), Atom)

# Test Domain
bounds = SVector(SA_F64[0.0, 1.0], SA_F64[0.0, 1.0], SA_F64[0.0,1.0])
bounds_type = SVector( "p", "p", "p" )
@test isa(Domain(bounds, bounds_type), Domain)

# Test Configuration
r = Atom(mass, position, velocity, symbol)
@test isa( Configuration{3, 1, 0, 0, 0, 0, Float64}( 
    SVector{3, Atom}(r, r, r), 
    SVector{0, Angle}(), 
    SVector{0, Bond}(),
    SVector{0, Improper}(),
    SVector{0, Dihedral}(), 
    SVector(:Ar),
    SVector(18.0),
    SVector(0.5), 
    SVector(38.0),
    Domain(bounds, bounds_type),
    "lj" ), Configuration)
