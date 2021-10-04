mutable struct Atom
    Mass            :: Real
    Position        :: Vector{<:Real}
    Velocity        :: Vector{<:Real}
    Type            :: Union{Symbol, Nothing}
end

function Atom(Mass::Real, Position:: Vector{<:Real}, Velocity::Vector{<:Real})
    return Atom(Mass, Position, Velocity, nothing)
end
