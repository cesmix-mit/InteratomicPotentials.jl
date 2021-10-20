struct Atom{T1<:Real, T2<:AbstractFloat}
    Mass            :: T1
    Position        :: SVector{3, T2}
    Velocity        :: SVector{3, T2}
    Type            :: Union{Symbol, Nothing}
end

function Atom(Mass :: Real, Position:: SVector{3, <:AbstractFloat}, Velocity::SVector{3, <:AbstractFloat})
    return Atom(Mass, Position, Velocity, nothing)
end
