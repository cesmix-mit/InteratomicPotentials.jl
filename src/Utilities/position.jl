####################### Positions #############################################
struct Position 
    x :: Float64
    y :: Float64 
    z :: Float64
    type ::Symbol
end
function Position(x, y, z)
    return Position(x, y, z, :nothing)
end
function norm(r::Position)
    return sqrt(r.x^2 + r.y^2 + r.z^2)
end

+(r1::Position, r2::Position) = (r1.type == r2.type ? Position(r1.x+r2.x, r1.y + r2.y, r1.z + r2.z, r1.type) : Position(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z, :nothing))
-(r1::Position, r2::Position) = (r1.type == r2.type ? Position(r1.x-r2.x, r1.y - r2.y, r1.z - r2.z, r1.type) : Position(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z, :nothing))
*(a::AbstractFloat, r::Position) = Position(a*r.x, a*r.y, a*r.z, r.type)
*(r::Position, a::AbstractFloat) = Position(a*r.x, a*r.y, a*r.z, r.type)

+(a::AbstractFloat, r::Position) = Position(a + r.x, a + r.y, a + r.z, r.type)
+(a::Vector{Float64}, r::Position) = Position(a[1] + r.x, a[2] + r.y, a[3] + r.z, r.type)
+(r::Position, a::Vector{Float64}) = Position(a[1] + r.x, a[2] + r.y, a[3] + r.z, r.type)
-(a::Vector{Float64}, r::Position) = Position(a[1] - r.x, a[2] - r.y, a[3] - r.z, r.type)

vec(r::Position) = [r.x, r.y, r.z]
vec(r::Vector{Position}) = [[ri.x, ri.y, ri.z] for ri in r]

function get_interparticle_distance(ri::Position, rj::Position)
    return norm(ri - rj)
end
function get_interparticle_distance(r::Vector{Position})
    n = length(r)
    distances = zeros(0)    
    for i = 1:n
        ri = r[i]
        for j = (i+1):n 
            d = get_interparticle_distance(ri, r[j])
            push!(distances, d, d)
        end
    end
    return distances
end