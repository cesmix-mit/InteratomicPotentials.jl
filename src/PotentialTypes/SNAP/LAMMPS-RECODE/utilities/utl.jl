include("compute_ui.jl")
include("compute_z.jl")
include("compute_b.jl")
include("compute_yi.jl")
include("compute_dij.jl")

function get_num_coeffs(twojmax :: Int, n_elems :: Int)
    J = twojmax 
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer!")
    end 
    return num_coeffs * n_elems * n_elems * n_elems
end