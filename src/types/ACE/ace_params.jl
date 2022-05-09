struct ACEParams <: BasisParameters
    species 
    body_order        :: Int
    polynomial_degree :: Int
    wL      
    csp
    r0 
    rcutoff
end

length(ace_params::ACEParams) = length(get_rpi(ace_params))

function get_rpi(ace_params::ACEParams)
    ACE1.rpi_basis(;
    species = ace_params.species,
    N       = ace_params.body_order - 1,
    maxdeg  = ace_params.polynomial_degree,
    D       = ACE1.SparsePSHDegree(; wL = ace_params.wL, csp = ace_params.csp),
    r0      = ace_params.r0, 
    rin     = 0.65*ace_params.r0,
    rcut    = ace_params.rcutoff,
    pin     = 0,
    )
end