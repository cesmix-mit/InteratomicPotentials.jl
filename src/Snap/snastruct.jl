

mutable struct SnaStruct
    twojmax::Int32
    ncoeff::Int32
    ncoeffall::Int32
    nperdim::Int32
    idxb_max::Int32
    idxu_max::Int32
    idxz_max::Int32
    idxcg_max::Int32
    ntypes::Int32
    nelements::Int32    
    ndoubles::Int32   
    ntriples::Int32   
    beta_max::Int32   
    bnormflag::Int32
    chemflag::Int32    
    quadraticflag::Int32
    switchflag::Int32
    bzeroflag::Int32
    wselfallflag::Int32
    rcutfacflag::Int32
    twojmaxflag::Int32 

    wself::Float64
    rmin0::Float64
    rfac0::Float64
    rcutfac::Float64
    rcutmax::Float64    

    map::Vector{Int32}  
    idx_max::Vector{Int32} 
    idxz::Vector{Int32}
    idxz_block::Vector{Int32}
    idxb::Vector{Int32}
    idxb_block::Vector{Int32}
    idxu_block::Vector{Int32}
    idxcg_block::Vector{Int32}

    rcutsq::Vector{Float64}    
    radelem::Vector{Float64}
    wjelem::Vector{Float64} 
    bzero::Vector{Float64}
    coeffelem::Vector{Float64}          
    rootpqarray::Vector{Float64} 
    cglist::Vector{Float64}

    SnaStruct() = new();
end

function idxbcount(twojmax::Int32)
    count = 0
    for j1 = 0:twojmax        
        for j2 = 0:j1
            for j = (j1-j2):2:min(twojmax, j1+j2) 
                if (j >= j1) 
                    count = count + 1
                end
            end
        end
    end
    return count
end

function idxucount(twojmax::Int32)
    count = 0
    for j = 0:twojmax        
        for mb = 0:j
            for ma = 0:j
                count = count + 1
            end
        end
    end
    return count
end

function idxzcount(twojmax::Int32)
    count = 0
    for j1 = 0:twojmax        
        for j2 = 0:j1
            for j = (j1-j2):2:min(twojmax, j1+j2) 
                for mb = 0:j
                    if (2*mb <= j)
                        for ma = 0:j
                            count = count + 1
                        end
                    end
                end
            end
        end
    end
    return count
end

function idxcgcount(twojmax::Int32)
    count = 0
    for j1 = 0:twojmax        
        for j2 = 0:j1
            for j = (j1-j2):2:min(twojmax, j1+j2) 
                for mb = 0:j1
                    for ma = 0:j2
                        count = count + 1
                    end
                end
            end
        end
    end
    return count
end

function initializesna(ntypes::Int32, twojmax::Int32, chemflag::Int32)

jdim = twojmax + 1
jdimpq = twojmax + 2
if (chemflag == 0)    
    nelements = 1;
else
    nelements = ntypes 
end

sna = SnaStruct()

sna.idxb_max = idxbcount(twojmax)
sna.idxu_max = idxucount(twojmax)
sna.idxz_max = idxzcount(twojmax)
sna.idxcg_max = idxcgcount(twojmax)

sna.ntypes = ntypes
sna.twojmax = twojmax
sna.nelements = nelements;    
sna.ndoubles = nelements*nelements;  
sna.ntriples = nelements*nelements*nelements;  
sna.ncoeff = sna.idxb_max*sna.ntriples
sna.ncoeffall = sna.ncoeff + 1    

sna.map = zeros(Int32, ntypes+1)
sna.idxcg_block = zeros(Int32, jdim*jdim*jdim)
sna.idxz_block = zeros(Int32, jdim*jdim*jdim)
sna.idxb_block = zeros(Int32, jdim*jdim*jdim)
sna.idxu_block = zeros(Int32, jdim)
sna.idx_max = zeros(Int32, 5)
sna.idxz = zeros(Int32, sna.idxz_max*10)
sna.idxb = zeros(Int32, sna.idxb_max*3)

sna.rcutsq = zeros(Float64, (ntypes)*(ntypes))
sna.radelem = zeros(Float64, ntypes+1)
sna.wjelem = zeros(Float64, ntypes+1)
sna.bzero = zeros(Float64, jdim)
sna.rootpqarray = zeros(Float64, jdimpq*jdimpq)
sna.cglist = zeros(Float64, sna.idxcg_max)

return sna

end


