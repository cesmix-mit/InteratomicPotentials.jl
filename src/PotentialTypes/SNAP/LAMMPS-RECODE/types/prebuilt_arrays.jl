struct SNA_BINDICES
    j1 :: Int
    j2 :: Int
    j  :: Int
end
struct SNA_ZINDICES 
    j1 :: Int
    j2 :: Int 
    j  :: Int
    ma1min :: Int 
    ma2max :: Int 
    mb1min :: Int 
    mb2max :: Int 
    na     :: Int 
    nb     :: Int 
    jju    :: Int
end
struct PrebuiltArrays{T}
    # Array indices, sizes
    idxcg_block :: Array{Int}
    idxu_block  :: Array{Int}
    idxb_block  :: Array{Int}
    idxz_block  :: Array{Int}
    idxcg_max :: Int
    idxu_max  :: Int
    idxb_max  :: Int
    idxz_max  :: Int

    idxb      :: Vector{SNA_BINDICES}
    idxz      :: Vector{SNA_ZINDICES}
    ncoeff    :: Int 


    # Clebsch Gordan
    rootpqarray :: Array{T}
    cglist      :: Vector{T}

    # bzero 
    wself       :: T
    bzero       :: Vector{T}
end


# /* ---------------------------------------------------------------------- */
## Initialize

function initialize_prebuilt_arrays(twojmax::Int, n_elements::Int, bzero_flag::Bool, bnorm_flag::Bool)
    jdim = twojmax+1
    ncount = 0

    wself = 1.0
    bzero = zeros(jdim)
    if (bzero_flag) 
        www = wself*wself*wself;
        for j = 0:twojmax 
            if (bnorm_flag)
                bzero[j+1] = www;
            else
                bzero[j+1] = www*(j+1);
            end
        end
    end


    idxcg_block = zeros(Int, jdim, jdim, jdim)
    idxu_block = zeros(Int, jdim)
    idxb_block = zeros(Int, jdim, jdim, jdim)
    idxz_block = zeros(Int, jdim, jdim, jdim)

    idxcg_count = 0;
    for j1 = 0:twojmax
        for j2 = 0:j1
            j = j1-j2
            while j <= minimum([twojmax, j1+j2]) 
                idxcg_block[j1+1,j2+1,j+1] = idxcg_count;
                for m1 = 0:j1 
                    for m2 = 0:j2 
                        idxcg_count+=1;
                    end
                end
                j+=2;
            end
        end
    end
    idxcg_max = idxcg_count
    #   // index list for uarray
    #   // need to include both halves

    idxu_count = 0;

    for j = 0:twojmax 
        idxu_block[j+1] = idxu_count;
        for mb = 0:j 
            for ma = 0:j 
                idxu_count+=1;
            end
        end
    end
  
    idxu_max = idxu_count;

    # // index list for beta and B

    idxb_count = 0;
    for j1 = 0:twojmax 
        for j2 = 0:j1 
            j = j1 - j2
            while j <= minimum([twojmax, j1+j2])
                if (j >= j1) 
                    idxb_count+=1;
                end
                j+=2;
            end
        end
    end

    idxb_max = idxb_count;
    idxb = Vector{SNA_BINDICES}(undef, idxb_max)
    for k = 1:idxb_max
        idxb[k] = SNA_BINDICES(0, 0, 0)
    end
    idxb_count = 0;
    for j1 = 0:twojmax 
        for j2 = 0:j1 
            j = j1-j2
            while j <= minimum([twojmax, j1+j2])
                if (j >= j1) 
                    idxb[idxb_count+1] = SNA_BINDICES(j1, j2, j)
                    idxb_count+=1;
                end
                j +=2;
            end
        end
    end        

    # // reverse index list for beta and b

    idxb_count = 0;
    for j1 = 0:twojmax 
        for j2 = 0:j1 
            j = j1 - j2
            while j <= minimum([twojmax, j1+j2])
                if (j >= j1) 
                    idxb_block[j1+1,j2+1,j+1] = idxb_count;
                    idxb_count+=1;
                end
                j +=2
            end
        end
    end
        
    # // index list for zlist
    idxz_count = 0;

    for j1 = 0:twojmax 
        for j2 = 0:j1 
            j = j1 - j2 
            while j <= minimum([twojmax, j1+j2])
                mb = 0
                while 2*mb <= j 
                    for ma = 0:j
                        idxz_count+=1;
                    end
                    mb += 1;
                end
                j+=2;
            end
        end
    end

    idxz_max = idxz_count;
    idxz = Vector{SNA_ZINDICES}(undef, idxz_max)
    for k = 1:idxz_max
        idxz[k] = SNA_ZINDICES(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    end

    idxz_count = 0;
    for j1 = 0:twojmax
        for j2 = 0:j1
            j = j1-j2 
            while j <= minimum([twojmax, j1+j2])
                idxz_block[j1+1,j2+1,j+1] = idxz_count;

                # // find right beta[jjb] entry
                # // multiply and divide by j+1 factors
                # // account for multiplicity of 1, 2, or 3
                mb = 0
                while 2*mb <= j
                    for ma = 0:j 
                        ma1min = maximum([ 0, (2 * ma - j - j2 + j1) ÷ 2]);
                        ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) ÷ 2;
                        na = minimum( [j1, (2 * ma - j + j2 + j1) ÷ 2] ) - ma1min + 1;
                        mb1min = maximum([0, (2 * mb - j - j2 + j1) ÷ 2]);
                        mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) ÷ 2;
                        nb = minimum([ j1, (2 * mb - j + j2 + j1) ÷ 2] ) - mb1min + 1;
                        # // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

                        jju = idxu_block[j+1] + (j+1)*mb + ma;

                        idxz[idxz_count+1] = SNA_ZINDICES(j1, j2, j, 
                                        ma1min, ma2max, mb1min, mb2max, 
                                        na, nb, jju)
                        idxz_count+=1;
                    end
                    mb +=1;
                end
                j+=2;
            end
        end
    end

    # Initialize Root array 
    rootpqarray = zeros(twojmax+2, twojmax+2)        # rootpqarray
    rootpqarray = init_rootpqarray(twojmax, rootpqarray)

    # Initialize Clebsch Gordan
    cglist      = zeros(idxcg_max)                  # cglist 
    cglist      = init_clebsch_gordan(twojmax, cglist)

    
    # Calculate number of coefficients
    ncoeff = ncount * n_elements * n_elements * n_elements

    return PrebuiltArrays(
    idxcg_block,
    idxu_block,
    idxb_block,
    idxz_block,
    idxcg_max,
    idxu_max,
    idxb_max,
    idxz_max,
    idxb,
    idxz,
    ncoeff,
    rootpqarray,
    cglist, 
    wself,
    bzero)
end


function deltacg(j1::Int, j2::Int, j::Int)
    sfaccg = factorial( (j1+j2+j) ÷ 2 + 1  )
    return sqrt( 
            factorial( (j1+j2 -j) ÷ 2 ) * 
            factorial( (j1-j2+j) ÷ 2  )  * 
            factorial( (-j1+j2+j) ÷ 2 ) / sfaccg
    )
end

function init_clebsch_gordan(twojmax::Int, cglist::Vector{T}) where T<: AbstractFloat
    idxcg_count = 0;
    for j1 = 0:twojmax
        for j2 = 0:j1
            j = j1 - j2
            while j <= minimum([twojmax, j1+j2])
                for m1 = 0:j1 
                    aa2 = 2 * m1 - j1;
    
                    for m2 = 0:j2   
                        bb2 = 2 * m2 - j2;
                        m = (aa2 + bb2 + j) ÷ 2 ;
        
                        if (m < 0) || (m > j) 
                            cglist[idxcg_count+1] = 0.0;
                            idxcg_count+=1;
                            continue
                        else
                            # println("j $j")
                            sum = 0.0;
                            min_z = maximum([0, maximum([-(j-j2+aa2) ÷ 2, -(j-j1-bb2) ÷ 2])])
                            max_z = minimum([ (j1 + j2 - j) ÷ 2, minimum([ (j1-aa2) ÷ 2, (j2+bb2) ÷ 2 ])])
                            for z = min_z:max_z
                                ifac = ( (z % 2)==0 ) ? 1.0 : -1.0;
                                sum += ifac /
                                    (factorial(Int(z)) *
                                    factorial( (j1 + j2 - j) ÷ 2 - z) *
                                    factorial( (j1 - aa2) ÷ 2 - z) *
                                    factorial( (j2 + bb2) ÷ 2 - z) *
                                    factorial( (j - j2 + aa2) ÷ 2 + z) *
                                    factorial( (j - j1 - bb2) ÷ 2 + z)) ;
                            end
            
                            cc2 = 2 * m - j;
                            dcg = deltacg(j1, j2, j);
                            sfaccg = sqrt(factorial( (j1 + aa2) ÷ 2)   *
                                        factorial( (j1 - aa2) ÷ 2)  *
                                        factorial( (j2 + bb2) ÷ 2)  *
                                        factorial( (j2 - bb2) ÷ 2)  *
                                        factorial( (j  + cc2) ÷ 2)  *
                                        factorial( (j  - cc2) ÷ 2) *
                                        (j + 1));
                            
                            # println("j1 $j1, j2 $j2, m1 $m1, m2 $m2, j $j, m, $m")
                        
                            #     println("dcg $dcg, sfaccg $sfaccg")
                            # cg = float(CG(j1, j2, j, m1, m2, m))
                            # println("cg $cg")
                            cglist[idxcg_count+1] = sum * dcg * sfaccg;
                            # cglist[idxcg_count+1] = cg
                            idxcg_count+=1;
                        end
                    end
                end
                j +=2;
            end
        end
    end
    return cglist
end

function init_rootpqarray(twojmax::Int, rootpqarray::Matrix{T}) where T<:AbstractFloat
    for p = 1:twojmax
        for q = 1:twojmax
            rootpqarray[p+1, q+1] = sqrt(p / q)
        end
    end
    return rootpqarray
end

