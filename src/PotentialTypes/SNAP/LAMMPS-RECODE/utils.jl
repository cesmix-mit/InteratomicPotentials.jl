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
# /* ---------------------------------------------------------------------- */

# void SNA::compute_ncoeff()
# {
#   int ncount;

#   ncount = 0;

#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2;
#            j <= MIN(twojmax, j1 + j2); j += 2)
#         if (j >= j1) ncount++;

#   ndoubles = nelements*nelements;
#   ntriples = nelements*nelements*nelements;
#   if (chem_flag)
#     ncoeff = ncount*ntriples;
#   else
#     ncoeff = ncount;
# }

function compute_ncoeff_build_list(twojmax, n_elements, chem_flag)
    jdim = twojmax+1
    ncount = 0
    idxcg_count = 0
    idxu_count = 0
    idxb_count = 0
    idxz_count = 0

    idxcg_block = zeros(Int, jdim, jdim, jdim)
    idxu_block = zeros(Int, jdim)
    idxb_block = zeros(Int, jdim, jdim, jdim)
    idxz_block = zeros(Int, jdim, jdim, jdim)
    for j1 = 0:twojmax
        for j2 = 0:j1
            for j = (j1-j2):2:minimum([twojmax, j1+j2])
                if j >= j1 
                    ncount += 1

                    idxb_block[j1+1, j2+1, j+1] = idxb_count
                    idxb_count += 1
                    
                end

                #idxcg
                idxcg_block[j1+1, j2+1, j+1] = idxcg_count
                for m1 = 0:j1
                    for m2 = 0:j2 
                        idxcg_count += 1
                    end
                end

                #idxz
                for mb = 0:floor(0.5*j)
                    for ma = 0:j;
                      idxz_count +=1
                    end
                end
                
            end
        end

        #idxu_block 
        idxu_block[j1+1] = idxu_count
        for mb = 0:j1 
            for ma = 0:j1
                idxu_count += 1
            end
        end

    end

    idxcg_max = idxcg_count
    idxu_max = idxu_count
    idxb_max = idxb_count
    idxz_max = idxz_count 

    idxb = Vector{SNA_BINDICES}(undef, idxb_max)
    idxz = Vector{SNA_ZINDICES}(undef, idxz_max)
    ## Second loop for idxb
    idxb_count = 0
    idxz_count = 0
    for j1 = 0:twojmax
        for j2 = 0:j1
            for j = (j1-j2):2:minimum([twojmax, j1+j2])
                if j >= j1 
                    idxb_count += 1
                    idxb[idxb_count] = SNA_BINDICES(j1, j2, j)
                end

                idxz_block[j1+1, j2+1, j+1] = idxz_count;

                # // find right beta[jjb] entry
                # // multiply and divide by j+1 factors
                # // account for multiplicity of 1, 2, or 3

                for mb = 0:floor(0.5*j)
                    for ma = 0:j 
                        idxz_count += 1
                        ma1min = maximum([0, (2 * ma - j - j2 + j1)/2])
                        ma2max = (2 * ma - j - (2*ma1min - j1) + j2) / 2
                        na     = minimum([j1, (2 * ma - j + j2 + j1) / 2]) - ma1min + 1 
                        mb1min = maximum([0, (2 * mb - j - j2 + j1) / 2])
                        mb2max = (2 * mb - j - (2*mb1min - j1) + j2)/2
                        nb     =  minimum([j1, (2*mb - j + j2 + j1) / 2]) - mb1min + 1

                        jju = idxu_block[j+1] + (j+1)*mb + ma 
                        idxz[idxz_count] = SNA_ZINDICES(j1, j2, j, 
                                            ma1min, ma2max, na, 
                                            mb1min, mb2max, nb, jju)
                    end
                end
                
            end
        end
    end

    # Calculate number of coefficients
    ndoubles = n_elements * n_elements
    ntriples = n_elements * n_elements * n_elements
    if chem_flag
        ncoeff = ncount * ntriples
    else
        ncoeff = ncount 
    end
    return ncoeff, idxb, idxz, idxcg
end

# /* ---------------------------------------------------------------------- */


# void SNA::build_indexlist()
# {

#   // index list for cglist

#   int jdim = twojmax + 1;
#   memory->create(idxcg_block, jdim, jdim, jdim,
#                  "sna:idxcg_block");

#   int idxcg_count = 0;
#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
#         idxcg_block[j1][j2][j] = idxcg_count;
#         for (int m1 = 0; m1 <= j1; m1++)
#           for (int m2 = 0; m2 <= j2; m2++)
#             idxcg_count++;
#       }
#   idxcg_max = idxcg_count;

#   // index list for uarray
#   // need to include both halves

#   memory->create(idxu_block, jdim,
#                  "sna:idxu_block");

#   int idxu_count = 0;

#   for (int j = 0; j <= twojmax; j++) {
#     idxu_block[j] = idxu_count;
#     for (int mb = 0; mb <= j; mb++)
#       for (int ma = 0; ma <= j; ma++)
#         idxu_count++;
#   }
#   idxu_max = idxu_count;

#   // index list for beta and B

#   int idxb_count = 0;
#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
#         if (j >= j1) idxb_count++;

#   idxb_max = idxb_count;
#   idxb = new SNA_BINDICES[idxb_max];

#   idxb_count = 0;
#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
#         if (j >= j1) {
#           idxb[idxb_count].j1 = j1;
#           idxb[idxb_count].j2 = j2;
#           idxb[idxb_count].j = j;
#           idxb_count++;
#         }

#   // reverse index list for beta and b

#   memory->create(idxb_block, jdim, jdim, jdim,
#                  "sna:idxb_block");
#   idxb_count = 0;
#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
#         if (j >= j1) {
#           idxb_block[j1][j2][j] = idxb_count;
#           idxb_count++;
#         }
#       }

#   // index list for zlist

#   int idxz_count = 0;

#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
#         for (int mb = 0; 2*mb <= j; mb++)
#           for (int ma = 0; ma <= j; ma++)
#             idxz_count++;

#   idxz_max = idxz_count;
#   idxz = new SNA_ZINDICES[idxz_max];

#   memory->create(idxz_block, jdim, jdim, jdim,
#                  "sna:idxz_block");

#   idxz_count = 0;
#   for (int j1 = 0; j1 <= twojmax; j1++)
#     for (int j2 = 0; j2 <= j1; j2++)
#       for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
#         idxz_block[j1][j2][j] = idxz_count;

#         // find right beta[jjb] entry
#         // multiply and divide by j+1 factors
#         // account for multiplicity of 1, 2, or 3

#         for (int mb = 0; 2*mb <= j; mb++)
#           for (int ma = 0; ma <= j; ma++) {
#             idxz[idxz_count].j1 = j1;
#             idxz[idxz_count].j2 = j2;
#             idxz[idxz_count].j = j;
#             idxz[idxz_count].ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
#             idxz[idxz_count].ma2max = (2 * ma - j - (2 * idxz[idxz_count].ma1min - j1) + j2) / 2;
#             idxz[idxz_count].na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count].ma1min + 1;
#             idxz[idxz_count].mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
#             idxz[idxz_count].mb2max = (2 * mb - j - (2 * idxz[idxz_count].mb1min - j1) + j2) / 2;
#             idxz[idxz_count].nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count].mb1min + 1;
#             // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

#             const int jju = idxu_block[j] + (j+1)*mb + ma;
#             idxz[idxz_count].jju = jju;

#             idxz_count++;
#           }
#       }
# }
