
# /* ----------------------------------------------------------------------
#    compute Yi from Ui without storing Zi, looping over zlist indices
# ------------------------------------------------------------------------- */

function compute_yi(beta::AbstractFloat, snap::SNAP, runtime_arrays::RuntimeArrays)

    for ielem1 = 0:(snap.n_elements-1) 
        for j = 0:snap.twojmax
            jju = snap.prebuilt_arrays.idxu_block[j+1];
            mb = 0
            while 2*mb < j
                for ma = 0:j 
                    runtime_arrays.ylist_r[ielem1*idxu_max+jju+1] = 0.0;
                    runtime_arrays.ylist_i[ielem1*idxu_max+jju+1] = 0.0;
                    jju+=1;
                end
                mb += 1;
            end # end loop over ma, mb
        end # loop over j
    end

    for elem1 = 0:(snap.n_elements-1)
        for elem2 = 0:(snap.n_elements-1)
            for jjz = 1:snap.prebuilt_arrays.idxz_max
                j1     = snap.prebuilt_arrays.idxz[jjz].j1;
                j2     = snap.prebuilt_arrays.idxz[jjz].j2;
                j      = snap.prebuilt_arrays.idxz[jjz].j;
                ma1min = snap.prebuilt_arrays.idxz[jjz].ma1min;
                ma2max = snap.prebuilt_arrays.idxz[jjz].ma2max;
                na     = snap.prebuilt_arrays.idxz[jjz].na;
                mb1min = snap.prebuilt_arrays.idxz[jjz].mb1min;
                mb2max = snap.prebuilt_arrays.idxz[jjz].mb2max;
                nb     = snap.prebuilt_arrays.idxz[jjz].nb;

                cgblock = snap.prebuilt_arrays.cglist[snap.prebuilt_arrays.idxcg_block[j1+1, j2+1, j+1]+1:end];

                ztmp_r = 0.0;
                ztmp_i = 0.0;

                jju1 = snap.prebuilt_arrays.idxu_block[j1+1] + (j1 + 1) * mb1min;
                jju2 = snap.prebuilt_arrays.idxu_block[j2+1] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                for ib = 1:nb # (int ib = 0; ib < nb; ib++) {

                    suma1_r = 0.0;
                    suma1_i = 0.0;

                    u1_r = runtime_arrays.ulisttot_r[elem1*idxu_max+jju1+1:end];
                    u1_i = runtime_arrays.ulisttot_i[elem1*idxu_max+jju1+1:end];
                    u2_r = runtime_arrays.ulisttot_r[elem2*idxu_max+jju2+1:end];
                    u2_i = runtime_arrays.ulisttot_i[elem2*idxu_max+jju2+1:end];

                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;

                    for ia = 1:na # (int ia = 0; ia < na; ia++) {
                        suma1_r += cgblock[icga+1] * (u1_r[ma1+1] * u2_r[ma2+1] - u1_i[ma1+1] * u2_i[ma2+1]);
                        suma1_i += cgblock[icga+1] * (u1_r[ma1+1] * u2_i[ma2+1] + u1_i[ma1+1] * u2_r[ma2+1]);
                        ma1+=1;
                        ma2-=1;
                        icga += j2;
                    end # loop over ia


                    ztmp_r += cgblock[icgb+1] * suma1_r;
                    ztmp_i += cgblock[icgb+1] * suma1_i;

                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                end # loop over ib

                # // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
                # // find right y_list[jju] and beta[jjb] entries
                # // multiply and divide by j+1 factors
                # // account for multiplicity of 1, 2, or 3

                if (snap.bnorm_flag) 
                    ztmp_i /= j+1;
                    ztmp_r /= j+1;
                end

                jju = snap.prebuilt_arrays.idxz[jjz].jju;
                for elem3 = 0:(snap.n_elements-1) #(int elem3 = 0; elem3 < nelements; elem3++) {
                    # // pick out right beta value
                    if (j >= j1) 
                        jjb = snap.prebuilt_arrays.idxb_block[j1+1][j2+1][j+1];
                        itriple = ((elem1 * snap.n_elements + elem2) * snap.n_nelements + elem3) * snap.prebuilt_arrays.idxb_max + jjb;
                        if (j1 == j) 
                            if (j2 == j) 
                                betaj = 3*beta[itriple+1];
                            else 
                                betaj = 2*beta[itriple+1];
                            end
                        else betaj = beta[itriple+1];
                        end
                    elseif (j >= j2) 
                        jjb = snap.prebuilt_arrays.idxb_block[j+1][j2+1][j1+1];
                        itriple = ((elem3 * snap.n_elements + elem2) * snap.n_elements + elem1) * snap.prebuilt_arrays.idxb_max + jjb;
                        if (j2 == j) 
                            betaj = 2*beta[itriple+1];
                        else 
                            betaj = beta[itriple+1];
                        end
                    else 
                        jjb = idxb_block[j2+1][j+1][j1+1];
                        itriple = ((elem2 * snap.n_elements + elem3) * snap.n_elements + elem1) * snap.prebuilt_arrays.idxb_max + jjb;
                        betaj = beta[itriple+1];
                    end

                    if (~bnorm_flag) && (j1 > j)
                        betaj *= (j1 + 1) / (j + 1.0);
                    end
                    runtime_arrays.ylist_r[elem3 * idxu_max + jju + 1] += betaj * ztmp_r;
                    runtime_arrays.ylist_i[elem3 * idxu_max + jju + 1] += betaj * ztmp_i;
                end
            end # loop over jjz
        end
    end
end