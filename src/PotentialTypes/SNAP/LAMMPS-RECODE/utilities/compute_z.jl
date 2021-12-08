function compute_zi(snap::SNAPParams, runtime_arrays::RuntimeArrays)
    idouble = 0;
    n_elements = length(snap.elements)
    for elem1 = 0:(n_elements-1)
        for elem2 = 0:(n_elements-1)

            for jjz = 1:snap.prebuilt_arrays.idxz_max
                j1 = snap.prebuilt_arrays.idxz[jjz].j1;
                j2 = snap.prebuilt_arrays.idxz[jjz].j2;
                j = snap.prebuilt_arrays.idxz[jjz].j;
                ma1min = snap.prebuilt_arrays.idxz[jjz].ma1min;
                ma2max = snap.prebuilt_arrays.idxz[jjz].ma2max;
                na = snap.prebuilt_arrays.idxz[jjz].na;
                mb1min = snap.prebuilt_arrays.idxz[jjz].mb1min;
                mb2max = snap.prebuilt_arrays.idxz[jjz].mb2max;
                nb = snap.prebuilt_arrays.idxz[jjz].nb;

                cgblock = snap.prebuilt_arrays.cglist[snap.prebuilt_arrays.idxcg_block[j1+1, j2+1, j+1]+1:end];
                runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max + jjz] = 0.0;
                runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max + jjz] = 0.0;

                jju1 = snap.prebuilt_arrays.idxu_block[j1+1] + (j1 + 1) * mb1min;
                jju2 = snap.prebuilt_arrays.idxu_block[j2+1] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                for ib = 0:(nb-1)

                    suma1_r = 0.0;
                    suma1_i = 0.0;
                    u1_r = runtime_arrays.ulisttot_r[(elem1*snap.prebuilt_arrays.idxu_max+jju1+1):end];
                    u1_i = runtime_arrays.ulisttot_i[(elem1*snap.prebuilt_arrays.idxu_max+jju1+1):end];
                    u2_r = runtime_arrays.ulisttot_r[(elem2*snap.prebuilt_arrays.idxu_max+jju2+1):end];
                    u2_i = runtime_arrays.ulisttot_i[(elem2*snap.prebuilt_arrays.idxu_max+jju2+1):end];

                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;

                    for ia = 0:(na-1)
                        suma1_r += cgblock[icga+1] * (u1_r[ma1+1] * u2_r[ma2+1] - u1_i[ma1+1] * u2_i[ma2+1]);
                        suma1_i += cgblock[icga+1] * (u1_r[ma1+1] * u2_i[ma2+1] + u1_i[ma1+1] * u2_r[ma2+1]);
                        ma1+=1;
                        ma2-=1;
                        icga += j2;
                    end # loop over ia
                    runtime_arrays.zlist_r[jjz] += cgblock[icgb+1] * suma1_r;
                    runtime_arrays.zlist_i[jjz] += cgblock[icgb+1] * suma1_i;

                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                end # loop over ib
                if (snap.bnorm_flag) 
                    runtime_arrays.zlist_r[jjz] /= (j+1);
                    runtime_arrays.zlist_i[jjz] /= (j+1);
                end
            end # loop over jjz
            idouble+=1;
        end
    end
end
