# /* ----------------------------------------------------------------------
#    compute Bi by summing conj(Ui)*Zi
# ------------------------------------------------------------------------- */

function compute_bi(snap::SNAPParams, runtime_arrays::RuntimeArrays) 
    #   // for j1 = 0,...,twojmax
    #   //   for j2 = 0,twojmax
    #   //     for j = |j1-j2|,Min(twojmax,j1+j2),2
    #   //        b(j1,j2,j) = 0
    #   //        for mb = 0,...,jmid
    #   //          for ma = 0,...,j
    #   //            b(j1,j2,j) +=
    #   //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)
    n_elements = length(snap.elements)
    for ielem = 1:n_elements
        itriple = 0;
        idouble = 0;
        for elem1 = 0:(n_elements-1)
            for elem2 = 0:(n_elements-1)
                zptr_r = runtime_arrays.zlist_r[(idouble*snap.prebuilt_arrays.idxz_max+1):end];
                zptr_i = runtime_arrays.zlist_i[(idouble*snap.prebuilt_arrays.idxz_max+1):end];
    
                for elem3 = 0:(n_elements-1)
                    for jjb = 0:(snap.prebuilt_arrays.idxb_max-1)
                        j1 = snap.prebuilt_arrays.idxb[jjb+1].j1;
                        j2 = snap.prebuilt_arrays.idxb[jjb+1].j2;
                        j  = snap.prebuilt_arrays.idxb[jjb+1].j;

                        jjz = snap.prebuilt_arrays.idxz_block[j1+1,j2+1,j+1];
                        jju = snap.prebuilt_arrays.idxu_block[j+1];
                        sumzu = 0.0
                        mb = 0

                        while 2*mb < j
                            for ma = 0:j 
                                sumzu += (runtime_arrays.ulisttot_r[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_r[jjz+1] +
                                        runtime_arrays.ulisttot_i[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_i[jjz+1])
                                jjz+=1
                                jju+=1
                            end
                            mb +=1;
                        end # loop over ma, mb
                        # // For j even, handle middle column
                        if (j % 2 == 0) 
                            mb = j ÷ 2;
                            ma = 0
                            # Strict inequality
                            while ma < mb
                                sumzu += runtime_arrays.ulisttot_r[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_r[jjz+1] +
                                        runtime_arrays.ulisttot_i[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_i[jjz+1];
                                jjz+=1;
                                jju+=1;
                                ma+=1
                            end
                            sumzu += 0.5 * (runtime_arrays.ulisttot_r[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_r[jjz+1] +
                                        runtime_arrays.ulisttot_i[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_i[jjz+1]);

                        end # if jeven

                        runtime_arrays.blist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1] = 2.0 * sumzu;
                    end
                    itriple+=1;
                end
                idouble+=1;
            end
        end
    
        # // apply bzero shift
    
        if (snap.bzero_flag) 
            if (~snap.wselfall_flag) 
                itriple = (ielem-1)*( n_elements+1)*n_elements+(ielem-1);
                for jjb = 0:(snap.prebuilt_arrays.idxb_max-1)
                    j = snap.prebuilt_arrays.idxb[jjb+1].j;
                    runtime_arrays.blist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1] -= snap.prebuilt_arrays.bzero[j+1];
                end # loop over JJ
            else 
                itriple = 0;
                for elem1=0:(n_elements-1)
                    for elem2 = 0:(n_elements-1)
                        for elem3 = 0:(n_elements-1)
                            for jjb = 0:(snap.prebuilt_arrays.idxb_max-1)
                                j = snap.prebuilt_arrays.idxb[jjb+1].j;
                                runtime_arrays.blist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1] -= snap.prebuilt_arrays.bzero[j+1];
                            end # loop over JJ
                            itriple+=1;
                        end # loop over elem3
                    end # loop over elem1,elem2
                end
            end
        end
    end
end
    
    
# /* ----------------------------------------------------------------------
#    calculate derivative of Bi w.r.t. atom j
#    variant using indexlist for j1,j2,j
#    variant using symmetry relation
# ------------------------------------------------------------------------- */

function compute_dbidrj(jelem::Int, snap::SNAPParams, runtime_arrays::RuntimeArrays)

    # // for j1 = 0,...,twojmax
    # //   for j2 = 0,twojmax
    # //     for j = |j1-j2|,Min(twojmax,j1+j2),2
    # //        zdb = 0
    # //        for mb = 0,...,jmid
    # //          for ma = 0,...,j
    # //            zdb +=
    # //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
    # //        dbdr(j1,j2,j) += 2*zdb
    # //        zdb = 0
    # //        for mb1 = 0,...,j1mid
    # //          for ma1 = 0,...,j1
    # //            zdb +=
    # //              Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
    # //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j1+1)
    # //        zdb = 0
    # //        for mb2 = 0,...,j2mid
    # //          for ma2 = 0,...,j2
    # //            zdb +=
    # //              Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)
    # //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j2+1)

    # elem3 = elem_duarray;
    n_elements = length(snap.elements)
    elem3 = jelem-1;
    for jjb = 0:(snap.prebuilt_arrays.idxb_max-1) 
        j1 = snap.prebuilt_arrays.idxb[jjb+1].j1;
        j2 = snap.prebuilt_arrays.idxb[jjb+1].j2;
        j  = snap.prebuilt_arrays.idxb[jjb+1].j;

        # // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

        for elem1 = 0:(n_elements-1) 
            for elem2 = 0:(n_elements-1)
                # println("j1 $j1, j2 $j2, j $j")

                jjz = snap.prebuilt_arrays.idxz_block[j1+1,j2+1,j+1];
                jju = snap.prebuilt_arrays.idxu_block[j+1];
                idouble = elem1*n_elements+elem2;
                itriple = (elem1*n_elements+elem2)*n_elements+elem3;
                dbdr = runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :];
                zptr_r = runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max+1:end];
                zptr_i = runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max+1:end];

                sumzdu_r = zeros(3)
                
                mb = 0
                while 2*mb < j 
                    for ma = 0:j
                        dudr_r = runtime_arrays.dulist_r[jju, :];
                        dudr_i = runtime_arrays.dulist_i[jju, :];
                        for k = 1:3 
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end # loop over ma mb
                    mb+=1;
                end

                # // For j even, handle middle column

                if (j % 2 == 0) 
                    mb = floor(Int, j / 2);
                    for ma = 0:(mb-1) 
                        dudr_r = runtime_arrays.dulist_r[jju+1, :];
                        dudr_i = runtime_arrays.dulist_i[jju+1, :];
                        for k = 1:3
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
            
                    dudr_r = runtime_arrays.dulist_r[jju+1, :];
                    dudr_i = runtime_arrays.dulist_i[jju+1, :];
                    for k = 1:3 
                        sumzdu_r[k] +=
                            (dudr_r[k] * zptr_r[jjz+1] +
                            dudr_i[k] * zptr_i[jjz+1]) * 0.5;
                    end
                    # // jjz++;
                    # // jju++;
                end # if jeven

                for k = 1:3
                    dbdr[k] += 2.0 * sumzdu_r[k];
                end
                runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :] = dbdr

                # // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
                j1fac   = (j + 1) / (j1 + 1.0);
                idouble = elem1*n_elements+elem2;
                itriple = (elem3*n_elements+elem2)*n_elements+elem1;
                dbdr    = runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :];
                jjz     = snap.prebuilt_arrays.idxz_block[j+1,j2+1,j1+1];
                jju     = snap.prebuilt_arrays.idxu_block[j1+1];
                zptr_r  = runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max+1:end];
                zptr_i  = runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max+1:end];

                for k = 1:3 
                    sumzdu_r[k] = 0.0;
                end

                mb = 0
                while 2*mb < j1
                    for ma = 0:j1 
                        dudr_r = runtime_arrays.dulist_r[jju+1, :];
                        dudr_i = runtime_arrays.dulist_i[jju+1, :];
                        for k = 1:3 
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
                    mb+=1
                end # loop over ma mb

                # // For j1 even, handle middle column

                if (j1 % 2 == 0) 
                    mb = floor(Int, j1 / 2);
                    for ma = 0:(mb-1)
                        dudr_r = runtime_arrays.dulist_r[jju+1, :];
                        dudr_i = runtime_arrays.dulist_i[jju+1, :];
                        for k = 1:3
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
            
                    dudr_r = runtime_arrays.dulist_r[jju+1, :];
                    dudr_i = runtime_arrays.dulist_i[jju+1, :];
                    for k = 1:3
                        sumzdu_r[k] +=
                            (dudr_r[k] * zptr_r[jjz+1] +
                            dudr_i[k] * zptr_i[jjz+1]) * 0.5;
                    end
                    # // jjz++;
                    # // jju++;
                end # if j1even

                for k = 1:3
                    if (snap.bnorm_flag)
                        dbdr[k] += 2.0 * sumzdu_r[k];
                    else
                        dbdr[k] += 2.0 * sumzdu_r[k] * j1fac;
                    end
                end
                runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :] = dbdr
                # // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)

                j2fac = (j + 1) / (j2 + 1.0);

                idouble = elem2*n_elements+elem1;
                itriple = (elem1*n_elements+elem3)*n_elements+elem2;
                dbdr = runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :];
                jjz = snap.prebuilt_arrays.idxz_block[j+1,j1+1,j2+1];
                jju = snap.prebuilt_arrays.idxu_block[j2+1];
                zptr_r = runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max+1:end];
                zptr_i = runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max+1:end];

                for k = 1:3
                    sumzdu_r[k] = 0.0;
                end

                mb = 0
                while 2*mb < j2 
                    for ma = 0:j2 
                        dudr_r = runtime_arrays.dulist_r[jju+1, :];
                        dudr_i = runtime_arrays.dulist_i[jju+1, :];
                        for k = 1:3
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
                    mb +=1;
                end # loop over ma mb

                # // For j2 even, handle middle column

                if (j2 % 2 == 0) 
                    mb = floor(Int, j2 / 2);
                    for ma = 0:(mb-1)
                        dudr_r = runtime_arrays.dulist_r[jju+1, :];
                        dudr_i = runtime_arrays.dulist_i[jju+1, :];
                        for k = 1:3 
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
                    dudr_r = runtime_arrays.dulist_r[jju+1, :];
                    dudr_i = runtime_arrays.dulist_i[jju+1, :];
                    for k = 1:3
                        sumzdu_r[k] +=
                            (dudr_r[k] * zptr_r[jjz+1] +
                            dudr_i[k] * zptr_i[jjz+1]) * 0.5;
                    end
                    # // jjz++;
                    # // jju++;
                end # if j2even

                for k = 1:3 
                    if (snap.bnorm_flag)
                        dbdr[k] += 2.0 * sumzdu_r[k];
                    else
                        dbdr[k] += 2.0 * sumzdu_r[k] * j2fac;
                    end
                end
                runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1, :] = dbdr;
            end
        end
    end

end
