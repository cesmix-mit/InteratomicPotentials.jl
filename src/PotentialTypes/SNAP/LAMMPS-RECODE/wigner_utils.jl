
# ----------------------------------------------------------------------
# compute Wigner U-functions for one neighbor
# ------------------------------------------------------------------------- */
include("utils.jl") 

function compute_ui(jnum::Int, ielem::Int, 
                    snap::SNAPParams, 
                    runtime_arrays :: RuntimeArrays)
    for j = 1:jnum
        x = snap.prebuilt_arrays.rij[j, 1]
        y = snap.prebuilt_arrays.rij[j, 2]
        z = snap.prebuilt_arrays.rij[j, 3]
        rsq = x*x + y*y + z*z 
        r   = sqrt(rsq) 

        theta0 = (r - snap.rmin0) * snap.rfac0 * pi / (snap.prebuilt_arrays.rcutij[j] - snap.rmin0)
        z0 = r / tan(theta0)
        println("$x, $y, $z, $rsq, $r, $theta0, $z0")
        compute_uarray(x, y, z, z0, r, j, 
                        snap, 
                        runtime_arrays)
        if snap.chem_flag
            add_uarraytot(r, snap.prebuilt_arrays.rcutij[j], snap.prebuilt_arrays.wj[j],
             j, snap.prebuilt_arrays.element[j], snap, runtime_arrays)
        else
            add_uarraytot(r, snap.prebuilt_arrays.rcutij[j], snap.prebuilt_arrays.wj[j],
             j, 1, snap, runtime_arrays)
        end
    end
end


function compute_zi(snap::SNAPParams, runtime_arrays::RuntimeArrays)
    idouble = 0;
    for elem1 = 0:(snap.n_elements-1)
        for elem2 = 0:(snap.n_elements-1)

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

                cgblock = snap.prebuilt_arrays.cglist .+ snap.prebuilt_arrays.idxcg_block[j1+1, j2+1, j+1];

                runtime_arrays.zlist_r[jjz] = 0.0;
                runtime_arrays.zlist_i[jjz] = 0.0;

                jju1 = snap.prebuilt_arrays.idxu_block[j1+1] + (j1 + 1) * mb1min;
                jju2 = snap.prebuilt_arrays.idxu_block[j2+1] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                for ib = 0:(nb-1)

                    suma1_r = 0.0;
                    suma1_i = 0.0;
                    ##### STOPPED HERE #####
                    u1_r = runtime_arrays.ulisttot_r[(elem1*snap.prebuilt_arrays.idxu_max+jju1+1):end];
                    u1_i = runtime_arrays.ulisttot_i[(elem1*snap.prebuilt_arrays.idxu_max+jju1+1):end];
                    println((elem2*snap.prebuilt_arrays.idxu_max+jju2+1))
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

                    runtime_arrays.zlist_r[jjz+1] += cgblock[icgb+1] * suma1_r;
                    runtime_arrays.zlist_r[jjz+1] += cgblock[icgb+1] * suma1_i;

                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                end # loop over ib
                if (snap.bnorm_flag) 
                    runtime_arrays.zlist_r[jjz+1] /= (j+1);
                    runtime_arrays.zlist_r[jjz+1] /= (j+1);
                end
            end # loop over jjz
            idouble+=1;
        end
    end
end

# /* ----------------------------------------------------------------------
#    compute Yi from Ui without storing Zi, looping over zlist indices
# ------------------------------------------------------------------------- */

function compute_yi(beta::AbstractFloat, snap::SNAPParams, runtime_arrays::RuntimeArrays)

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

                cgblock = snap.prebuilt_arrays.cglist .+ snap.prebuilt_arrays.idxcg_block[j1+1][j2+1][j+1];

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

# /* ----------------------------------------------------------------------
#    compute dEidRj
# ------------------------------------------------------------------------- */

function compute_deidrj(dedr::Vector{AbstractFloat}, jelem::Int)

    for k = 1:3 
        dedr[k] = 0.0;
    end

    for j = 0:snap.twojmax #(int j = 0; j <= twojmax; j++) {
        jju = idxu_block[j+1];

        mb = 0
        while 2*mb < j
            for ma = 0:j 

                dudr_r = dulist_r[jju+1];
                dudr_i = dulist_i[jju+1];
                jjjmambyarray_r = ylist_r[(jelem-1)*idxu_max+jju+1];
                jjjmambyarray_i = ylist_i[(jelem-1)*idxu_max+jju+1];

                for k = 1:3
                    dedr[k] +=
                    dudr_r[k] * jjjmambyarray_r +
                    dudr_i[k] * jjjmambyarray_i;
                end
                jju+=1;
            end
            mb += 1
        end # loop over ma mb

        # // For j even, handle middle column

        if (j%2 == 0) 
            mb = floor(Int, j/2);
            for ma = 0:(mb-1)
                dudr_r = dulist_r[jju+1, :];
                ddudr_i = dulist_i[jju+1, :];
                jjjmambyarray_r = ylist_r[(jelem-1)*idxu_max+jju+1];
                jjjmambyarray_i = ylist_i[(jelem-1)*idxu_max+jju+1];

                for k=1:3
                    dedr[k] +=
                        dudr_r[k] * jjjmambyarray_r +
                        dudr_i[k] * jjjmambyarray_i;
                    
                end
                jju+=1;
            end

            dudr_r = dulist_r[jju+1];
            dudr_i = dulist_i[jju+1];
            jjjmambyarray_r = ylist_r[(jelem-1)*idxu_max+jju+1];
            jjjmambyarray_i = ylist_i[(jelem-1)*idxu_max+jju+1];

            for k = 1:3
                dedr[k] +=
                    (dudr_r[k] * jjjmambyarray_r +
                    dudr_i[k] * jjjmambyarray_i)*0.5;
                # // jju++;
            end
        end # if j even

    end # loop over j

    for k = 1:3
        dedr[k] *= 2.0;
    end
end

# /* ----------------------------------------------------------------------
#    compute Bi by summing conj(Ui)*Zi
# ------------------------------------------------------------------------- */

function compute_bi(ielem::Int, snap::SNAPParams, runtime_arrays::RuntimeArrays) 
#   // for j1 = 0,...,twojmax
#   //   for j2 = 0,twojmax
#   //     for j = |j1-j2|,Min(twojmax,j1+j2),2
#   //        b(j1,j2,j) = 0
#   //        for mb = 0,...,jmid
#   //          for ma = 0,...,j
#   //            b(j1,j2,j) +=
#   //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

    itriple = 0;
    idouble = 0;
    for elem1 = 0:(snap.n_elements-1)
        for elem2 = 0:(snap.n_elements-1)

            zptr_r = runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max+1:end];
            zptr_i = runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max+1:end];

            for elem3 = 0:(snap.n_elements-1)
                for jjb = 0:(snap.prebuilt_arrays.idxb_max-1)
                    j1 = snap.prebuilt_arrays.idxb[jjb+1].j1;
                    j2 = snap.prebuilt_arrays.idxb[jjb+1].j2;
                    j  = snap.prebuilt_arrays.idxb[jjb+1].j;

                    jjz = snap.prebuilt_arrays.idxz_block[j1+1][j2+1][j+1];
                    jju = snap.prebuilt_arrays.idxu_block[j+1];
                    sumzu = 0.0;
                    mb = 0
                    while 2*mb < j
                        for ma = 0:j 
                            sumzu += runtime_arrays.ulisttot_r[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_r[jjz+1] +
                                    runtime_arrays.ulisttot_i[elem3*snap.prebuilt_arrays.idxu_max+jju+1] * zptr_i[jjz+1];
                            jjz+=1;
                            jju+=1;
                        end
                        mb +=1;
                    end # loop over ma, mb

                    # // For j even, handle middle column

                    if (j % 2 == 0) 
                        mb = floor(Int, j / 2);
                        for ma = 0:(mb-1) 
                            sumzu += runtime_arrays.ulisttot_r[elem3*idxu_max+jju+1] * zptr_r[jjz+1] +
                                    runtime_arrays.ulisttot_i[elem3*idxu_max+jju+1] * zptr_i[jjz+1];
                            jjz+=1;
                            jju+=1;
                        end
                        sumzu += 0.5 * (runtime_arrays.ulisttot_r[elem3*idxu_max+jju+1] * zptr_r[jjz+1] +
                                    runtime_arrays.ulisttot_i[elem3*idxu_max+jju+1] * zptr_i[jjz+1]);
                    end # if jeven

                    runtime_arrays.blist[itriple*idxb_max+jjb+1] = 2.0 * sumzu;
                end
                itriple+=1;
            end
            idouble+=1;
        end
    end

    # // apply bzero shift

    if (snap.bzero_flag) 
        if (~snap.wselfall_flag) 
            itriple = (ielem*nelements+ielem)*nelements+ielem;
            for jjb = 0:(snap.idxb_max-1)
                j = idxb[jjb+1].j;
                runtime_arrays.blist[itriple*idxb_max+jjb+1] -= bzero[j+1];
            end # loop over JJ
        else 
            itriple = 0;
            for elem1=0:(snap.n_elements-1)
                for elem2 = 0:(snap.n_elements-1)
                    for elem3 = 0:(snap.n_elements-1)
                        for jjb = 0:(snap.prebuilt_arrays.idxb_max-1)
                            j = idxb[jjb+1].j;
                            runtime_arrays.blist[itriple*idxb_max+jjb+1] -= bzero[j+1];
                        end # loop over JJ
                        itriple+=1;
                    end # loop over elem3
                end # loop over elem1,elem2
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
    elem3 = jelem-1;
    for jjb = 0:(snap.prebuilt_arrays.idxb_max-1) 
        j1 = snap.prebuilt_arrays.idxb[jjb+1].j1;
        j2 = snap.prebuilt_arrays.idxb[jjb+1].j2;
        j  = snap.prebuilt_arrays.idxb[jjb+1].j;

        # // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

        for elem1 = 0:(snap.n_elements-1) 
            for elem2 = 0:(snap.n_elements-1)

                jjz = snap.prebuilt_arrays.idxz_block[j1+1][j2+1][j+1];
                jju = snap.prebuilt_arrays.idxu_block[j+1];
                idouble = elem1*nelements+elem2;
                itriple = (elem1*nelements+elem2)*nelements+elem3;
                dbdr = runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1];
                zptr_r = runtime_arrays.zlist_r[idouble*snap.prebuilt_arrays.idxz_max+1:end];
                zptr_i = runtime_arrays.zlist_i[idouble*snap.prebuilt_arrays.idxz_max+1:end];

                for k = 1:3
                    sumzdu_r[k] = 0.0;
                end
                
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
            
                    dudr_r = runtime_arrays.dulist_r[jju+1];
                    dudr_i = runtime_arrays.dulist_i[jju+1];
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
                runtime_arrays.dblist[itriple*snap.prebuilt_arrays.idxb_max+jjb+1] = dbdr

                # // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
                j1fac   = (j + 1) / (j1 + 1.0);
                idouble = elem1*nelements+elem2;
                itriple = (elem3*nelements+elem2)*nelements+elem1;
                dbdr    = runtime_arrays.dblist[itriple*idxb_max+jjb+1, :];
                jjz     = snap.prebuilt_arrays.idxz_block[j+1][j2+1][j1+1];
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
                        dudr_r = runtime_arrays.dulist_r[jju+1];
                        dudr_i = runtime_arrays.dulist_i[jju+1];
                        for k = 1:3
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
            
                    dudr_r = runtime_arrays.dulist_r[jju+1];
                    dudr_i = runtime_arrays.dulist_i[jju+1];
                    for k = 1:3
                        sumzdu_r[k] +=
                            (dudr_r[k] * zptr_r[jjz+1] +
                            dudr_i[k] * zptr_i[jjz+1]) * 0.5;
                    end
                    # // jjz++;
                    # // jju++;
                end # if j1even

                for k = 1:3
                    if (bnorm_flag)
                        dbdr[k] += 2.0 * sumzdu_r[k];
                    else
                        dbdr[k] += 2.0 * sumzdu_r[k] * j1fac;
                    end
                end
                runtime_arrays.dblist[itriple*idxb_max+jjb+1, :] = dbdr
                # // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)

                j2fac = (j + 1) / (j2 + 1.0);

                idouble = elem2*nelements+elem1;
                itriple = (elem1*nelements+elem3)*nelements+elem2;
                dbdr = runtime_arrays.dblist[itriple*idxb_max+jjb+1, :];
                jjz = snap.prebuilt_arrays.idxz_block[j+1][j1+1][j2+1];
                jju = snap.prebuilt_arrays.idxu_block[j2];
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
                end # loop over ma mb

                # // For j2 even, handle middle column

                if (j2 % 2 == 0) 
                    mb = floor(Int, j2 / 2);
                    for ma = 0:(mb-1)
                        dudr_r = dulist_r[jju+1, :];
                        dudr_i = dulist_i[jju+1, :];
                        for k = 1:3 
                            sumzdu_r[k] +=
                                dudr_r[k] * zptr_r[jjz+1] +
                                dudr_i[k] * zptr_i[jjz+1];
                        end
                        jjz+=1;
                        jju+=1;
                    end
                    dudr_r = runtime_arrays.dulist_r[jju+1];
                    dudr_i = runtime_arrays.dulist_i[jju+1];
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
                runtime_arrays.dblist[itriple*idxb_max+jjb+1, :] = dbdr;
            end
        end
    end

end


####################################################################################


function compute_uarray(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                      z0::AbstractFloat, r::AbstractFloat, jj::Int, 
                      snap::SNAPParams, runtime_arrays::RuntimeArrays)

    # compute Cayley-Klein parameters for unit quaternion

    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;
    # // VMK Section 4.8.2


    ulist_r = runtime_arrays.ulist_r_ij[jj, :];
    ulist_i = runtime_arrays.ulist_i_ij[jj, :];
    ulist_r[1] = 1.0;
    ulist_i[1] = 0.0;

    for j = 1:(snap.twojmax) 
        jju = snap.prebuilt_arrays.idxu_block[j+1]+1;
        jjup = snap.prebuilt_arrays.idxu_block[j]+1;

        # fill in left side of matrix layer from previous layer
        mb = 0
        while 2*mb < j
            ulist_r[jju] = 0.0;
            ulist_i[jju] = 0.0;

            for ma = 0:j-1 
                rootpq = snap.prebuilt_arrays.rootpqarray[Int(j - ma+1), Int(j - mb+1)];
                ulist_r[jju] +=
                    rootpq *
                    (a_r * ulist_r[jjup] +
                    a_i * ulist_i[jjup]);
                ulist_i[jju] +=
                    rootpq *
                    (a_r * ulist_i[jjup] -
                    a_i * ulist_r[jjup]);

                rootpq = snap.prebuilt_arrays.rootpqarray[Int(ma + 1+1), Int(j - mb+1)];
                ulist_r[jju+1] =
                    -rootpq *
                    (b_r * ulist_r[jjup] +
                    b_i * ulist_i[jjup]);
                ulist_i[jju+1] =
                    -rootpq *
                    (b_r * ulist_i[jjup] -
                    b_i * ulist_r[jjup]);
                jju+=1;
                jjup+=1;
            end
            jju+=1;
            mb +=1;
        end
 

        #  // copy left side to right side with inversion symmetry VMK 4.4(2)
        #  // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

        jju = snap.prebuilt_arrays.idxu_block[j+1];
        jjup = jju+(j+1)*(j+1) - 1;
        mbpar = 1;
        mb = 0;
        while 2*mb < j
            mapar = mbpar;
            for ma = 0:j 
                if (mapar == 1) 
                    ulist_r[jjup+1] = ulist_r[jju+1];
                    ulist_i[jjup+1] = -ulist_i[jju+1];
                else 
                    ulist_r[jjup+1] = -ulist_r[jju+1];
                    ulist_i[jjup+1] = ulist_i[jju+1];
                end
                mapar = -mapar;
                jju+=1;
                jjup-=1;
            end
            mbpar = -mbpar;
            mb += 1;
        end
    end
    runtime_arrays.ulist_r_ij[jj, :] = ulist_r
    runtime_arrays.ulist_i_ij[jj, :] = ulist_i
end
        
function add_uarraytot(r::AbstractFloat, rcut::AbstractFloat, wj::AbstractFloat,
                        jj::Int, jelem :: Int,
                        snap::SNAPParams,
                        runtime_arrays::RuntimeArrays)
    sfac = compute_sfac(snap, r, rcut);
    sfac *= wj;
    ulist_r = runtime_arrays.ulist_r_ij[jj, :];
    ulist_i = runtime_arrays.ulist_i_ij[jj, :];

    for j = 0:snap.twojmax 
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        for mb = 0:j 
            for ma = 0:j
                runtime_arrays.ulisttot_r[(jelem-1)*snap.prebuilt_arrays.idxu_max+jju+1] +=
                    sfac * ulist_r[jju+1];
                runtime_arrays.ulisttot_i[(jelem-1)*snap.prebuilt_arrays.idxu_max+jju+1] +=
                    sfac * ulist_i[jju+1];
                jju+=1
            end
        end
    end
end
      
  

        

# /* ----------------------------------------------------------------------
# Compute derivatives of Wigner U-functions for one neighbor
# see comments in compute_uarray()
# ------------------------------------------------------------------------- */

function compute_du(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                       z0::AbstractFloat, r::AbstractFloat, dz0dr::AbstractFloat,
                       wj::AbstractFloat, rcut::AbstractFloat, jj::Int, jelem::Int,
                       snap::SNAPParams, runtime_arrays::RuntimeArrays)

    rinv = 1.0/r;
    ux   = x * rinv;
    uy   = y * rinv;
    uz   = z * rinv;

    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r   = z0 * r0inv;
    a_i   = -z * r0inv;
    b_r   = y * r0inv;
    b_i   = -x * r0inv;

    dr0invdr = -r0inv*r0inv*r0inv * (r + z0 * dz0dr);

    dr0inv = [dr0invdr * ux, dr0invdr * uy, dr0invdr * uz];
    dz0    = [dz0dr * ux, dz0dr * uy, dz0dr * uz];

    da_r   = [dz0[k] * r0inv + z0 * dr0inv[k] for k in range(3)]
    da_i   = [-z * dr0inv[k] for k in range(3)]

    da_i[3] += -r0inv;

    db_r   = [ y * dr0inv[k] for k in range(3)]
    db_i   = [ -x * dr0inv[k] for k in range(3)]

    db_i[1] += -r0inv;
    db_r[2] += r0inv;

    ulist_r = runtime_arrays.ulist_r_ij[jj];
    ulist_i = runtime_arrays.ulist_i_ij[jj];

    for j = 1:snap.twojmax 
        jju = idxu_block[j+1];
        jjup = idxu_block[j];
        mb = 0
        while 2*mb < j
            
            for ma = 0:(j-1) 
                rootpq = snap.prebuilt_arrays.rootpqarray[j - ma+1][j - mb+1];
                for k = 1:3 
                    runtime_arrays.dulist_r[jju+1, k] +=
                        rootpq * (da_r[k] * ulist_r[jjup+1] +
                        da_i[k] * ulist_i[jjup+1] +
                        a_r * runtime_arrays.dulist_r[jjup+1, k] +
                        a_i * runtime_arrays.dulist_i[jjup+1, k]);
                    runtime_arrays.dulist_i[jju+1, k] +=
                        rootpq * (da_r[k] * ulist_i[jjup+1] -
                        da_i[k] * ulist_r[jjup+1] +
                        a_r * runtime_arrays.dulist_i[jjup+1, k] -
                        a_i * runtime_arrays.dulist_r[jjup+1, k]);
                end
     
                rootpq = snap.prebuilt_arrays.rootpqarray[ma + 1 + 1][j - mb + 1];
                for k = 1:3 
                    runtime_arrays.dulist_r[jju+1+1, k] =
                        -rootpq * (db_r[k] * ulist_r[jjup+1] +
                        db_i[k] * ulist_i[jjup+1] +
                        b_r * runtime_arrays.dulist_r[jjup+1, k] +
                        b_i * runtime_arrays.dulist_i[jjup+1, k]);
                    runtime_arrays.dulist_i[jju+1+1, k] =
                        -rootpq * (db_r[k] * ulist_i[jjup+1] -
                        db_i[k] * ulist_r[jjup+1] +
                        b_r * runtime_arrays.dulist_i[jjup+1, k] -
                        b_i * runtime_arrays.dulist_r[jjup+1, k]);
                end
     
                jju+=1;
                jjup+=1;
            end
            jju+=1
            mb +=1;
        end
    
 

        #  // copy left side to right side with inversion symmetry VMK 4.4(2)
        #  // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

        jju = snap.prebuilt_arrays.idxu_block[j+1];
        jjup = jju+(j+1)*(j+1)-1;
        mbpar = 1;
        mb = 0
        while 2*mb < j
            mapar = mbpar;
            for ma = 0:j 
                if (mapar == 1) 
                    for k = 1:3
                        runtime_arrays.dulist_r[jjup+1, k] = runtime_arrays.dulist_r[jju+1, k];
                        runtime_arrays.dulist_i[jjup+1, k] = -runtime_arrays.dulist_i[jju+1, k];
                    end
                else 
                    for k = 1:3 
                        runtime_arrays.dulist_r[jjup+1, k] = -runtime_arrays.dulist_r[jju+1, k];
                        runtime_arrays.dulist_i[jjup+1, k] = runtime_arrays.dulist_i[jju+1, k];
                    end
                end
                mapar = -mapar;
                jju+=1;
                jjup-=1;
            end
            mbpar = -mbpar;
            mb +=1;
        end
    end

    sfac = compute_sfac(r, rcut);
    dsfac = compute_dsfac(r, rcut);

    sfac *= wj;
    dsfac *= wj;
    for j = 0:snap.twojmax 
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        mb = 0;
        while 2*mb < j
            for ma = 0:j 
                runtime_arrays.dulist_r[jju+1][1] = dsfac * ulist_r[jju+1] * ux +
                                    sfac * runtime_arrays.dulist_r[jju+1, 1];
                runtime_arrays.dulist_i[jju+1][1] = dsfac * ulist_i[jju+1] * ux +
                                    sfac * runtime_arrays.dulist_i[jju_1, 1];
                runtime_arrays.dulist_r[jju+1][2] = dsfac * ulist_r[jju+1] * uy +
                                    sfac * runtime_arrays.dulist_r[jju+1, 2];
                runtime_arrays.dulist_i[jju+1][2] = dsfac * ulist_i[jju+1] * uy +
                                    sfac * runtime_arrays.dulist_i[jju+1, 2];
                runtime_arrays.dulist_r[jju+1][3] = dsfac * ulist_r[jju+1] * uz +
                                    sfac * runtime_arrays.dulist_r[jju+1, 3];
                runtime_arrays.dulist_i[jju+1][3] = dsfac * ulist_i[jju+1] * uz +
                                    sfac * runtime_arrays.dulist_i[jju+1, 3];
                jju++1;
            end
            mb +=1;
        end
    end
end

function compute_duidrj(rij::Vector{<:AbstractFloat}, wj::AbstractFloat, 
                        rcut::AbstractFloat, jj::Int, 
                        jelem::Int, rmin0::AbstractFloat,
                        snap::SNAPParams, runtime_arrays::RuntimeArrays)
    x = rij[1];
    y = rij[2];
    z = rij[3];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);
    rscale0 = rfac0 * pi / (rcut - snap.rmin0);
    theta0 = (r - snap.rmin0) * rscale0;
    cs = cos(theta0);
    sn = sin(theta0);
    z0 = r * cs / sn;
    dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

    compute_duarray(x, y, z,
                    z0, r, dz0dr,
                    wj, rcut, jj, jelem,
                    snap, runtime_arrays)
end