
# ----------------------------------------------------------------------
# compute Wigner U-functions for one neighbor
# ------------------------------------------------------------------------- */
function compute_ui(
    ielem::Int,
    snap::SNAPParams, 
    runtime_arrays :: RuntimeArrays)

    if ~snap.chem_flag
        ielem = 1
    end



    zero_uarraytot(ielem, snap, runtime_arrays)

    num_of_interactions = length(runtime_arrays.indij)
    for ind = 1:num_of_interactions
        jj = runtime_arrays.indij[ind][2]
        rij = runtime_arrays.rij[ind]
        rcutij = runtime_arrays.rcutij[ind]
        wj = runtime_arrays.wj[ind]
        elementj = runtime_arrays.element[ind]

        x = rij[1]
        y = rij[2]
        z = rij[3]
        rsq = x*x + y*y + z*z 
        r   = sqrt(rsq) 
        theta0 = (r - snap.rmin0) * snap.rfac0 * pi / (rcutij - snap.rmin0) 
        z0 = r / tan(theta0)
        compute_uarray(x, y, z, z0, r, ind+1, 
                snap, 
                runtime_arrays)
        
        if (snap.chem_flag)
            add_uarraytot(r, rcutij, wj,
                    ind+1, elementj, snap, runtime_arrays)
        else
            add_uarraytot(r, rcutij, wj, ind+1, 1, snap, runtime_arrays)
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
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        jjup = snap.prebuilt_arrays.idxu_block[j];

        # fill in left side of matrix layer from previous layer
        mb = 0
        while 2*mb <= j
            ulist_r[jju+1] = 0.0;
            ulist_i[jju+1] = 0.0;

            for ma = 0:(j-1) 
                rootpq = snap.prebuilt_arrays.rootpqarray[j - ma+1, j - mb+1];
                ulist_r[jju+1] +=
                    rootpq *
                    (a_r * ulist_r[jjup+1] +
                    a_i * ulist_i[jjup+1]);
                ulist_i[jju+1] +=
                    rootpq *
                    (a_r * ulist_i[jjup+1] -
                    a_i * ulist_r[jjup+1]);
                rootpq = snap.prebuilt_arrays.rootpqarray[ma + 1+1, j - mb+1];
                ulist_r[jju+1+1] =
                        -rootpq *
                        (b_r * ulist_r[jjup+1] +
                        b_i * ulist_i[jjup+1]);
                ulist_i[jju+1+1] =
                        -rootpq *
                        (b_r * ulist_i[jjup+1] -
                        b_i * ulist_r[jjup+1]);
                jju+=1;
                jjup+=1;
            end
            jju+=1;
            mb +=1;
        end


        #  // copy left side to right side with inversion symmetry VMK 4.4(2)
        #  // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        jjup = jju+(j+1)*(j+1) - 1 ;
        mbpar = 1;
        mb = 0;
        while 2*mb <= j
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

function zero_uarraytot(ielem::Int, snap::SNAPParams, runtime_arrays::RuntimeArrays)
    num_elements = snap.prebuilt_arrays.num_elements
    for jelem = 0:(num_elements-1)
        for j = 0:snap.twojmax 
            jju = snap.prebuilt_arrays.idxu_block[j+1];
            for mb = 0:j 
                for ma = 0:j 
                    runtime_arrays.ulisttot_r[jelem*snap.prebuilt_arrays.idxu_max+jju+1] = 0.0;
                    runtime_arrays.ulisttot_i[jelem*snap.prebuilt_arrays.idxu_max+jju+1] = 0.0;

                    # // utot(j,ma,ma) = wself, sometimes
                    if ((jelem+1) == ielem) || (snap.wselfall_flag)
                        if (ma==mb)
                            runtime_arrays.ulisttot_r[jelem*snap.prebuilt_arrays.idxu_max+jju+1] = snap.prebuilt_arrays.wself; # ///// double check this
                        end
                    end
                    jju+=1;
                end
            end
        end
    end
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

function compute_duarray(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
     z0::AbstractFloat, r::AbstractFloat, dz0dr::AbstractFloat,
     wj::AbstractFloat, rcut::AbstractFloat, jj::Int,
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

    da_r   = [dz0[k] * r0inv + z0 * dr0inv[k] for k = 1:3]
    da_i   = [-z * dr0inv[k] for k = 1:3]

    da_i[3] += -r0inv;

    db_r   = [ y * dr0inv[k] for k = 1:3]
    db_i   = [ -x * dr0inv[k] for k = 1:3]

    db_i[1] += -r0inv;
    db_r[2] += r0inv;

    ulist_r = runtime_arrays.ulist_r_ij[jj, :];
    ulist_i = runtime_arrays.ulist_i_ij[jj, :];

    runtime_arrays.dulist_r[1, :] = [0.0, 0.0, 0.0]
    runtime_arrays.dulist_i[1, :] = [0.0, 0.0, 0.0]

    for j = 1:snap.twojmax 
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        jjup = snap.prebuilt_arrays.idxu_block[j];
        mb = 0
        while 2*mb <= j

            runtime_arrays.dulist_r[jju+1, :] = [0.0, 0.0, 0.0]
            runtime_arrays.dulist_i[jju+1, :] = [0.0, 0.0, 0.0]

            for ma = 0:(j-1) 
                rootpq = snap.prebuilt_arrays.rootpqarray[j - ma+1, j - mb+1];
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

                rootpq = snap.prebuilt_arrays.rootpqarray[ma + 1 + 1,j - mb + 1];
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
        while 2*mb <= j
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

    sfac = compute_sfac(snap, r, rcut);
    dsfac = compute_dsfac(snap, r, rcut);
    sfac *= wj;
    dsfac *= wj;
    for j = 0:snap.twojmax 
        jju = snap.prebuilt_arrays.idxu_block[j+1];
        mb = 0;
        while 2*mb <= j
            for ma = 0:j 
                runtime_arrays.dulist_r[jju+1,1] = dsfac * ulist_r[jju+1] * ux +
                                sfac * runtime_arrays.dulist_r[jju+1, 1];
                runtime_arrays.dulist_i[jju+1,1] = dsfac * ulist_i[jju+1] * ux +
                                sfac * runtime_arrays.dulist_i[jju+1, 1];
                runtime_arrays.dulist_r[jju+1,2] = dsfac * ulist_r[jju+1] * uy +
                                sfac * runtime_arrays.dulist_r[jju+1, 2];
                runtime_arrays.dulist_i[jju+1,2] = dsfac * ulist_i[jju+1] * uy +
                                sfac * runtime_arrays.dulist_i[jju+1, 2];
                runtime_arrays.dulist_r[jju+1,3] = dsfac * ulist_r[jju+1] * uz +
                                sfac * runtime_arrays.dulist_r[jju+1, 3];
                runtime_arrays.dulist_i[jju+1,3] = dsfac * ulist_i[jju+1] * uz +
                                sfac * runtime_arrays.dulist_i[jju+1, 3];
                jju+=1;
            end
            mb +=1;
        end
    end
end

function compute_duidrj(rij::SVector{3, <:AbstractFloat}, wj::AbstractFloat, 
        rcut::AbstractFloat, jj::Int, 
        snap::SNAPParams, runtime_arrays::RuntimeArrays)
    x = rij[1];
    y = rij[2];
    z = rij[3];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);
    rscale0 = snap.rfac0 * pi / (rcut - snap.rmin0);
    theta0 = (r - snap.rmin0) * rscale0;
    cs = cos(theta0);
    sn = sin(theta0);
    z0 = r * cs / sn;
    dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

    compute_duarray(x, y, z,
        z0, r, dz0dr,
        wj, rcut, jj, 
        snap, runtime_arrays)
end


function compute_sfac(snap::SNAPParams, r::AbstractFloat, rcut::AbstractFloat)
    if ~snap.switch_flag
        sfac = 1.0 
    else
        if r <= snap.rmin0
            sfac = 1.0 
        elseif r > rcut
            sfac =  0.0
        else
            rcutfac = pi / (rcut - snap.rmin0)
            sfac = 0.5 * ( cos( (r - snap.rmin0) * rcutfac ) + 1.0 )
        end
    end
    return sfac
end

function compute_dsfac(snap::SNAPParams, r::AbstractFloat, rcut::AbstractFloat)
    if ~snap.switch_flag
        dsfac = 0.0
    else
        if r <= snap.rmin0
            dsfac = 0.0 
        elseif rcut > rcut
            dsfac = 0.0
        else
            rcutfac = pi / (rcut - snap.rmin0)
            dsfac = -0.5 * ( sin( (r - snap.rmin0) * rcutfac ) * rcutfac )
        end
    end
    return dsfac
end