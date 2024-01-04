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
        ## Note, unlike many other while loops, this is strictly less than (see original code)
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
            ma = 0
            # Also strictly <
            while ma < mb
                dudr_r = dulist_r[jju+1, :];
                dudr_i = dulist_i[jju+1, :];
                jjjmambyarray_r = ylist_r[(jelem-1)*idxu_max+jju+1];
                jjjmambyarray_i = ylist_i[(jelem-1)*idxu_max+jju+1];

                for k=1:3
                    dedr[k] +=
                        dudr_r[k] * jjjmambyarray_r +
                        dudr_i[k] * jjjmambyarray_i;
                    
                end
                jju+=1;
                ma+=1;
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



