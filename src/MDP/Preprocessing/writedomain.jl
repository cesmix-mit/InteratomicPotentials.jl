#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function writedomain(dom,filename)

    tmp = [dom.triclinic];

    nsize = zeros(6,1);
    nsize[1] = length(tmp);
    nsize[2] = length(dom.boxlo[:]);
    nsize[3] = length(dom.boxhi[:]);
    nsize[4] = length(dom.boxtilt[:]);
    nsize[5] = length(dom.pbc[:]);
    nsize[6] = length(dom.bcs[:]);

    fileID = open(filename,"w");
    write(fileID,Float64(length(nsize[:])));
    write(fileID,Float64.(nsize[:]));
    write(fileID,Float64.(tmp[:]));
    write(fileID,Float64.(dom.boxlo[:]));
    write(fileID,Float64.(dom.boxhi[:]));
    write(fileID,Float64.(dom.boxtilt[:]));
    write(fileID,Float64.(dom.pbc[:]));
    write(fileID,Float64.(dom.bcs[:]));    
    close(fileID);

end
