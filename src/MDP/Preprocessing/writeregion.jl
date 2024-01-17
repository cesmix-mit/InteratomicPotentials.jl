#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function writeregion(reg,filename)

    tmp = [reg.triclinic];

    nsize = zeros(4,1);
    nsize[1] = length(tmp);
    nsize[2] = length(reg.boxlo[:]);
    nsize[3] = length(reg.boxhi[:]);
    nsize[4] = length(reg.boxtilt[:]);

    fileID = open(filename,"w");
    write(fileID,Float64(length(nsize[:])));
    write(fileID,Float64.(nsize[:]));
    write(fileID,Float64.(tmp[:]));
    write(fileID,Float64.(reg.boxlo[:]));
    write(fileID,Float64.(reg.boxhi[:]));
    write(fileID,Float64.(reg.boxtilt[:]));
    close(fileID);

end
