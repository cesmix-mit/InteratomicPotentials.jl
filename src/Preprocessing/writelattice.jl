#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function writelattice(lat,filename)

    tmp = [lat.style size(lat.basis,1) lat.spaceflag lat.scale];
    basis = transpose(lat.basis);

    nsize = zeros(11,1);
    nsize[1] = length(tmp);
    nsize[2] = length(lat.origin[:]);
    nsize[3] = length(lat.orientx[:]);
    nsize[4] = length(lat.orienty[:]);
    nsize[5] = length(lat.orientz[:]);
    nsize[6] = length(lat.spacing[:]);
    nsize[7] = length(lat.a1[:]);
    nsize[8] = length(lat.a2[:]);
    nsize[9] = length(lat.a3[:]);
    nsize[10] = length(basis[:]);
    nsize[11] = length(lat.type[:]);

    fileID = open(filename,"w");
    write(fileID,Float64(length(nsize[:])));
    write(fileID,Float64.(nsize[:]));
    write(fileID,Float64.(tmp[:]));
    write(fileID,Float64.(lat.origin[:]));
    write(fileID,Float64.(lat.orientx[:]));
    write(fileID,Float64.(lat.orienty[:]));
    write(fileID,Float64.(lat.orientz[:]));
    write(fileID,Float64.(lat.spacing[:]));
    write(fileID,Float64.(lat.a1[:]));
    write(fileID,Float64.(lat.a2[:]));
    write(fileID,Float64.(lat.a3[:]));
    write(fileID,Float64.(basis[:]));
    write(fileID,Float64.(lat.type[:]));
    close(fileID);

end
