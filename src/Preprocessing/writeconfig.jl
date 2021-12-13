#****************************************************************************
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project
#
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function writeconfig(config, filename)

# disp('Writing configurations into a binary file...');

nconfigs = config.nconfigs; # number of configurations
dim = config.dim; # physical dimension
ncx = config.ncx; # number of compoments of x
ncv = config.ncv; # number of compoments of v
nce = config.nce; # number of compoments of e
ncf = config.ncf; # number of compoments of f
ncq = config.ncq; # number of compoments of q

tmp = [nconfigs dim ncq ncv nce ncf];

nsize = zeros(20,1);
nsize[1] = length(tmp[:]);
nsize[2] = length(config.natom[:]);
nsize[3] = length(config.a[:]);
nsize[4] = length(config.b[:]);
nsize[5] = length(config.c[:]);
nsize[6] = length(config.e[:]);
nsize[7] = length(config.t[:]);
nsize[8] = length(config.x[:]);
nsize[9] = length(config.q[:]);
nsize[10] = length(config.v[:]);
nsize[11] = length(config.f[:]);
nsize[12] = length(config.we[:]);
nsize[13] = length(config.wf[:]);

fileID = open(filename,"w");

write(fileID,Float64(length(nsize[:])));
write(fileID,Float64.(nsize[:]));
write(fileID,Float64.(tmp[:]));
write(fileID,Float64.(config.natom));
write(fileID,Float64.(config.a));
write(fileID,Float64.(config.b));
write(fileID,Float64.(config.c));
write(fileID,Float64.(config.e));
write(fileID,Float64.(config.t));
write(fileID,Float64.(config.x));
write(fileID,Float64.(config.q));
write(fileID,Float64.(config.v));
write(fileID,Float64.(config.f));
write(fileID,Float64.(config.we));
write(fileID,Float64.(config.wf));
close(fileID);

end
