#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function readLAMMPS(tmp)
    
config = CONFIGStruct();        

config.dim = 3;
config.ncq = 0;
config.ncv = 0;
config.nce = 1;
config.ncf = 3;
config.ncx = 3;


k = 1;
sz1, sz2 = size(tmp);
while (k < sz1)
    n = tmp[k,1]; # number of atoms
    if k ==1
        config.natom = reshape([n],1,1);
    else
        config.natom = [config.natom n];
    end
    b = transpose(reshape(tmp[(k+1):(k + n),:], n, 7));
    if k==1
        config.t = reshape(b[1,:],1,n); 
        config.x = b[2:4,:]; 
        config.f = b[5:7,:]; 
    else
        config.t = [config.t reshape(b[1,:],1,n)]; 
        config.x = [config.x b[2:4,:]]; 
        config.f = [config.f b[5:7,:]]; 
    end
    b = tmp[(k+n+1):(k+n+2),:];
    if k == 1
        config.e = reshape([b[1,1]],1,1);
        config.a = reshape([b[2,1] 0 0],3,1);
        config.b = reshape([0 b[2,2] 0],3,1);
        config.c = reshape([0 0 b[2,3]],3,1);    
    else
        config.e = [config.e reshape([b[1,1]],1,1)];
        config.a = [config.a [b[2,1] 0 0]'];
        config.b = [config.b [0 b[2,2] 0]'];
        config.c = [config.c [0 0 b[2,3]]'];
    end
    k = k + n + 3;
end

config.q = reshape([],0,2);
config.v = reshape([],0,2);

return config 

end
