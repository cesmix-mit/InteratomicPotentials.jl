#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function checkconfig(x, ximages, B2C, C2B)

nd = size(B2C,1);
n = size(x,2);

xx = B2C*x;
pp = B2C*ximages;

if nd ==2
    for i = 1:n
        xt = xx[:,i];
        if (0<=xt[1]) && (xt[1]<=1) && (0<=xt[2]) && (xt[2]<=1) 
        else        
            xp = xt + pp;        
            for j = 1:9
                xq = xp[:,j];
                if (0<=xq[1]) && (xq[1]<=1) && (0<=xq[2]) && (xq[2]<=1) 
                    xx[:,i] = xq;                
                    break;
                end
            end
        end
    end    
else    
    for i = 1:n
        xt = xx[:,i];
        if (0<=xt[1]) && (xt[1]<=1) && (0<=xt[2]) && (xt[2]<=1) && (0<=xt[3]) && (xt[3]<=1)
        else        
            xp = xt .+ pp;  # periodic images      
            for j = 1:27
                xq = xp[:,j];
                if (0<=xq[1]) && (xq[1]<=1) && (0<=xq[2]) && (xq[2]<=1) && (0<=xq[3]) && (xq[3]<=1)
                    xx[:,i] = xq;                
                    break;
                end
            end
        end
    end
end

x = C2B*xx;

return x

end
