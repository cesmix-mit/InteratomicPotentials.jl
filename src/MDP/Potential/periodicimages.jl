#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function periodicimages(pbc::Array{Int64,1}, a::Array{Float64,1}, b::Array{Float64,1}, c=nothing)
    # a = (a1, a2, a3), b = (b1,b2,b3), c=(c1,c2,c3) are three vectors defining the simulation box
    # r is the distance from the original parallelepiped to the bounding parallelepiped 
    
    if c === nothing
        a = a[:];
        b = b[:];
        nd = length(a);
        px = zeros(nd,3^nd);
    
        px[:,1] = [0.0; 0.0];
        px[:,2] = a;
        px[:,3] = -a;
        px[:,4] = b;
        px[:,5] = -b;
        px[:,6] = a+b;
        px[:,7] = a-b;
        px[:,8] = -a+b;
        px[:,9] = -a-b;
    
        if (pbc[1]==0) && (pbc[2]==0)
            indx = 1;    
        elseif (pbc[1]==1) && (pbc[2]==0)
            indx = [1 2 3];    
        elseif (pbc[1]==0) && (pbc[2]==1)
            indx = [1 4 5];       
        elseif (pbc[1]==1) && (pbc[2]==1)
            indx = [1 2 3 4 5 6 7 8 9];        
        end
        if length(indx)==1
            ximages = px[:,indx];    
        else
            ximages = px[:,indx[:]];    
        end
    else    
        a = a[:];
        b = b[:];
        c = c[:];
        nd = length(a);
        px = zeros(nd,3^nd);
    
        px[:,1] = [0.0; 0.0; 0.0];
        px[:,2] = a;
        px[:,3] = -a;
        px[:,4] = b;
        px[:,5] = -b;
        px[:,6] = a+b;
        px[:,7] = a-b;
        px[:,8] = -a+b;
        px[:,9] = -a-b;
    
        px[:,10] = c;
        px[:,11] = c+a;
        px[:,12] = c-a;
        px[:,13] = c+b;
        px[:,14] = c-b;
        px[:,15] = c+a+b;
        px[:,16] = c+a-b;
        px[:,17] = c-a+b;
        px[:,18] = c-a-b;
    
        px[:,19] = -c;
        px[:,20] = a-c;
        px[:,21] = -a-c;
        px[:,22] = b-c;
        px[:,23] = -b-c;
        px[:,24] = a+b-c;
        px[:,25] = a-b-c;
        px[:,26] = -a+b-c;
        px[:,27] = -a-b-c;
    
        if (pbc[1]==0) && (pbc[2]==0) && (pbc[3]==0)
            indx = 1;    
        elseif (pbc[1]==1) && (pbc[2]==0) && (pbc[3]==0)
            indx = [1 2 3];    
        elseif (pbc[1]==0) && (pbc[2]==1) && (pbc[3]==0)
            indx = [1 4 5];       
        elseif (pbc[1]==0) && (pbc[2]==0) && (pbc[3]==1)
            indx = [1 10 19];           
        elseif (pbc[1]==1) && (pbc[2]==1) && (pbc[3]==0)
            indx = reshape(collect(1:9), 1, 9);        
        elseif (pbc[1]==1) && (pbc[2]==0) && (pbc[3]==1)
            indx = ones(1,9) + [0 1 2 9 10 11 18 19 20];            
        elseif (pbc[1]==0) && (pbc[2]==1) && (pbc[3]==1)
            indx = ones(1,9) + [0 3 4 9 12 13 18 21 22];                
        elseif (pbc[1]==1) && (pbc[2]==1) && (pbc[3]==1)
            indx = reshape(collect(1:27), 1, 27);        
        end    
        if length(indx)==1
            ximages = px[:,indx];    
        else
            ximages = px[:,indx[:]];    
        end

        return ximages;
    end
    
    end
    
    
    