
function neighborlist(r, rcut)

dim, N = size(r);

r1 = reshape(r,(dim, N, 1));
r2 = reshape(r,(dim, 1, N));
r12 = r1 .- r2;
d = reshape(sqrt.(sum(r12.*r12,dims=1)), (N,N)); # distances between atoms

A = (d .<= rcut);  # adjacency marix
index = 1:(N+1):N*N; # indices of diagonal entries
A[index] .= 0.0;  # set diagonal of A to zero

# pairs of atoms (i,j) within the cut-off radius rcut
ind = findall(A.==1);
ai = getindex.(ind, 2);
aj = getindex.(ind, 1);

# number of neighbors for each atom i
numneigh = accumarray(ai, ones(length(ai),1)); 

ai = Int64.(ai);
aj = Int64.(aj);
numneigh = Int64.(numneigh);

return ai, aj, numneigh

end


