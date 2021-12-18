function tallysingle(ei, fi, ai)

# energies per atom i
e = accumarray(ai, ei);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fi,1);

# forces
f = zeros(dim,N);
for d = 1:dim
    # tally forces on each dimension
    f[d,:] = reshape(accumarray(ai, fi[d,:]),(1, N)); 
end

return e, f

end

function tallypair(eij, fij, ai, aj)

# energies per atom i
e = accumarray(ai, eij);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
# f = zeros(dim,N);
# for d = 1:dim
#     # tally forces on each dimension
#     f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
# end

return e, f


end

function tallytriplet(eijk, fij, fik, ai, aj, ak)

# energies per atom i
e = accumarray(ai, eijk);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
f = f + accumarray2d(ai, fik) - accumarray2d(ak, fik);
# f = zeros(dim,N);
# for d = 1:dim
#     # tally forces on each dimension
#     f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
#     f[d,:] = reshape(accumarray(ai, fik[d,:]),(1, N)) - reshape(accumarray(ak, fik[d,:]),(1, N)); 
# end

return e, f

end

function tallyquadlet(eijkl, fij, fik, fil, ai, aj, ak, al)

# energies per atom i
e = accumarray(ai, eijkl);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
f = f + accumarray2d(ai, fik) - accumarray2d(ak, fik);
f = f + accumarray2d(ai, fil) - accumarray2d(al, fil);
# f = zeros(dim,N);
# for d = 1:dim
#     # tally forces on each dimension
#     f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
#     f[d,:] = reshape(accumarray(ai, fik[d,:]),(1, N)) - reshape(accumarray(ak, fik[d,:]),(1, N)); 
#     f[d,:] = reshape(accumarray(ai, fil[d,:]),(1, N)) - reshape(accumarray(al, fil[d,:]),(1, N)); 
# end

return e, f

end



















