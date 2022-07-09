function fullneighborlist(x, a, b, c, pbc, rcutmax)

# offset simulation box for periodic boundary conditions
boxoffset = [rcutmax; rcutmax; rcutmax]; 

# bounding box
B2R, R2B, v, w, vr, wr, pimages = prebox(pbc[:], boxoffset, a[:], b[:], c[:]);

y, alist = atomlist(x, pimages, wr, B2R);

neighi, neighj, neighnum = neighborlist(y, rcutmax);

n = size(x,2);
neighi = neighi[neighi .<= n];
ne = length(neighi);
neighlist = neighj[1:ne];
neighnum = neighnum[1:n];
neighnum = [0; cumsum(neighnum)];

return y, alist, neighlist, neighnum

end

