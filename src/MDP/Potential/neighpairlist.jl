function neighpairlist(x, ilist, neighlist, neighnum, rcutsq)

inum = length(ilist);

pairnum = zeros(eltype(ilist), inum);
pairlist = [];
for ii = 1:inum
    i = ilist[ii]; # atom i
    n1 = neighnum[i];
    n2 = neighnum[i+1];    

    g = neighlist[(n1+1):n2]; # a list of neighbors around atom i
    d = x[:,g] .- x[:,i]; # distance between neighbors and atom i
    dsq = sum(d.*d, dims=1);
    dsq = dsq[:];
    ind = findall(dsq .<= rcutsq); # find neighbors within rcut from atom i
    # a list of neighbors within rcut from atom i
    pairnum[ii] = length(ind);       
    pairlist = [pairlist; g[ind]];      
end

pairnum = [0; cumsum(pairnum)];
pairlist = eltype(ilist).(pairlist)

return pairlist, pairnum

end

function neighpairlist(x, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)

#typej = unique(typej);
inum = length(ilist);

pairnum = zeros(eltype(ilist), inum);
pairlist = [];
for ii = 1:inum
    i = ilist[ii]; # atom i
    n1 = neighnum[i];
    n2 = neighnum[i+1];    

    # cut-off radius between atom of type i and atom of type j
    itype = atomtype[i];
    g = neighlist[(n1+1):n2]; # a list of neighbors around atom i            
    j = alist[g];    
    jtype = atomtype[j]; # types of neighbors around atom i
    rsq = rcutsq[jtype .+ (itype-1)*(ntypes)];

    d = x[:,g] .- x[:,i]; # distance between neighbors and atom i
    dsq = sum(d.*d, dims=1);
    dsq = dsq[:]
    ind = findall(dsq .<= rsq); # find neighbors within rcut from atom i        
        
    # a list of neighbors with typej within rcut from atom i     
    pairnum[ii] = length(ind);       
    pairlist = [pairlist; g[ind]];      
end

pairnum = [0; cumsum(pairnum)];
pairlist = eltype(ilist).(pairlist)

return pairlist, pairnum

end
    
function neighpairlisttypej(x, ilist, alist, atomtype, typej, neighlist, neighnum, rcutsq)

typej = unique(typej);
inum = length(ilist);

pairnum = zeros(eltype(ilist), inum);
pairlist = [];
for ii = 1:inum
    i = ilist[ii]; # atom i
    n1 = neighnum[i];
    n2 = neighnum[i+1];    
    g = neighlist[(n1+1):n2]; # a list of neighbors around atom i    
    j = alist[g];
    atype = atomtype[j]; # types of neighbors around atom i
    ind1 = findall(typej .== atype);  # choose neighbors of typej
    
    d = x[:,g] - x[:,i]; # distance between neighbors and atom i
    dsq = sum(d.*d, dims=1);
    dsq = dsq[:]
    ind2 = findall(dsq .<= rcutsq); # find neighbors within rcut from atom i        
    
    # a list of neighbors with typej within rcut from atom i 
    ind = intersect(ind1, ind2);
    pairnum[ii] = length(ind);       
    pairlist = [pairlist; g[ind]];      
end
pairnum = [0; cumsum(pairnum)];

pairlist = eltype(ilist).(pairlist)

return pairlist, pairnum

end




