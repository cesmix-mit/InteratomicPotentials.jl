
function neighpairs(x, q, pairlist, pairnum, atomtype, ilist, alist)

dim = size(x,1);
ncq = size(q,1);
inum = length(ilist);
ijnum = sum(pairnum[ilist.+1]-pairnum[ilist]);


ai = zeros(eltype(ilist),ijnum);
aj = zeros(eltype(ilist),ijnum);
ti = zeros(eltype(ilist),ijnum);
tj = zeros(eltype(ilist),ijnum);
xij = zeros(dim,ijnum);
qi = zeros(ncq,ijnum);
qj = zeros(ncq,ijnum);

count = 0;
for ii=1:inum
    i = ilist[ii]; # atom i
    n1 = pairnum[i];
    n2 = pairnum[i+1];    
    n = n2-n1; 
    g = pairlist[(n1+1):n2]; # a list of neighbors around atom i    
    j = alist[g];    
    ind = (count+1):(count+n);
    count = count + n;
        
    ai[ind] .= i;
    aj[ind] .= j;    
    ti[ind] .= atomtype[i];
    tj[ind] .= atomtype[j];
    xij[:,ind] = x[:,g] .- x[:,i];
    
    if ncq>0
        qi[ind] .= q[i];
        qj[ind] .= q[j];
    end
end

return xij, ai, aj, ti, tj, qi, qj

end


function neighpairs(x, pairlist, pairnum, atomtype, ilist, alist)

    dim = size(x,1);
    inum = length(ilist);
    ijnum = sum(pairnum[ilist.+1]-pairnum[ilist]);
        
    ai = zeros(eltype(ilist),ijnum);
    aj = zeros(eltype(ilist),ijnum);
    ti = zeros(eltype(ilist),ijnum);
    tj = zeros(eltype(ilist),ijnum);
    xij = zeros(dim,ijnum);
    
    count = 0;
    for ii=1:inum
        i = ilist[ii]; # atom i
        n1 = pairnum[i];
        n2 = pairnum[i+1];    
        n = n2-n1; 
        g = pairlist[(n1+1):n2]; # a list of neighbors around atom i    
        j = alist[g];    
        ind = (count+1):(count+n);
        count = count + n;
            
        ai[ind] .= i;
        aj[ind] .= j;    
        ti[ind] .= atomtype[i];
        tj[ind] .= atomtype[j];
        xij[:,ind] = x[:,g] .- x[:,i];
        
    end
    
    return xij, ai, aj, ti, tj
    
end
    
    
        
        
