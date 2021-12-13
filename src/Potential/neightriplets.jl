function  neightriplets(x, q, tripletlist, tripletnum, atomtype, ilist, alist)

    dim = size(x,1);
    ncq = size(q,1);
    inum = length(ilist);
    ijknum = sum(tripletnum[ilist.+1]-tripletnum[ilist]);
    
    ai = zeros(Int64,ijknum);
    aj = zeros(Int64,ijknum);
    ak = zeros(Int64,ijknum);
    ti = zeros(Int64,ijknum);
    tj = zeros(Int64,ijknum);
    tk = zeros(Int64,ijknum);
    xij = zeros(dim,ijknum);
    xik = zeros(dim,ijknum);
    qi = zeros(ncq,ijknum);
    qj = zeros(ncq,ijknum);
    qk = zeros(ncq,ijknum);
    
    count = 0;
    for ii=1:inum
        i = ilist[ii]; # atom i
        n1 = tripletnum[i];
        n2 = tripletnum[i+1];        
        n = n2-n1; 
        gj = tripletlist[(n1+1):n2,1]; # a list of neighbors around atom i    
        j = alist[gj];    
        gk = tripletlist[(n1+1):n2,2]; # a list of neighbors around atom i    
        k = alist[gk];    
        
        ind = (count+1):(count+n);
        count = count + n;
        
        ai[ind] = i;
        aj[ind] = j;
        ak[ind] = k;    
        ti[ind] = atomtype[i];
        tj[ind] = atomtype[j];
        tk[ind] = atomtype[k];
        xij[:,ind] = x[:,gj] - x[:,i];
        xik[:,ind] = x[:,gk] - x[:,i];    
        if ncq>0
            qi[ind] = q[i];
            qj[ind] = q[j];
            qk[ind] = q[k];
        end
    end
    
    return xij, xik, ai, aj, ak, ti, tj, tk, qi, qj, qk
    
    end
    