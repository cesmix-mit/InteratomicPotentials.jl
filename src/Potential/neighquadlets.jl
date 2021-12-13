function neighquadlets(x, q, quadletlist, quadletnum, atomtype, ilist, alist)

    dim = size(x,1);
    ncq = size(q,1);
    inum = length(ilist);
    ijklnum = sum(quadletnum(ilist+1)-quadletnum(ilist));
    
    ai = zeros(Int64,ijklnum);
    aj = zeros(Int64,ijklnum);
    ak = zeros(Int64,ijklnum);
    al = zeros(Int64,ijklnum);
    ti = zeros(Int64,ijklnum);
    tj = zeros(Int64,ijklnum);
    tk = zeros(Int64,ijklnum);
    tl = zeros(Int64,ijklnum);
    xij = zeros(dim,ijklnum);
    xik = zeros(dim,ijklnum);
    xil = zeros(dim,ijklnum);
    qi = zeros(ncq,ijklnum);
    qj = zeros(ncq,ijklnum);
    qk = zeros(ncq,ijklnum);
    ql = zeros(ncq,ijklnum);
    
    count = 0;
    for ii=1:inum
      i = ilist[ii]; # atom i
      n1 = quadletnum[i];
      n2 = quadletnum[i+1];    
      n = n2-n1; 
      gj = quadletlist((n1+1):n2,1); # a list of neighbors around atom i    
      j = alist[gj];    
      gk = quadletlist((n1+1):n2,2); # a list of neighbors around atom i    
      k = alist[gk];    
      gl = quadletlist((n1+1):n2,3); # a list of neighbors around atom i    
      l = alist[gl];    
      
      ind = (count+1):(count+n);
      count = count + n;
      
      ai[ind] = i;
      aj[ind] = j;
      ak[ind] = k;
      al[ind] = l;
      ti[ind] = atomtype[i];
      tj[ind] = atomtype[j];
      tk[ind] = atomtype[k];
      tl[ind] = atomtype[l];
      xij[:,ind] = x[:,gj] - x[:,i];
      xik[:,ind] = x[:,gk] - x[:,i];
      xil[:,ind] = x[:,gl] - x[:,i];
      if ncq>0
          qi[ind] = q[i];
          qj[ind] = q[j];
          qk[ind] = q[k];
          ql[ind] = q[l];
      end
    end
    
    return xij, xik, xil, ai, aj, ak, al, ti, tj, tk, tl, qi, qj, qk, ql
    
    end
    
    