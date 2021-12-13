function empiricalpotential(x, q, t, a, b, c, pbc, rcutmax, eta, kappa, potential)

# compute full neighbor list
y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
if !isempty(q)
    q = q[:,alist];
end

dim, N = size(x);
e = zeros(N);
f = zeros(dim,N);    
for i = 1:length(potential)        
    
    style, bondtype, bondedatomtypes, potentialfunc, mu, rcut = getpotential(potential[i]);    
    potentialfunc = getfield(Main, Symbol(potentialfunc));
    
    # single potentials
    if (style == "single") & (bondtype == "nonbonded")            
        ilist = Array(1:N);                        
     elseif (style == "single") & (bondtype == "bonded")     
        typei = bondedatomtypes[1];
        ilist = findatomtype(Array(1:N), t, typei);      
    end
    if (style == "single")
        xi, qi, ai, ti = neighsingles(x, q, t, ilist);        
        ei, fi = potentialfunc(xi, qi, ti, mu, eta, kappa);
        ea, fa = tallysingle(ei, fi, ai);            
    end

    # pair potentials
    if (style == "pair") & (bondtype == "nonbonded")            
        ilist = Int64.(Array(1:N));                        
        pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);
     elseif (style == "pair") & (bondtype == "bonded")     
        typei = bondedatomtypes[1];
        typej = bondedatomtypes[2];            

        ilist = findatomtype(Array(1:N), t, typei);      
        pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, typej, neighlist, neighnum, rcut*rcut);                                                
    end
    if (style == "pair")
        xij, ai, aj, ti, tj, qi, qj = neighpairs(y, q, pairlist, pairnum, t, ilist, alist);            
        eij, fij = potentialfunc(xij, qi, qj, ti, tj, mu, eta, kappa);                
        ea, fa = tallypair(eij, fij, ai, aj);            
    end

    # triplet potentials
    if (style == "triplet") & (bondtype == "nonbonded")            
        ilist = Array(1:N);                     
        pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);        
     elseif (style == "triplet") & (bondtype == "bonded")     
        typei = bondedatomtypes[1];
        typej = bondedatomtypes[2];            
        typek = bondedatomtypes[3];            

        ilist = findatomtype(Array(1:N), t, typei);      
        pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
    end
    if (style == "triplet")
        tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
        xij, xik, ai, aj, ak, ti, tj, tk, qi, qj, qk = neightriplets(y, q, tripletlist, tripletnum, atomtype, ilist, alist);
        eij, fij, fik = potentialfunc(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);        
        ea, fa = tallytriplet(eij, fij, fik, ai, aj, ak);            
    end
    
    # quadlet potentials
    if (style == "quadlet") & (bondtype == "nonbonded")            
        ilist = Array(1:N);                        
        pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);        
     elseif (style == "quadlet") & (bondtype == "bonded")     
        typei = bondedatomtypes[1];
        typej = bondedatomtypes[2];            
        typek = bondedatomtypes[3];      
        typel = bondedatomtypes[4];      

        ilist = findatomtype(Array(1:N), t, typei);      
        pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej typek typel], neighlist, neighnum, rcut*rcut);                                                
    end
    if (style == "quadlet")
        quadletlist, quadletnum = neighquadletlist(pairlist, pairnumsum, ilist);        
        xij, xik, xil, ai, aj, ak, al, ti, tj, tk, tl, qi, qj, qk, ql = 
                    neighquadlets(x, q, quadletlist, quadletnum, atomtype, ilist, alist);
        eij, fij, fik, fil = potentialfunc(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa);
        ea, fa = tallyquadlet(eij, fij, fik, fil, ai, aj, ak, al);            
    end
        
    e = e .+ ea[:];
    f = f .+ fa;                    
end    

etot = sum(e);

return etot, e, f

end





