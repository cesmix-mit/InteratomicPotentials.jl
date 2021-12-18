function makejk(n)
    
m = (n-1)*n/2;
indj = zeros(Int64,m);
indk = zeros(Int64,m);
k1 = 1;
for i = 1:(n-1)
    ind = k1:(k1+(n-i)-1);
    indj[ind] = i;
    indk[ind] = (i+1):n;
    k1 = k1 + (n-i);
end

return indj, indk
    
end

function neightripletlist(pairlist, pairnumsum, ilist)

pairnum = pairnumsum[2:end]-pairnumsum[1:end-1];
tripletnum  = (pairnum-1).*pairnum/2;

inum = length(ilist);
tripletlistj = [];
tripletlistk = [];

for ii = 1:inum
    n1 = pairnumsum[ii];
    n2 = pairnumsum[ii+1];        
    g = pairlist[(n1+1):n2]; # a list of neighbors around atom i
    
    indj, indk = makejk(n);        
    tripletlistj = [tripletlistj; g[indj]];      
    tripletlistk = [tripletlistk; g[indk]];       
end

tripletlist = [tripletlistj[:] tripletlistk[:]];
tripletnum = [0; cumsum(tripletnum)];

tripletlist = Int64.(tripletlist)
tripletnum  = Int64.(tripletnum)

return tripletlist, tripletnum

end
    
    