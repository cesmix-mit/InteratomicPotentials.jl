function makejkl(n)
    
m = (n-2)*(n-1)*n/6;
indj = zeros(Int64,m);
indk = zeros(Int64,m);
indl = zeros(Int64,m);
count = 0;
for j = 1:(n-2)
    for k = (j+1):(n-1)
        for l = (k+1):n                
            count = count+1;
            indj[count] = j;
            indk[count] = k;
            indl[count] = l;
        end
    end
end

return indj, indk, indl

end

function neighquadletlist(pairlist, pairnumsum, ilist)

pairnum = pairnumsum[2:end]-pairnumsum[1:end-1];
quadletnum  = (pairnum-2).*(pairnum-1).*pairnum/6;

inum = length(ilist);
quadletlistj = [];
quadletlistk = [];
quadletlistl = [];

for ii = 1:inum
    n1 = pairnumsum[ii];
    n2 = pairnumsum[ii+1];    
    g = pairlist[(n1+1):n2]; # a list of neighbors around atom i
    
    indj, indk, indl = makejkl(n);            
    quadletlistj = [quadletlistj; g[indj]];      
    quadletlistk = [quadletlistk; g[indk]];      
    quadletlistl = [quadletlistl; g[indl]];      
end

quadletlist = [quadletlistj[:] quadletlistk[:] quadletlistl[:]];
quadletnum = [0; cumsum(quadletnum)];

quadletlist = Int64.(quadletlist)
quadletnum  = Int64.(quadletnum)

return quadletlist, quadletnum

end



