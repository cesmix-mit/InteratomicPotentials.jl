function ljpotential(xij, qi, qj, ti, tj, mu, eta, kappa)

A = mu[1];
B = mu[2];

r2 = xij[1,:].*xij[1,:] .+ xij[2,:].*xij[2,:] .+ xij[3,:].*xij[3,:];    
#r1 = sqrt.(r2);
r4 = r2.*r2;
r6 = r2.*r4;
r12 = r6.*r6;
eij = (A./r12 - B./r6);        

r8 = r6.*r2;
r14 = r12.*r2;
fij = xij;
for i = 1:3
    fij[i,:] = xij[i,:].*((6*B)./r8 - (12*A)./r14); 
end
# r7 = r6.*r2;
# r13 = r12.*r1;
# dr = xij./r1;
# fij = dr.*((6*B)./r7 - (12*A)./r13); 

return eij, fij

end


    