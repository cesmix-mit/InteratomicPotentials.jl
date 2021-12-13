function findatomtype(alist, atomtype, typei)

atype = atomtype[alist]; # atom types of alist
ilist = alist[atype .== typei];    # return atoms of type i

return ilist

end
