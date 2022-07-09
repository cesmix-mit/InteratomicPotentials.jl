function neighsingles(x, q, atomtype, ilist)

if !isempty(x)
    xi = x[:,ilist];
else
    xi = [];
end
if !isempty[q]
    qi = q[:,ilist];
else
    qi = [];
end

ai = ilist;
ti = atomtype[ilist];

return xi, qi, ai, ti

end

