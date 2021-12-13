function atomlist(xi, pimages, wr, B2R)

dim = size(xi,1);
n = size(xi,2);
m = size(pimages,2);
alist = Int64.(Array(1:n));
k = 0;

if dim==2
    for i = 1:n
        for j = 2:m
            xj = xi[:,i] + pimages[:,j];
            xs = B2R*xj;        
            if (wr[1,1] <= xs[1]) && (xs[1] <= wr[1,3]) &&  (wr[2,1] <= xs[2]) && (xs[2] <= wr[2,3])
                k = k + 1;
                xi = [xi xj];
                append!(alist,i);
            end
        end    
    end
else
    for i = 1:n
        for j = 2:m
            xj = xi[:,i] + pimages[:,j];
            xs = B2R*xj;        
            if (wr[1,1] <= xs[1]) && (xs[1] <= wr[1,7]) &&  (wr[2,1] <= xs[2]) && (xs[2] <= wr[2,7]) &&  (wr[3,1] <= xs[3]) && (xs[3] <= wr[3,7])
                k = k + 1;
                xi = [xi xj];
                append!(alist,i);                
            end
        end    
    end    
end

return xi, alist, n, k

end

