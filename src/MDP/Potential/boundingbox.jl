#****************************************************************************                               
#                     Molecular Dynamics Potentials [MDP]
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen [cuongng@mit.edu, exapde@gmail.com]
#****************************************************************************


function boundingbox2D(pbc::Array{Int64,1}, a::Array{Float64,1}, b::Array{Float64,1}, r::Array{Float64,1})
    # a = [a1, a2] and b = [b1,b2] are two vectors defining a parallelogram 
    # r is the distance from the original parallelogram to the bounding parallelogram 
    # w are four vertices defining the bounding parallelogram 
    
    if length(r)==1
        r = [r; r];
    end
    
    a = a[:]; b = b[:];
    
    a1 = a[1]; a2 = a[2];
    b1 = b[1]; b2 = b[2];
    norma = norm(a);
    normb = norm(b);
    
    # vertices of the parallelogram defined by a and b
    v1 = [0.0; 0.0];
    v2 = a;
    v3 = a+b;
    v4 = b;
    
    # a2*x - a1*y = 0 -> y = a2/a1 * x [1st principal axis]
    # a2*x - a1*y = da a parallel line to 1st principal axis
    # The distance between the two lines is da/norma. As da/norma = r, we get
    da = r[2]*norma;
    # a2*x - a1*y = da is a line which is below and parallel to the 1st axis at a distance r
    
    # apply the same formula to the second axis b
    db = r[1]*normb;
    # b2*x - b1*y = -db is a line which is left and parallel to the 2nd axis at a distance r
    
    # intersection of a2*x - a1*y = da and b2*x - b1*y = -db
    w1 = [a2 -a1; b2 -b1]\[da; -db]; # the 1st vertex of the bounding parallelogram
    
    # find e = [e1,e2] such that e1*a + e2*b = w1
    e = [a b]\w1;
    
    # distance between w1 and e[1]*a
    sb = norm(w1-e[1]*a);
    # distance between w1 and e[2]*b
    sa = norm(w1-e[2]*b);
    
    # # length of the bounding parallelogram along the 1st axis
    # l1 = norm[a] + 2*sa;
    # # length of the bounding parallelogram along the 2nd axis
    # l2 = norm[b] + 2*sb;
    
    # length of the bounding parallelogram along the 1st axis
    l1 = norm(a) + 2*sa*pbc[1];
    # length of the bounding parallelogram along the 2nd axis
    l2 = norm(b) + 2*sb*pbc[2];
    
    # the 1st vertex of  the bounding parallelepiped
    w1 = pbc[1]*e[1]*a + pbc[2]*e[2]*b;
    # the 2nd vertex of  the bounding parallelogram
    w2 = w1 + l1*a/norm(a);
    # the 3rd vertex of  the bounding parallelogram
    w3 = v3 - w1;
    # the 4th vertex of  the bounding parallelogram
    w4 = w1 + l2*b/norm(b);
    
    v = [v1 v2 v3 v4];
    w = [w1 w2 w3 w4];

    return v, w
    
end


function boundingbox3D(pbc::Array{Int64,1}, a::Array{Float64,1}, b::Array{Float64,1}, c::Array{Float64,1}, r::Array{Float64,1})
    # a = (a1, a2, a3), b = (b1,b2,b3), c=(c1,c2,c3) are three vectors defining a parallelepiped 
    # r is the distance from the original parallelepiped to the bounding parallelepiped 
    # w are 8 vertices defining the bounding parallelepiped 
    
    if length(r)==1
        r = [r; r; r];
    end
    
    a = a[:]; b = b[:]; c = c[:];
    
    a1 = a[1]; a2 = a[2]; a3 = a[3];
    b1 = b[1]; b2 = b[2]; b3 = b[3];
    c1 = c[1]; c2 = c[2]; c3 = c[3];
    norma = norm(a);
    normb = norm(b);
    normc = norm(c);
    
    # vertices of the parallelogram defined by a, b, and c
    v1 = [0; 0; 0];
    v2 = a;
    v3 = a+b;
    v4 = b;
    v5 = c;
    v6 = a+c;
    v7 = a+b+c;
    v8 = b+c;
    
    # the 1st plane defined by a and b
    p1 = [a2*b3-a3*b2; a3*b1-a1*b3; a1*b2-a2*b1];
    normp1 = sqrt(p1[1]*p1[1] + p1[2]*p1[2] + p1[3]*p1[3]);
    
    # Since the distance between (0,0,0) and the 1st plane p1[1]*x + p1[2]*y + p1[3]*z = d1
    # is equal to r[3], we have
    d1 = -sign(p1[3])*r[3]*normp1;
    
    # the 2nd plane defined by b and c
    p2 = [b2*c3-b3*c2; b3*c1-b1*c3; b1*c2-b2*c1];
    normp2 = sqrt(p2[1]*p2[1] + p2[2]*p2[2] + p2[3]*p2[3]);
    
    # Since the distance between (0,0,0) and the 2nd plane p2[1]*x + p2[2]*y + p2[3]*z = d2
    # is equal to r[1], we have
    d2 = -sign(p2[1])*r[1]*normp2;
    
    # the 3rd plane defined by c and a
    p3 = [c2*a3-c3*a2; c3*a1-c1*a3; c1*a2-c2*a1];
    normp3 = sqrt(p3[1]*p3[1] + p3[2]*p3[2] + p3[3]*p3[3]);
    
    # Since the distance between (0,0,0) and the 2nd plane p3(1)*x + p3(2)*y + p3(3)*z = d3
    # is equal to r[2], we have
    d3 = -sign(p3[2])*r[2]*normp3;
    
    # intersection of the above three planes    
    A = reshape([p1; p2; p3],3,3)';
    w1 = A\[d1; d2; d3];
    
    # find e = (e1,e2,e3) such that e1*a + e2*b + e3*c = w1
    e = [a b c]\w1;
    
    # distance between w1 and the point e2*b + e3*c
    sa = norm(w1 - e[2]*b - e[3]*c);
    
    # distance between w1 and the point e1*a + e3*c
    sb = norm(w1 - e[1]*a - e[3]*c);
    
    # distance between w1 and the point e1*a + e2*b
    sc = norm(w1 - e[1]*a - e[2]*b);
     
    # length of the bounding parallelepiped along the 1st axis
    la = norma + 2*sa*pbc[1];
    # length of the bounding parallelepiped along the 2nd axis
    lb = normb + 2*sb*pbc[2];
    # length of the bounding parallelepiped along the 3rd axis
    lc = normc + 2*sc*pbc[3];
    
    # the 1st vertex of  the bounding parallelepiped
    w1 = pbc[1]*e[1]*a + pbc[2]*e[2]*b + pbc[3]*e[3]*c;    
    # the 2nd vertex of  the bounding parallelepiped
    w2 = w1 + la*a/norma;
    # the 4th vertex of  the bounding parallelepiped
    w4 = w1 + lb*b/normb;
    # the 3rd vertex of  the bounding parallelepiped
    w3 = w2 + w4 - w1;
    # the 5th vertex of  the bounding parallelepiped
    w5 = w1 + lc*c/normc;
    # the 6th vertex of  the bounding parallelepiped
    w6 = w5 + la*a/norma;
    # the 8th vertex of  the bounding parallelepiped
    w8 = w5 + lb*b/normb;
    # the 7th vertex of  the bounding parallelepiped
    w7 = w6 + w8 - w5;
    
    v = [v1 v2 v3 v4 v5 v6 v7 v8];
    w = [w1 w2 w3 w4 w5 w6 w7 w8];    

    return v, w

end


function boundingbox(pbc::Array{Int64,1}, rmax::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1}, c=nothing)
    # pbc periodic boundary conditions
    # a = (a1, a2, a3), b = (b1,b2,b3), c=(c1,c2,c3) are three vectors defining a parallelepiped 
    # rmax is the distance from the original parallelepiped to the bounding parallelepiped 
    # v are 8 vertices defining the original parallelepiped 
    # w are 8 vertices defining the bounding parallelepiped 
    
    if c===nothing
        v, w = boundingbox2D(pbc, a, b, rmax);
    else    
        v, w = boundingbox3D(pbc, a, b, c, rmax);
    end

    return v, w
end

    




    
    