function accumarray1d(A, subs, val)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
    return A
end

function accumarray2d(A, subs, val)
    for i = 1:length(val)
        A[:,subs[i]] = A[:,subs[i]] .+ val[:,i]
    end
    return A
end

function accumarray3d(A, subs, val)
    for i = 1:length(val)
        A[:,subs[i],:] = A[:,subs[i],:] .+ val[:,i,:]
    end
    return A
end

function accumarray(subs, val)
    n = maximum(subs)
    A = zeros(eltype(val), n)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
    A
end

function accumarray2d(subs, val)
    n = maximum(subs)
    m = size(val,1)
    A = zeros(eltype(val), (m, n))
    for i = 1:length(val)
        A[:,subs[i]] = A[:,subs[i]] .+ val[:,i]
    end
    return A
end

function accumarray3d(subs, val)
    n = maximum(subs)
    m = size(val,1)
    p = size(val,3)
    A = zeros(eltype(val), (m, n, p))
    for i = 1:length(val)
        A[:,subs[i],:] = A[:,subs[i],:] .+ val[:,i,:]
    end
    return A
end
