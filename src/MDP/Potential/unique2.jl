  function unique2(vec::AbstractVector)
      n = length(vec)
      id = unique(vec)
      nUnique = length(id)
      ia = zeros(Int64, nUnique)
      ic = zeros(Int64, n)
      for i = 1:nUnique
          @inbounds ia[i] = findfirst(isequal(id[i]), vec)
      end
      for i = 1:n
          @inbounds ic[i] = findfirst(isequal(vec[i]), id)
      end
      return id, ia, ic
  end
  
  