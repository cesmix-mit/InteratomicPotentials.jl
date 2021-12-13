function accumarray(subs, val, sz=(maximum(subs),))
      A = zeros(eltype(val), sz...)
      for i = 1:length(val)
          @inbounds A[subs[i]] += val[i]
      end
      A
end
