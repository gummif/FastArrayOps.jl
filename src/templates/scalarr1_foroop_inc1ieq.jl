@inbounds for i = ix:ix-1+n
    x[i] = $op( a, y[i] )
end
