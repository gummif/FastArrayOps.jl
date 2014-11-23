@inbounds for i = ix:ix-1+n
    x[i] = $op( y[i], z[i] )
end 
