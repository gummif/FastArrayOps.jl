incx = abs(incx)
@inbounds for i = ix:incx:ix-1+n*incx
    x[i] = $op( a, x[i] )
end
