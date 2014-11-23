incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = $op( a, y[d+i] )
end
