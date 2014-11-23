incx = abs(incx)
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = $op( y[dy+i], z[dz+i] )
end
