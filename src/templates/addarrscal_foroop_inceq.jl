incx = abs(incx)
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = y[dy+i] + z[dz+i]*a
end
