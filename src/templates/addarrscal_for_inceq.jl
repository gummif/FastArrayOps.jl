incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = x[i] + y[d+i]*a
end
