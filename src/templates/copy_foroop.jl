incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = y[iy+i*incy]
end
