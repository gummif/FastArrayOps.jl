incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
incz < 0 && (iz = iz+(n-1)*abs(incz))
@inbounds for i = 0:n-1
    x[ix+i*incx] = $op( y[iy+i*incy], z[iz+i*incz] )
end
