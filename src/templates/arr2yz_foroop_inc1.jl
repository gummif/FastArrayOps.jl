dy = iy - ix
dz = iz - ix
@inbounds for i = ix:ix-1+n
    x[i] = $op( y[dy+i], z[dz+i] )
end
