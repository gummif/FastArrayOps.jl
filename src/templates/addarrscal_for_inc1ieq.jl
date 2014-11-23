@inbounds for i = ix:ix-1+n
    x[i] = x[i] + y[i]*a
end
