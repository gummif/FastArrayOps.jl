px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
