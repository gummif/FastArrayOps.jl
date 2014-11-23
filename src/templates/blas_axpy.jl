px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
