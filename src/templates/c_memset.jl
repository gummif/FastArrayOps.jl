a::Int32 = a
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
ccall(:memset, Ptr{Void}, (Ptr{Void}, Int32, Csize_t),
    px, a, n*selty)
