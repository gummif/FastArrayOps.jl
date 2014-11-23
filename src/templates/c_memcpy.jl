selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
