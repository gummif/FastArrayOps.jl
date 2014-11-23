0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
