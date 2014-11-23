(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
