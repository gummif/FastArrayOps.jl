module FastArrayOps

import Base.LinAlg: BlasReal, BlasComplex, BlasFloat, BlasInt, BlasChar
export fast_scale!, unsafe_fast_scale!

## scale by scalar
# fast_scale!(x, ix, a, n, incx)
# fast_scale!(x, ix, y, iy, a, n, incx, incy)
## scale by array
# fast_scale!(x, ix, y, iy, n, incx, incy)
# fast_scale!(x, ix, y, iy, z, iz, n, incx, incy, incz)
## add scalar
# fast_add!(x, ix, a, n, incx)
# fast_add!(x, ix, y, iy, a, n, incx, incy)
## add array
# fast_add!(x, ix, y, iy, n, incx, incy)
# fast_add!(x, ix, y, iy, z, iz, n, incx, incy, incz)
## add array times constant
# fast_add!(x, ix, a, y, iy, n, incx, incy)
# fast_add!(x, ix, y, iy, a, z, iz, n, incx, incy, incz)
## copy
# fast_copy!(x, ix, y, iy, n, incx, incy)

const NLIM_SCALE = 13
for (f, isunsafe) in ( (:fast_scale!, false), (:unsafe_fast_scale!, true) )
@eval begin

function ($f){T<:BlasFloat}(x::Array{T}, ix::Integer, a::T, n::Integer, incx::Integer=1)
	if !($isunsafe)
		0 < incx || throw(ArgumentError("non-positive increment"))
		0 < ix || throw(BoundsError())
		ix-1+n*incx <= length(x) || throw(BoundsError())
	end
	if n < $NLIM_SCALE && incx == 1
        @inbounds for i = ix:ix-1+n
            x[i] *= a
        end
    else
    	BLAS.scal!(n, a, pointer(x, ix), incx)
    end
    return x
end

function ($f){T<:BlasFloat}(x::Array{T}, ix::Integer, y::Array{T}, iy::Integer, a::T, n::Integer, incx::Integer=1, incy::Integer=1)
	if !($isunsafe)
		(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
		(0 < ix && 0 < iy) || throw(BoundsError())
		ix-1+n*abs(incx) <= length(x) || throw(BoundsError())
		iy-1+n*abs(incy) <= length(y) || throw(BoundsError())
	end
	# or scale by 0.0 then axpy #BLAS.axpy!(n, one(T), pointer(x, ix), incx)
	BLAS.blascopy!(n, pointer(y, iy), incy, pointer(x, ix), incx)
	BLAS.scal!(n, a, pointer(x, ix), abs(incx))
	return x
end

function ($f){T<:BlasFloat}(x::Array{T}, ix::Integer, y::Array{T}, iy::Integer, n::Integer, incx::Integer=1, incy::Integer=1)
	if !($isunsafe)
		(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
		(0 < ix && 0 < iy) || throw(BoundsError())
		ix-1+n*abs(incx) <= length(x) || throw(BoundsError())
		iy-1+n*abs(incy) <= length(y) || throw(BoundsError())
	end
	# TBMV,  GBMV, SBMV
	# scale 0
	# alpha=0, beta=0
	#BLAS.sbmv!('U', k::Integer,
    #                  alpha::($elty), A::StridedMatrix{$elty}, x::StridedVector{$elty}, 
    #                  beta::($elty), y::StridedVector{$elty})
	return x
end

end # eval begin
end # for


end # module

