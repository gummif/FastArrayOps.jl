module FastArrayOps
import Base.LinAlg: BlasReal, BlasComplex, BlasFloat, BlasInt, BlasChar
const libblas = Base.libblas_name

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
## add array times scalar
# fast_add!(x, ix, a, y, iy, n, incx, incy)
# fast_add!(x, ix, y, iy, a, z, iz, n, incx, incy, incz)
## copy array
# fast_copy!(x, ix, y, iy, n, incx, incy)


## Generic for loop functions

function fast_gen_scal{T<:BlasFloat}(x::Array{T}, a::T, n::Int, ix::Int, incx::Int)
    @inbounds for i = ix:incx:n
        x[i] *= a
    end
    return 0
end
function fast_gen_scal_oop{T<:BlasFloat}(x::Array{T}, y::Array{T}, a::T, nel::Int, ix::Int, incx::Int, iy::Int, incy::Int)
    incx < 0 && (ix = ix+(nel-1)*abs(incx))
    incy < 0 && (iy = iy+(nel-1)*abs(incy))
    @inbounds for i = 0:nel-1
        x[ix+i*incx] = y[iy+i*incy]*a
    end
    return 0
end
function fast_inc1_scal_oop{T<:BlasFloat}(x::Array{T}, y::Array{T}, a::T, nel::Int, ix::Int, iy::Int)
    if iy == ix
        @inbounds for i = ix:nel+ix-1
            x[i] = y[i]*a
        end
    else
        d::Int = iy-ix
        @inbounds for i = ix:nel+ix-1
            x[i] = y[d+i]*a
        end
    end
    return 0
end

## cutoff constants

const NLIM_SCALE = 13
const NLIM_SCALE_OOP1 = 80
const NLIM_SCALE_OOP2 = 100000

## FAO scale methods

for (f, isunsafe) in ( (:fast_scale!, false), (:unsafe_fast_scale!, true) )
for (fscal, fcopy, elty) in ((:dscal_,:dcopy_,:Float64), 
                             (:sscal_,:scopy_,:Float32),
                             (:zscal_,:zcopy_,:Complex128), 
                             (:cscal_,:ccopy_,:Complex64))
@eval begin

function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int, incx::Int=1)
    if !($isunsafe)
        0 < incx || throw(ArgumentError("non-positive increment"))
        0 < ix || throw(BoundsError())
        ix-1+n*incx <= length(x) || throw(BoundsError())
    end
    if n < $NLIM_SCALE*incx #|| n < (incx-1)*2000  # && incx == 1  #
        fast_gen_scal(x, a, ix-1+n*incx, ix, incx)  #ix:ix-1+n
        #fast_gen_scal(x, a, ix-1+n, ix, incx)
    else
        #BLAS.scal!(n, a, pointer(x, ix), incx)
        px = pointer(x, ix)
        ccall(($(string(fscal)),libblas), Void,
                  (Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}),
                  &n, &a, px, &incx)
    end
    return x
end

function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int, incx::Int=1, incy::Int=1)
    if !($isunsafe)
        (0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
        (0 < ix && 0 < iy) || throw(BoundsError())
        ix-1+n*abs(incx) <= length(x) || throw(BoundsError())
        iy-1+n*abs(incy) <= length(y) || throw(BoundsError())
    end
    mul = max(abs(incx),abs(incy))
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2 #*mul
        if incx*incy == 1
            fast_inc1_scal_oop(x, y, a, n, ix, iy)
        else
            fast_gen_scal_oop(x, y, a, n, ix, incx, iy, incy)
        end
    else
        #BLAS.blascopy!(n, pointer(y, iy), incy, pointer(x, ix), incx)
        #BLAS.scal!(n, a, pointer(x, ix), abs(incx))
        px = pointer(x, ix)
        py = pointer(y, iy)
        absincx = abs(incx)
        ccall(($(string(fcopy)),libblas), Void,
                (Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}),
                 &n, py, &incy, px, &incx)
        ccall(($(string(fscal)),libblas), Void,
                  (Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}),
                  &n, &a, px, &absincx)
    end
    return x
end

function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int, incx::Int=1, incy::Int=1)
    if !($isunsafe)
        (0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
        (0 < ix && 0 < iy) || throw(BoundsError())
        ix-1+n*abs(incx) <= length(x) || throw(BoundsError())
        iy-1+n*abs(incy) <= length(y) || throw(BoundsError())
    end
    # TBMV,  GBMV, SBMV
    # scale 0
    # alpha=0, beta=0
    #BLAS.sbmv!('U', k::Int,
    #                  alpha::($elty), A::StridedMatrix{$elty}, x::StridedVector{$elty}, 
    #                  beta::($elty), y::StridedVector{$elty})
    return x
end

end # eval begin
end # for
end # for


end # module

