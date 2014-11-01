module FastArrayOps
import Base.LinAlg: BlasReal, BlasComplex, BlasFloat, BlasInt, BlasChar
const libblas = Base.libblas_name

export fast_scale!, unsafe_fast_scale!, fast_copy!, unsafe_fast_copy!
export @fast_check1, @fast_check2, nmax2nel, nel2nmax, fast_args2range, fast_range2args

## CONSTANTS

const NLIM_SCALE = 13
const NLIM_SCALE_OOP1 = 80
const NLIM_SCALE_OOP2 = 100000
const NLIM_SCALEARR = typemax(Int)
const NLIM_COPY1 = 80
const NLIM_COPY2 = 100000

const FAO_ZERO = zero(Int)
const FAO_ONE = one(Int)

## UTILS

function nmax2nel(i::Int, inc::Int, nmax::Int)
    @assert 0 < i
    nmax < i && return FAO_ZERO
    return div(nmax - i, abs(inc)) + FAO_ONE
end
function nel2nmax(i::Int, inc::Int, nel::Int)
    @assert 0 < i
    nel < 0 && return i - FAO_ONE
    return i + (nel- FAO_ONE)*abs(inc)
end
function fast_args2range(i::Int, inc::Int, n::Int)
    @assert 0 < i
    r = i:abs(inc):nel2nmax(i, abs(inc), n)
    if inc < 0
        r = reverse(r)
    end
    return r
end
function fast_range2args(r::Range)
    inc = step(r)
    if inc > 0
        return (first(r), inc, length(r))
    else
        return (last(r), inc, length(r))
    end
end

## MACROS
include("macros.jl")


## FAO scale methods

for (fscal, fcopy, elty) in ((:dscal_,:dcopy_,:Float64), 
                             (:sscal_,:scopy_,:Float32),
                             (:zscal_,:zcopy_,:Complex128), 
                             (:cscal_,:ccopy_,:Complex64))
for (f, isunsafe) in ( (:fast_scale!, false), (:unsafe_fast_scale!, true) )
@eval begin


# x = a*x
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, a::$elty, n::Int)
    $isunsafe || @fast_check1(x, ix, incx, n)
    if n < $NLIM_SCALE*incx
        @scale_for(x, ix, incx, a, n, $FAO_ONE)
    else
        @scale_blas($(string(fscal)), $(elty), x, ix, incx, a, n)
    end
    return x
end

# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || @fast_check1(x, ix, $FAO_ONE, n)
    if n < $NLIM_SCALE #*incx
        @scale_for_inc1(x, ix, a, n, $FAO_ONE)
    else
        @scale_blas($(string(fscal)), $(elty), x, ix, $FAO_ONE, a, n)
    end
    return x
end

# x = a*y
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incy, n)
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        @scale_foroop(x, ix, incx, y, iy, incy, a, n, $FAO_ZERO, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, incx, y, iy, incy, n)
        @scale_blas($(string(fscal)), $(elty), x, ix, incx, a, n)
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incx, n)
    mul = abs(incx)
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        @scale_foroop_inceq(x, ix, incx, y, iy, a, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, incx, y, iy, incx, n)
        @scale_blas($(string(fscal)), $(elty), x, ix, incx, a, n)
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        @scale_foroop_inc1(x, ix, y, iy, a, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
        @scale_blas($(string(fscal)), $(elty), x, ix, $FAO_ONE, a, n)
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        @scale_foroop_inc1ieq(x, ix, y, a, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
        @scale_blas($(string(fscal)), $(elty), x, ix, $FAO_ONE, a, n)
    end
    return x
end

# x = x.*y
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incy, n)
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_SCALEARR #*mul # || n*mul > $NLIM_SCALEARR
        @scalearr_for(x, ix, incx, y, iy, incy, n, $FAO_ZERO, $FAO_ONE)
    else
        @vecmult_blas($(string(:dtbmv_)), $(elty), x, ix, incx, y, iy, incy, n)
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incx, n)
    mul = abs(incx)
    if n < $NLIM_SCALEARR #*mul # || n*mul > $NLIM_SCALEARR
        @scalearr_for_inceq(x, ix, incx, y, iy, n, $FAO_ONE)
    else
        @vecmult_blas($(string(:dtbmv_)), $(elty), x, ix, incx, y, iy, incx, n)
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
    if n < $NLIM_SCALEARR
        @scalearr_for_inc1(x, ix, y, iy, n, $FAO_ONE)
    else
        @vecmult_blas($(string(:dtbmv_)), $(elty), x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
    if n < $NLIM_SCALEARR
        @scalearr_for_inc1ieq(x, ix, y, n, $FAO_ONE)
    else
        @vecmult_blas($(string(:dtbmv_)), $(elty), x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
    end
    return x
end

# x = y.*z
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, n::Int)
    @fast_check3(x, ix, incx, y, iy, incx, z, iz, incz, n)
    @vecmultoop_blas($(string(:dsbmv_)), $(elty), x, ix, incx, y, iy, incy, z, iz, incz, n)
    return x
end



end # eval begin
end # for


for (f, isunsafe) in ( (:fast_copy!, false), (:unsafe_fast_copy!, true) )
@eval begin

# x = y
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incy, n)
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        @copy_foroop(x, ix, incx, y, iy, incy, n, $FAO_ZERO, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, incx, y, iy, incy, n)
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, incx, y, iy, incx, n)
    mul = abs(incx)
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        @copy_foroop_inceq(x, ix, incx, y, iy, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, incx, y, iy, incx, n)
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
    if n < $NLIM_COPY1 || n > $NLIM_COPY2
        @copy_foroop_inc1(x, ix, y, iy, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, $FAO_ONE, y, iy, $FAO_ONE, n)
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || @fast_check2(x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
    if n < $NLIM_COPY1 || n > $NLIM_COPY2
        @copy_foroop_inc1ieq(x, ix, y, n, $FAO_ONE)
    else
        @copy_blas($(string(fcopy)), $(elty), x, ix, $FAO_ONE, y, ix, $FAO_ONE, n)
    end
    return x
end


end # eval begin
end # for

end # for







end # module

