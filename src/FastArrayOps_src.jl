module FastArrayOps
import Base.LinAlg: BlasReal, BlasComplex, BlasFloat, BlasInt, BlasChar
const libblas = Base.libblas_name
export fast_scale!,     unsafe_fast_scale!, 
       fast_add!,       unsafe_fast_add!,
       fast_addscal!,   unsafe_fast_addscal!,
       fast_copy!,      unsafe_fast_copy!,
       fast_fill!,      unsafe_fast_fill!
export fast_check1, fast_check2, fast_check3, nmax2nel, nel2nmax, fast_args2range, fast_range2args

# WARNING: FastArrayOps.jl gets overwritten by FastArrayOps_src.jl when running make.jl

## CONSTANTS

const NLIM_SCALE = 13
const NLIM_SCALE_OOP1 = 80
const NLIM_SCALE_OOP2 = 100000
const NLIM_SCALEARR = typemax(Int)
const NLIM_SCALEARR_OOP = typemax(Int)
const NLIM_ADD = typemax(Int)
const NLIM_ADD_OOP = typemax(Int)
const NLIM_ADDARR = 13
const NLIM_ADDARR_OOP1 = 30
const NLIM_ADDARR_OOP2 = 100000
const NLIM_ADDARRSCAL = 13
const NLIM_ADDARRSCAL_OOP1 = 30
const NLIM_ADDARRSCAL_OOP2 = 100000
const NLIM_COPY1 = 80
const NLIM_COPY2 = 100000
const NLIM_FILL = 13


## UTILS

function nmax2nel(i::Int, inc::Int, nmax::Int)
    @assert 0 < i
    nmax < i && return 0
    return div(nmax - i, abs(inc)) + 1
end
function nel2nmax(i::Int, inc::Int, nel::Int)
    @assert 0 < i
    nel < 0 && return i - 1
    return i + (nel- 1)*abs(inc)
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

function fast_check1(x, ix, incx, n)
    ##!fast_check1
    return 0
end
function fast_check2(x, ix, incx, y, iy, incy, n)
    ##!fast_check2
    return 0
end
function fast_check3(x, ix, incx, y, iy, incy, z, iz, incz, n)
    ##!fast_check3
    return 0
end

# utils for 0 fill
inttype(::Type{Float32}) = Int32
inttype(::Type{Float64}) = Int64
inttype(::Type{Complex64}) = Int64
inttype(::Type{Complex128}) = Int128
function fast_reinterpret1{T<:Number}(::Type{T}, a::Array)
    @assert length(a) == 1
    ccall(:jl_reshape_array, Array{T,1}, (Any, Any, Any), Array{T,1}, a, (1,))
end
# are float 0s zero bits?
function zerobits()
    v = reinterpret(inttype(Float32), convert(Float32, 0))
    v += reinterpret(inttype(Float64), convert(Float64, 0))
    v += fast_reinterpret1(inttype(Complex64), [convert(Complex64, 0)])[1]
    v += fast_reinterpret1(inttype(Complex128), [convert(Complex128, 0)])[1]
    if v == 0
        return true
    else
        return false
    end
end
const ZEROFLOAT = zerobits()


for (fscal, fcopy, faxpy, ftbmv, fsbmv, elty) in (
                        (:dscal_, :dcopy_, :daxpy_, :dtbmv_, :dsbmv_, :Float64), 
                        (:sscal_, :scopy_, :saxpy_, :stbmv_, :ssbmv_, :Float32),
                        (:zscal_, :zcopy_, :zaxpy_, :ztbmv_, :zsbmv_, :Complex128), 
                        (:cscal_, :ccopy_, :caxpy_, :ctbmv_, :csbmv_, :Complex64))

## SCALE METHODS

for (f, isunsafe) in ( (:fast_scale!, false), (:unsafe_fast_scale!, true) )
@eval begin

# x = a*x
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check1
    end
    if n < $NLIM_SCALE*incx
        ##!scalarr1_for#*
    else
        ##!blas_scale
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set1_inc1
        ##!fast_check1
    end
    if n < $NLIM_SCALE #*incx
        ##!scalarr1_for_inc1#*
    else
        ##!set1_inc1
        ##!blas_scale
    end
    return x
end

# x = a*y
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        ##!scalarr1_foroop#*
    else
        ##!blas_copy
        ##!blas_scale
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    mul = abs(incx)
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        ##!scalarr1_foroop_inceq#*
    else
        ##!set2_inceq
        ##!blas_copy
        ##!blas_scale
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        ##!scalarr1_foroop_inc1#*
    else
        ##!set2_inc1
        ##!c_memcpy
        ##!blas_scale
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        ##!scalarr1_foroop_inc1ieq#*
    else
        ##!set2_inc1ieq
        ##!c_memcpy
        ##!blas_scale
    end
    return x
end

# x = x.*y
# ========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    ##!arr2xy_for#*
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    ##!arr2xy_for_inceq#*
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    ##!arr2xy_for_inc1#*
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    ##!arr2xy_for_inc1ieq#*
    return x
end

# x = y.*z
# ========

# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, n::Int)
    $isunsafe || begin
        ##!fast_check3
    end
    ##!arr2yz_foroop#*
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        ##!set3_inceq
        ##!fast_check3
    end
    ##!arr2yz_foroop_inceq#*
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        ##!set3_inc1
        ##!fast_check3
    end
    ##!arr2yz_foroop_inc1#*
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, n::Int)
    $isunsafe || begin
        ##!set3_inc1ieq
        ##!fast_check3
    end
    ##!arr2yz_foroop_inc1ieq#*
    return x
end


end # eval begin
end # for


## ADD METHODS

for (f, isunsafe) in ( (:fast_add!, false), (:unsafe_fast_add!, true) )
@eval begin

# x = x + a
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check1
    end
    ##!scalarr1_for#+
    return x
end

# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set1_inc1
        ##!fast_check1
    end
    ##!scalarr1_for_inc1#+
    return x
end

# x = y + a
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    ##!scalarr1_foroop#+
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    ##!scalarr1_foroop_inceq#+
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    ##!scalarr1_foroop_inc1#+
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    ##!scalarr1_foroop_inc1ieq#+
    return x
end

# x = x + y
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    if n < $NLIM_ADDARR #*mul # || n*mul > $NLIM_SCALEARR
        ##!arr2xy_for#+
    else
        a = 1
        ##!blas_axpy
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    if n < $NLIM_ADDARR
        ##!arr2xy_for_inceq#+
    else
        a = 1
        ##!set2_inceq
        ##!blas_axpy
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    if n < $NLIM_ADDARR
        ##!arr2xy_for_inc1#+
    else
        a = 1
        ##!set2_inc1
        ##!blas_axpy
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    if n < $NLIM_ADDARR
        ##!arr2xy_for_inc1ieq#+
    else
        a = 1
        ##!set2_inc1ieq
        ##!blas_axpy
    end
    return x
end

# x = y + z
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, n::Int)
    $isunsafe || begin
        ##!fast_check3
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2 #*mul # || n*mul > $NLIM_SCALEARR
        ##!arr2yz_foroop#+
    else
        a = 1
        ##!blas_copy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        ##!set3_inceq
        ##!fast_check3
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        ##!arr2yz_foroop_inceq#+
    else
        a = 1
        ##!set3_inceq
        ##!blas_copy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        ##!set3_inc1
        ##!fast_check3
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        ##!arr2yz_foroop_inc1#+
    else
        a = 1
        ##!set3_inc1
        ##!c_memcpy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, n::Int)
    $isunsafe || begin
        ##!set3_inc1ieq
        ##!fast_check3
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        ##!arr2yz_foroop_inc1ieq#+
    else
        a = 1
        ##!set3_inc1ieq
        ##!c_memcpy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end

end # eval begin
end # for


## ADDSCAL METHODS

for (f, isunsafe) in ( (:fast_addscal!, false), (:unsafe_fast_addscal!, true) )
@eval begin

# x = x + a*y
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    if n < $NLIM_ADDARRSCAL #*mul # || n*mul > $NLIM_SCALEARR
        ##!addarrscal_for
    else
        ##!blas_axpy
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    if n < $NLIM_ADDARRSCAL
        ##!addarrscal_for_inceq
    else
        ##!set2_inceq
        ##!blas_axpy
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    if n < $NLIM_ADDARRSCAL
        ##!addarrscal_for_inc1
    else
        ##!set2_inc1
        ##!blas_axpy
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    if n < $NLIM_ADDARRSCAL
        ##!addarrscal_for_inc1ieq
    else
        ##!set2_inc1ieq
        ##!blas_axpy
    end
    return x
end

# x = y + a*z
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check3
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2 #*mul # || n*mul > $NLIM_SCALEARR
        ##!addarrscal_foroop
    else
        ##!blas_copy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set3_inceq
        ##!fast_check3
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        ##!addarrscal_foroop_inceq
    else
        ##!set3_inceq
        ##!blas_copy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set3_inc1
        ##!fast_check3
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        ##!addarrscal_foroop_inc1
    else
        ##!set3_inc1
        ##!c_memcpy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        ##!set3_inc1ieq
        ##!fast_check3
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        ##!addarrscal_foroop_inc1ieq
    else
        ##!set3_inc1ieq
        ##!c_memcpy
        y, iy, incy = z, iz, incz
        ##!blas_axpy
    end
    return x
end

end # eval begin
end # for


## COPY METHODS

for (f, isunsafe) in ( (:fast_copy!, false), (:unsafe_fast_copy!, true) )
@eval begin

# x = y
# =====
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || begin
        ##!fast_check2
    end
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        ##!copy_foroop
    else
        ##!blas_copy
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inceq
        ##!fast_check2
    end
    mul = abs(incx)
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        ##!copy_foroop_inceq
    else
        ##!set2_inceq
        ##!blas_copy
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        ##!set2_inc1
        ##!fast_check2
    end
    ##!set2_inc1
    ##!c_memcpy
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        ##!set2_inc1ieq
        ##!fast_check2
    end
    ##!set2_inc1ieq
    ##!c_memcpy
    return x
end


end # eval begin
end # for


## FILL METHODS

for (f, isunsafe) in ( (:fast_fill!, false), (:unsafe_fast_fill!, true) )
@eval begin

# x = a
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!fast_check1
    end
    ##!fill_foroop
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        ##!set1_inc1
        ##!fast_check1
    end
    if a == 0 && n > $NLIM_FILL && ZEROFLOAT
        ##!set1_inc1
        ##!c_memset
    else
        ##!fill_foroop_inc1
    end
    return x
end

end # eval begin
end # for

end # for


end # module

