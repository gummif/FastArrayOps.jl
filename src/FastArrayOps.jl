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
    # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    return 0
end
function fast_check2(x, ix, incx, y, iy, incy, n)
    # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    return 0
end
function fast_check3(x, ix, incx, y, iy, incy, z, iz, incz, n)
    # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
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
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    if n < $NLIM_SCALE*incx
        # src: scalarr1_for.jl
incx = abs(incx)
@inbounds for i = ix:incx:ix-1+n*incx
    x[i] = *( a, x[i] )
end
    else
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set1_inc1.jl
incx = 1
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    if n < $NLIM_SCALE #*incx
        # src: scalarr1_for_inc1.jl
@inbounds for i = ix:ix-1+n
    x[i] = *( a, x[i] )
end
    else
        # src: set1_inc1.jl
incx = 1
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end

# x = a*y
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        # src: scalarr1_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = *( a, y[iy+i*incy] )
end 
    else
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    mul = abs(incx)
    if n < $NLIM_SCALE_OOP1*mul || n*mul > $NLIM_SCALE_OOP2
        # src: scalarr1_foroop_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = *( a, y[d+i] )
end
    else
        # src: set2_inceq.jl
incy = incx
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        # src: scalarr1_foroop_inc1.jl
d = iy - ix
@inbounds for i = ix:ix-1+n
    x[i] = *( a, y[d+i] )
end
    else
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_SCALE_OOP1 || n > $NLIM_SCALE_OOP2
        # src: scalarr1_foroop_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = *( a, y[i] )
end
    else
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        # src: blas_scale.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
ccall(($(string(fscal)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), px, &(incx))
    end
    return x
end

# x = x.*y
# ========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || begin
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: arr2xy_for.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = *( x[ix+i*incx], y[iy+i*incy] )
end 
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: arr2xy_for_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = *( x[i], y[d+i] )
end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: arr2xy_for_inc1.jl
d = iy - ix
@inbounds for i = ix:ix-1+n
    x[i] = *( x[i], y[d+i] )
end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: arr2xy_for_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = *( x[i], y[i] )
end
    return x
end

# x = y.*z
# ========

# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, n::Int)
    $isunsafe || begin
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    # src: arr2yz_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
incz < 0 && (iz = iz+(n-1)*abs(incz))
@inbounds for i = 0:n-1
    x[ix+i*incx] = *( y[iy+i*incy], z[iz+i*incz] )
end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        # src: set3_inceq.jl
incy = incx
incz = incx
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    # src: arr2yz_foroop_inceq.jl
incx = abs(incx)
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = *( y[dy+i], z[dz+i] )
end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        # src: set3_inc1.jl
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    # src: arr2yz_foroop_inc1.jl
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:ix-1+n
    x[i] = *( y[dy+i], z[dz+i] )
end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, n::Int)
    $isunsafe || begin
        # src: set3_inc1ieq.jl
iy = ix
iz = ix
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    # src: arr2yz_foroop_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = *( y[i], z[i] )
end 
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
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    # src: scalarr1_for.jl
incx = abs(incx)
@inbounds for i = ix:incx:ix-1+n*incx
    x[i] = +( a, x[i] )
end
    return x
end

# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set1_inc1.jl
incx = 1
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    # src: scalarr1_for_inc1.jl
@inbounds for i = ix:ix-1+n
    x[i] = +( a, x[i] )
end
    return x
end

# x = y + a
# =======
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: scalarr1_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = +( a, y[iy+i*incy] )
end 
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: scalarr1_foroop_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = +( a, y[d+i] )
end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: scalarr1_foroop_inc1.jl
d = iy - ix
@inbounds for i = ix:ix-1+n
    x[i] = +( a, y[d+i] )
end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: scalarr1_foroop_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = +( a, y[i] )
end
    return x
end

# x = x + y
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, n::Int)
    $isunsafe || begin
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR #*mul # || n*mul > $NLIM_SCALEARR
        # src: arr2xy_for.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = +( x[ix+i*incx], y[iy+i*incy] )
end 
    else
        a = 1
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR
        # src: arr2xy_for_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = +( x[i], y[d+i] )
end
    else
        a = 1
        # src: set2_inceq.jl
incy = incx
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR
        # src: arr2xy_for_inc1.jl
d = iy - ix
@inbounds for i = ix:ix-1+n
    x[i] = +( x[i], y[d+i] )
end
    else
        a = 1
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR
        # src: arr2xy_for_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = +( x[i], y[i] )
end
    else
        a = 1
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end

# x = y + z
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, n::Int)
    $isunsafe || begin
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2 #*mul # || n*mul > $NLIM_SCALEARR
        # src: arr2yz_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
incz < 0 && (iz = iz+(n-1)*abs(incz))
@inbounds for i = 0:n-1
    x[ix+i*incx] = +( y[iy+i*incy], z[iz+i*incz] )
end
    else
        a = 1
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        # src: set3_inceq.jl
incy = incx
incz = incx
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        # src: arr2yz_foroop_inceq.jl
incx = abs(incx)
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = +( y[dy+i], z[dz+i] )
end
    else
        a = 1
        # src: set3_inceq.jl
incy = incx
incz = incx
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, n::Int)
    $isunsafe || begin
        # src: set3_inc1.jl
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        # src: arr2yz_foroop_inc1.jl
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:ix-1+n
    x[i] = +( y[dy+i], z[dz+i] )
end
    else
        a = 1
        # src: set3_inc1.jl
incx = 1
incy = 1
incz = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, n::Int)
    $isunsafe || begin
        # src: set3_inc1ieq.jl
iy = ix
iz = ix
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARR_OOP1 || n > $NLIM_ADDARR_OOP2
        # src: arr2yz_foroop_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = +( y[i], z[i] )
end 
    else
        a = 1
        # src: set3_inc1ieq.jl
iy = ix
iz = ix
incx = 1
incy = 1
incz = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
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
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL #*mul # || n*mul > $NLIM_SCALEARR
        # src: addarrscal_for.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = x[ix+i*incx] + y[iy+i*incy]*a
end
    else
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL
        # src: addarrscal_for_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = x[i] + y[d+i]*a
end
    else
        # src: set2_inceq.jl
incy = incx
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL
        # src: addarrscal_for_inc1.jl
d = iy - ix
@inbounds for i = ix:ix-1+n
    x[i] = x[i] + y[d+i]*a
end
    else
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL
        # src: addarrscal_for_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = x[i] + y[i]*a
end
    else
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end

# x = y + a*z
# =========
# general
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, incy::Int, z::Array{$elty}, iz::Int, incz::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2 #*mul # || n*mul > $NLIM_SCALEARR
        # src: addarrscal_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
incz < 0 && (iz = iz+(n-1)*abs(incz))
@inbounds for i = 0:n-1
    x[ix+i*incx] = y[iy+i*incy] + z[iz+i*incz]*a
end
    else
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set3_inceq.jl
incy = incx
incz = incx
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        # src: addarrscal_foroop_inceq.jl
incx = abs(incx)
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = y[dy+i] + z[dz+i]*a
end
    else
        # src: set3_inceq.jl
incy = incx
incz = incx
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, z::Array{$elty}, iz::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set3_inc1.jl
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        # src: addarrscal_foroop_inc1.jl
dy = iy - ix
dz = iz - ix
@inbounds for i = ix:ix-1+n
    x[i] = y[dy+i] + z[dz+i]*a
end
    else
        # src: set3_inc1.jl
incx = 1
incy = 1
incz = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
    end
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, z::Array{$elty}, a::$elty, n::Int)
    $isunsafe || begin
        # src: set3_inc1ieq.jl
iy = ix
iz = ix
incx = 1
incy = 1
incz = 1
        # src: fast_check3.jl
(0 != incx && 0 != incy && 0 != incz) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy && 0 < iz) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
iz+(n-1)*abs(incz) <= length(z) || throw(BoundsError())
    end
    if n < $NLIM_ADDARRSCAL_OOP1 || n > $NLIM_ADDARRSCAL_OOP2
        # src: addarrscal_foroop_inc1ieq.jl
@inbounds for i = ix:ix-1+n
    x[i] = y[i] + z[i]*a
end
    else
        # src: set3_inc1ieq.jl
iy = ix
iz = ix
incx = 1
incy = 1
incz = 1
        # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
        y, iy, incy = z, iz, incz
        # src: blas_axpy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(faxpy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), &(a), py, &(incy), px, &(incx))
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
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    mul = max(abs(incx), abs(incy))
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        # src: copy_foroop.jl
incx < 0 && (ix = ix+(n-1)*abs(incx))
incy < 0 && (iy = iy+(n-1)*abs(incy))
@inbounds for i = 0:n-1
    x[ix+i*incx] = y[iy+i*incy]
end
    else
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
    end
    return x
end
# inceq
function ($f)(x::Array{$elty}, ix::Int, incx::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inceq.jl
incy = incx
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    mul = abs(incx)
    if n < $NLIM_COPY1*mul || n*mul > $NLIM_COPY2
        # src: copy_foroop_inceq.jl
incx = abs(incx)
d = iy - ix
@inbounds for i = ix:incx:ix+(n-1)*incx
    x[i] = y[d+i]
end
    else
        # src: set2_inceq.jl
incy = incx
        # src: blas_copy.jl
px = convert(Ptr{$(elty)},x) + (ix-1)*sizeof($(elty))
py = convert(Ptr{$(elty)},y) + (iy-1)*sizeof($(elty))
ccall(($(string(fcopy)),libblas), Void,
    (Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}, Ptr{$(elty)}, Ptr{BlasInt}),
    &(n), py, &(incy), px, &(incx))
    end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, iy::Int, n::Int)
    $isunsafe || begin
        # src: set2_inc1.jl
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: set2_inc1.jl
incx = 1
incy = 1
    # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
    return x
end
# inc1ieq
function ($f)(x::Array{$elty}, ix::Int, y::Array{$elty}, n::Int)
    $isunsafe || begin
        # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
        # src: fast_check2.jl
(0 != incx && 0 != incy) || throw(ArgumentError("zero increment"))
(0 < ix && 0 < iy) || throw(BoundsError())
ix+(n-1)*abs(incx) <= length(x) || throw(BoundsError())
iy+(n-1)*abs(incy) <= length(y) || throw(BoundsError())
    end
    # src: set2_inc1ieq.jl
iy = ix
incx = 1
incy = 1
    # src: c_memcpy.jl
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
py = convert(Ptr{$(elty)},y) + (iy-1)*selty
ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
    px, py, n*selty)
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
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    # src: fill_foroop.jl
incx = abs(incx)
@inbounds for i = ix:incx:ix-1+n*incx
    x[i] = a
end
    return x
end
# inc1
function ($f)(x::Array{$elty}, ix::Int, a::$elty, n::Int)
    $isunsafe || begin
        # src: set1_inc1.jl
incx = 1
        # src: fast_check1.jl
0 < incx || throw(ArgumentError("non-positive increment"))
0 < ix || throw(BoundsError())
ix+(n-1)*incx <= length(x) || throw(BoundsError())
    end
    if a == 0 && n > $NLIM_FILL && ZEROFLOAT
        # src: set1_inc1.jl
incx = 1
        # src: c_memset.jl
a::Int32 = a
selty = sizeof($(elty))
px = convert(Ptr{$(elty)},x) + (ix-1)*selty
ccall(:memset, Ptr{Void}, (Ptr{Void}, Int32, Csize_t),
    px, a, n*selty)
    else
        # src: fill_foroop_inc1.jl
@inbounds for i = ix:ix-1+n
    x[i] = a
end
    end
    return x
end

end # eval begin
end # for

end # for


end # module

