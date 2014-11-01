using FastArrayOps
using Base.Test

include("testingmacros.jl")


for nmax in 20:30,
    i in 1:20,
    inc in [-10:-1:-1,1:10]
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)) <= nmax
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)+1) > nmax
    @test nmax2nel(i, inc, nel2nmax(i, inc, nmax)) == nmax
    @test fast_range2args(fast_args2range(i, inc, nmax)) == (i, inc, nmax)
    @test fast_args2range(fast_range2args(i:inc:nmax)...) == i:inc:nmax
end # for


const TYPES = (Float64, Float32, Complex128, Complex64)
const IX = (1, 2)
const INC1 = (1, 2)
const INC2 = ((2, 2), (2, 2), (2, 3), (-1, 2), (-1, -3))
const ANUM = 1.123

println("fast_scale!(x, ix, incx, a, n) ...")
for dtype in TYPES,
    ix in IX,
    incx in INC1,
    nmax in (FastArrayOps.NLIM_SCALE-1, FastArrayOps.NLIM_SCALE+1)
    
    n = nmax2nel(ix, incx, nmax)
    @printtest1(dtype, nmax, ix, incx, n)

    x_exp, x0 = makesignals1(dtype, nmax)
    a = convert(dtype,ANUM)
    x_exp[fast_args2range(ix, incx, n)] *= a
    
    @testfunc_arr1a(fast_scale!, x_exp, x0, ix, incx, a, n)
end # for


println("fast_scale!(x, ix, incx, y, iy, incy, a, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (FastArrayOps.NLIM_SCALE_OOP1-1, FastArrayOps.NLIM_SCALE_OOP1+1)
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    a = convert(dtype,ANUM)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] = a*y0[ry]
    
    @testfunc_arr2a(fast_scale!, x_exp, x0, ix, incx, y0, iy, incy, a, n)
end # for


println("fast_copy!(x, ix, incx, y, iy, incy, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (FastArrayOps.NLIM_COPY1-1, FastArrayOps.NLIM_COPY1+1)
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] = y0[ry]
    
    @testfunc_arr2(fast_copy!, x_exp, x0, ix, incx, y0, iy, incy, n)
end # for





