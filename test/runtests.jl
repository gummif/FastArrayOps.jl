using FastArrayOps
using Base.Test

include("testingmacros.jl")


for nmax in 20:30,
    i in 1:20,
    inc in [-10:1:-1,1:10]
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)) <= nmax
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)+1) > nmax
    @test nmax2nel(i, inc, nel2nmax(i, inc, nmax)) == nmax
    @test fast_range2args(fast_args2range(i, inc, nmax)) == (i, inc, nmax)
    @test fast_args2range(fast_range2args(i:inc:nmax)...) == i:inc:nmax
end # for


const TYPES = (Float64, Float32, Complex128, Complex64)
const IX = (1, 2)
const INC1 = (1, 2)
const INC2 = ((1, 1), (2, 2), (2, 3), (-1, 2), (-1, -3))
const INC3 = ((1, 1, 1), (2, 2, 2), (2, 3, 1), (-1, 2, 3), (-1, -3, -2))
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

println("fast_scale!(x, ix, incx, y, iy, incy, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (15,22)  # NLIM is max Int
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] .*= y0[ry]
    
    @testfunc_arr2(fast_scale!, x_exp, x0, ix, incx, y0, iy, incy, n)
end # for

println("fast_scale!(x, ix, incx, y, iy, incy, z, iz, incz, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy, incz) in INC3,
    nmax in (15,22)  # NLIM is max Int
    
    iy = ix
    iz = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    nz = nmax2nel(iz, incz, nmax)
    n = min(min(nx, ny), nz)
    @printtest3(dtype, nmax, ix, incx, iy, incy, iz, incz, n)

    x_exp, x0, y0, z0 = makesignals3(dtype, nmax)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    rz = fast_args2range(iz, incz, n)
    x_exp[rx] = y0[ry].*z0[rz]
    
    @testfunc_arr3(fast_scale!, x_exp, x0, ix, incx, y0, iy, incy, z0, iz, incz, n)
end # for

println("fast_add!(x, ix, incx, a, n) ...")
for dtype in TYPES,
    ix in IX,
    incx in INC1,
    nmax in (15,22)  # NLIM is max Int
    
    n = nmax2nel(ix, incx, nmax)
    @printtest1(dtype, nmax, ix, incx, n)

    x_exp, x0 = makesignals1(dtype, nmax)
    a = convert(dtype,ANUM)
    x_exp[fast_args2range(ix, incx, n)] += a
    
    @testfunc_arr1a(fast_add!, x_exp, x0, ix, incx, a, n)
end # for

println("fast_add!(x, ix, incx, y, iy, incy, a, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (15,22)  # NLIM is max Int
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    a = convert(dtype,ANUM)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] = y0[ry] + a
    
    @testfunc_arr2a(fast_add!, x_exp, x0, ix, incx, y0, iy, incy, a, n)
end # for

println("fast_add!(x, ix, incx, y, iy, incy, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (FastArrayOps.NLIM_ADDARR-1, FastArrayOps.NLIM_ADDARR+1)
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] .+= y0[ry]
    
    @testfunc_arr2(fast_add!, x_exp, x0, ix, incx, y0, iy, incy, n)
end # for

println("fast_add!(x, ix, incx, y, iy, incy, z, iz, incz, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy, incz) in INC3,
    nmax in (FastArrayOps.NLIM_ADDARR_OOP1-1, FastArrayOps.NLIM_ADDARR_OOP1+1)
    
    iy = ix
    iz = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    nz = nmax2nel(iz, incz, nmax)
    n = min(min(nx, ny), nz)
    @printtest3(dtype, nmax, ix, incx, iy, incy, iz, incz, n)

    x_exp, x0, y0, z0 = makesignals3(dtype, nmax)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    rz = fast_args2range(iz, incz, n)
    x_exp[rx] = y0[ry] .+ z0[rz]
    
    @testfunc_arr3(fast_add!, x_exp, x0, ix, incx, y0, iy, incy, z0, iz, incz, n)
end # for

println("fast_addscal!(x, ix, incx, y, iy, incy, a, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy) in INC2,
    nmax in (FastArrayOps.NLIM_ADDARRSCAL-1, FastArrayOps.NLIM_ADDARRSCAL+1)
    
    iy = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    n = min(nx, ny)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)

    x_exp, x0, y0 = makesignals2(dtype, nmax)
    a = convert(dtype,ANUM)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] .+= a*y0[ry]
    
    @testfunc_arr2a(fast_addscal!, x_exp, x0, ix, incx, y0, iy, incy, a, n)
end # for

println("fast_addscal!(x, ix, incx, y, iy, incy, z, iz, incz, a, n) ...")
for dtype in TYPES,
    ix in IX,
    (incx, incy, incz) in INC3,
    nmax in (FastArrayOps.NLIM_ADDARRSCAL_OOP1-1, FastArrayOps.NLIM_ADDARRSCAL_OOP1+1)
    
    iy = ix
    iz = ix
    nx = nmax2nel(ix, incx, nmax)
    ny = nmax2nel(iy, incy, nmax)
    nz = nmax2nel(iz, incz, nmax)
    n = min(min(nx, ny), nz)
    @printtest3(dtype, nmax, ix, incx, iy, incy, iz, incz, n)

    x_exp, x0, y0, z0 = makesignals3(dtype, nmax)
    a = convert(dtype,ANUM)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    rz = fast_args2range(iz, incz, n)
    x_exp[rx] = y0[ry] .+ a*z0[rz]
    
    @testfunc_arr3a(fast_addscal!, x_exp, x0, ix, incx, y0, iy, incy, z0, iz, incz, a, n)
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

println("fast_fill!(x, ix, incx, a, n) ...")
for dtype in TYPES,
    ix in IX,
    incx in INC1,
    nmax in (FastArrayOps.NLIM_FILL-1, FastArrayOps.NLIM_FILL+1),
    aa in (ANUM, 0)
    
    n = nmax2nel(ix, incx, nmax)
    @printtest1(dtype, nmax, ix, incx, n)

    x_exp, x0 = makesignals1(dtype, nmax)
    a = convert(dtype,aa)
    x_exp[fast_args2range(ix, incx, n)] = a
    
    @testfunc_arr1a(fast_fill!, x_exp, x0, ix, incx, a, n)
end # for

println("\nTESTING: SUCCESS")



