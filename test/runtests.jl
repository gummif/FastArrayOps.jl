using FastArrayOps
using Base.Test

macro printtest1(dtype, nmax, ix, incx, n)
    quote
        nmax=$(esc(nmax))
        ix=$(esc(ix))
        incx=$(esc(incx))
        n=$(esc(n))
        $(esc(dtype))<:Float64 && println("  nmax: $nmax, ix: $ix, incx: $incx, n: $n")
    end
end
macro printtest2(dtype, nmax, ix, incx, iy, incy, n)
    quote
        nmax=$(esc(nmax))
        ix=$(esc(ix))
        incx=$(esc(incx))
        iy=$(esc(iy))
        incy=$(esc(incy))
        n=$(esc(n))
        $(esc(dtype))<:Float64 && println("  nmax: $nmax, ix: $ix, incx: $incx, iy: $iy, incy: $incy, n: $n")
    end
end

for nmax in 20:30
for i in 1:20
for inc in [-10:-1:-1,1:10]
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)) <= nmax
    @test nel2nmax(i, inc, nmax2nel(i, inc, nmax)+1) > nmax
    @test nmax2nel(i, inc, nel2nmax(i, inc, nmax)) == nmax
    @test fast_range2args(fast_args2range(i, inc, nmax)) == (i, inc, nmax)
    @test fast_args2range(fast_range2args(i:inc:nmax)...) == i:inc:nmax
end
end
end

const TYPES = (Float64, Float32, Complex128, Complex64)

## scale by scalar
anum = 1.1
println("fast_scale!(x, ix, incx, a, n) ...")
for dtype in TYPES
for ix in (1, 2)
for incx in (1, 2)
for nmax in (FastArrayOps.NLIM_SCALE-1, FastArrayOps.NLIM_SCALE+1)
    #incx = 1
    n = nmax2nel(ix, incx, nmax)
    @printtest1(dtype, nmax, ix, incx, n)
    # INPLACE
    x0 = dtype[1:nmax]
    a = convert(dtype,anum)
    x_exp = copy(x0)
    x_exp[fast_args2range(ix, incx, n)] *= a
    
    x = copy(x0)
    fast_scale!(x,ix,incx,a,n)
    @test_approx_eq x x_exp
    
    # compare kinds
    # inc1
    if incx == 1
        x = copy(x0)
        fast_scale!(x,ix,a,n)
        @test_approx_eq x x_exp 
    end
end # for
end # for
end # for
end # for

println("fast_scale!(x, ix, incx, y, iy, incy, a, n) ...")
for dtype in TYPES
for ix in (1, 2)
for incx in (1, 2)
for nmax in (FastArrayOps.NLIM_SCALE_OOP1-1, FastArrayOps.NLIM_SCALE_OOP1+1)
    iy = ix
    incy = incx
    n = nmax2nel(ix, incx, nmax)
    @printtest2(dtype, nmax, ix, incx, iy, incy, n)
    # OUT-OF-PLACE
    y0 = dtype[1:nmax]
    y_exp = y0
    x0 = rand(dtype,nmax)
    a = convert(dtype,anum)
    x_exp = copy(x0)
    rx = fast_args2range(ix, incx, n)
    ry = fast_args2range(iy, incy, n)
    x_exp[rx] = a*y0[ry]
    
    x = copy(x0)
    y = copy(y0)
    fast_scale!(x,ix,incx,y,iy,incy,a,n)
    @test_approx_eq x x_exp
    @test_approx_eq y y_exp
    
    # compare kinds
    # inceq
    if incx == incy
        x = copy(x0)
        y = copy(y0)
        fast_scale!(x,ix,incx,y,iy,a,n)
        @test_approx_eq x x_exp
        @test_approx_eq y y_exp
        # inc1
        if incx == 1
            x = copy(x0)
            y = copy(y0)
            fast_scale!(x,ix,y,iy,a,n)
            @test_approx_eq x x_exp
            @test_approx_eq y y_exp
            # inc1ieq
            if ix == iy
                x = copy(x0)
                y = copy(y0)
                fast_scale!(x,ix,y,a,n)
                @test_approx_eq x x_exp
                @test_approx_eq y y_exp
            end
        end  
    end
    
end # for
end # for
end # for
end # for





