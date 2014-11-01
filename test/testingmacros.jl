function makesignals1(dtype, nmax)
    x0 = rand(dtype, nmax)
    x_exp = copy(x0)
    return x_exp, x0 
end
function makesignals2(dtype, nmax)
    x0 = rand(dtype, nmax)
    y0 = rand(dtype, nmax)
    x_exp = copy(x0)
    return x_exp, x0, y0
end
function makesignals3(dtype, nmax)
    x0 = rand(dtype, nmax)
    y0 = rand(dtype, nmax)
    z0 = rand(dtype, nmax)
    x_exp = copy(x0)
    return x_exp, x0, y0, z0
end


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

# 1 Array, scalar
macro testfunc_arr1a(func, x_exp, x0, ix, incx, a, n)
    quote
        func = $(esc(func))
        x_exp = $(esc(x_exp))
        x0 = $(esc(x0))
        ix = $(esc(ix))
        incx = $(esc(incx))
        a = $(esc(a))
        n = $(esc(n))

        x = copy(x0)
        func(x,ix,incx,a,n)
        @test_approx_eq x x_exp
        # compare kinds
        # inc1
        if incx == 1
            x = copy(x0)
            func(x,ix,a,n)
            @test_approx_eq x x_exp 
        end
    end
end
# 2 Arrays, scalar
macro testfunc_arr2a(func, x_exp, x0, ix, incx, y0, iy, incy, a, n)
    quote
        func = $(esc(func))
        x_exp = $(esc(x_exp))
        x0 = $(esc(x0))
        ix = $(esc(ix))
        incx = $(esc(incx))
        y0 = $(esc(y0))
        iy = $(esc(iy))
        incy = $(esc(incy))
        a = $(esc(a))
        n = $(esc(n))

        x = copy(x0)
        y = copy(y0)
        func(x,ix,incx,y,iy,incy,a,n)
        @test_approx_eq x x_exp
        @test_approx_eq y y0
        
        # compare kinds
        # inceq
        if incx == incy
            x = copy(x0)
            y = copy(y0)
            func(x,ix,incx,y,iy,a,n)
            @test_approx_eq x x_exp
            @test_approx_eq y y0
            # inc1
            if incx == 1
                x = copy(x0)
                y = copy(y0)
                func(x,ix,y,iy,a,n)
                @test_approx_eq x x_exp
                @test_approx_eq y y0
                # inc1ieq
                if ix == iy
                    x = copy(x0)
                    y = copy(y0)
                    func(x,ix,y,a,n)
                    @test_approx_eq x x_exp
                    @test_approx_eq y y0
                end
            end  
        end
    end
end
# 2 Arrays
macro testfunc_arr2(func, x_exp, x0, ix, incx, y0, iy, incy, n)
    quote
        func = $(esc(func))
        x_exp = $(esc(x_exp))
        x0 = $(esc(x0))
        ix = $(esc(ix))
        incx = $(esc(incx))
        y0 = $(esc(y0))
        iy = $(esc(iy))
        incy = $(esc(incy))
        n = $(esc(n))

        x = copy(x0)
        y = copy(y0)
        func(x,ix,incx,y,iy,incy,n)
        @test_approx_eq x x_exp
        @test_approx_eq y y0
        
        # compare kinds
        # inceq
        if incx == incy
            x = copy(x0)
            y = copy(y0)
            func(x,ix,incx,y,iy,n)
            @test_approx_eq x x_exp
            @test_approx_eq y y0
            # inc1
            if incx == 1
                x = copy(x0)
                y = copy(y0)
                func(x,ix,y,iy,n)
                @test_approx_eq x x_exp
                @test_approx_eq y y0
                # inc1ieq
                if ix == iy
                    x = copy(x0)
                    y = copy(y0)
                    func(x,ix,y,n)
                    @test_approx_eq x x_exp
                    @test_approx_eq y y0
                end
            end  
        end
    end
end


