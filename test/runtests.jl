using FastArrayOps
using Base.Test

## scale by scalar
for dtype in (Float32, Float64, Complex64, Complex128)
for n in (FastArrayOps.NLIM_SCALE-1, FastArrayOps.NLIM_SCALE+1)
    # inplace
    x = dtype[1:n]
    a = convert(dtype,2.2)
    x_exp = x*a
    fast_scale!(x,1,1,a,n)
    @test_approx_eq x x_exp
    
    # out-of-place
    y = dtype[1:n]
    y_exp = copy(y)
    x = rand(dtype,n)
    a = convert(dtype,2.2)
    x_exp = y*a
    fast_scale!(x,1,1,y,1,1,a,n)
    @test_approx_eq x x_exp
    @test_approx_eq y y_exp
end
end

