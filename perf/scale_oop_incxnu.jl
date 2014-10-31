
# 

calc(::BM_TEST, ::FUNT, x, y, v, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, y, v, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, y, v, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    copy!(xv, yv)
    scale!(xv, v)
end
calc(::BM_JBase2, ::FUNT, x, y, v, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    copy!(xv, yv)
    scale!(v, xv)
end
calc(::BM_BLAS, ::FUNT, x, y, v, n, nel, ix, incx) = begin
    BLAS.blascopy!(nel, pointer(y,ix), incx, pointer(x,ix), incx)
    BLAS.scal!(nel, v, pointer(x,ix), incx)
end
calc(::BM_Broadcast, ::FUNT, x, y, v, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    broadcast!(.*,xv,yv,v)
end
calc(::BM_Forloop, ::FUNT, x, y, v, n, nel, ix, incx) = scaleloopgen(x, y, v, n, ix, incx)
function scaleloopgen(x::Array, y::Array, v, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] = y[i]*v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, v, n, nel, ix, incx) = unsafe_fast_scale!(x, ix, incx, y, ix, v, nel)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), rand(n), float64(1.000001), n, nmax2nel(2, 40, n), 2, 40)
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    s::(Vector{Float64}, Float64, Int, Int, Int, Int)
    x, y, v, n, nel, ix, incx = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, y, v, n, nel, ix, incx)
    return nothing
end

