
# assume ix = 1 and incx = 1

calc(::BM_TEST, ::FUNT, x, y, v, incx) = calc(BM_JBase1(), FUNT(), x, y, v, incx)
calc(::BM_JBase1, ::FUNT, x, y, v, incx) = begin
    copy!(x, y)
    scale!(x, v)
end
calc(::BM_JBase2, ::FUNT, x, y, v, incx) = begin
    copy!(x, y)
    scale!(v, x)
end
calc(::BM_BLAS, ::FUNT, x, y, v, incx) = begin
    n = length(x)
    BLAS.blascopy!(n, y, incx, x, incx)
    BLAS.scal!(n, v, x, incx)
end
calc(::BM_Broadcast, ::FUNT, x, y, v, incx) = broadcast!(.*, x, y, v)
calc(::BM_Forloop, ::FUNT, x, y, v, incx) = scaleloop1(x, y, v)
function scaleloop1(x::Array, y::Array, v)
    @inbounds for i=1:length(x)
        x[i] = y[i]*v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, v, incx) = unsafe_fast_scale!(x, 1, y, 1, v, length(x), incx)
#calc(::BM_FAO, x, v, incx) = perfscale(length(x), v, x, incx)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), rand(n), convert(Float64,1+1e-8), 1)
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    s::(Vector{Float64}, Float64, Int)
    x, y, v, incx = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, y, v, incx)
    return nothing
end

