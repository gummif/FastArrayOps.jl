
# 

calc(::BM_TEST, ::FUNT, x, y, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, y, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, y, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    copy!(xv, yv)
end
calc(::BM_JBase2, ::FUNT, x, y, n, nel, ix, incx) = calc(BM_Forloop(), FUNT(), x, y, n, nel, ix, incx)
calc(::BM_BLAS, ::FUNT, x, y, n, nel, ix, incx) = begin
    n = length(x)
    BLAS.blascopy!(nel, pointer(y,ix), incx, pointer(x,ix), incx)
end
calc(::BM_Broadcast, ::FUNT, x, y, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    broadcast!(+,xv,yv)
end
calc(::BM_Forloop, ::FUNT, x, y, n, nel, ix, incx) = copyarrloopgen(x, y, n, ix, incx)
function copyarrloopgen(x::Array, y::Array, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] = y[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, n, nel, ix, incx) = unsafe_fast_copy!(x, ix, incx, y, ix, nel)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), n, nmax2nel(2, 40, n), 2, 40)

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

