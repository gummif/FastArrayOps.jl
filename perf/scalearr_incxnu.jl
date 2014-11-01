
# 

calc(::BM_TEST, ::FUNT, x, y, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, y, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, y, n, nel, ix, incx) = begin
    x[1:2] .*= y[1:2]
    #xx=reshape(x,(1,length(x)))
    #xv = view(xx, 1, ix:incx:n)
    #yv = view(y, ix:incx:n)
    #scale!(xv,yv)
end
calc(::BM_JBase2, ::FUNT, x, y, n, nel, ix, incx) = calc(BM_Forloop(), FUNT(), x, y, n, nel, ix, incx)
calc(::BM_BLAS, ::FUNT, x, y, n, nel, ix, incx) = calc(BM_Forloop(), FUNT(), x, y, n, nel, ix, incx)
calc(::BM_Broadcast, ::FUNT, x, y, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    broadcast!(*,xv,xv,yv)
end
calc(::BM_Forloop, ::FUNT, x, y, n, nel, ix, incx) = scalearrloopgen(x, y, n, ix, incx)
function scalearrloopgen(x::Array, y::Array, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] *= y[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, n, nel, ix, incx) = unsafe_fast_scale!(x, ix, incx, y, ix, nel)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), n, nmax2nel(2, 40, n), 2, 40)

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

