
# 

calc(::BM_TEST, ::FUNT, x, y, z, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, y, z, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, y, z, n, nel, ix, incx) = begin
    x[1:2] .*= y[1:2]
end
calc(::BM_JBase2, ::FUNT, x, y, z, n, nel, ix, incx) = calc(BM_Forloop(), FUNT(), x, y, z, n, nel, ix, incx)
calc(::BM_BLAS, ::FUNT, x, y, z, n, nel, ix, incx) = calc(BM_Forloop(), FUNT(), x, y, z, n, nel, ix, incx)
calc(::BM_Broadcast, ::FUNT, x, y, z, n, nel, ix, incx) = begin
    xv = view(x, ix:incx:n)
    yv = view(y, ix:incx:n)
    zv = view(z, ix:incx:n)
    broadcast!(*,xv,yv,zv)
end
calc(::BM_Forloop, ::FUNT, x, y, z, n, nel, ix, incx) = scalearrloopgen(x, y, z, n, ix, incx)
function scalearrloopgen(x::Array, y::Array, z::Array, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] = y[i]*z[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, z, n, nel, ix, incx) = unsafe_fast_scale!(x, ix, incx, y, ix, z, ix, nel)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), rand(n), n, nmax2nel(2, 40, n), 2, 40)

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

