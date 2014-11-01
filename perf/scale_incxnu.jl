
# 

calc(::BM_TEST, ::FUNT, x, v, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, v, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, v, n, nel, ix, incx) = scale!(view(x, ix:incx:n),v)
calc(::BM_JBase2, ::FUNT, x, v, n, nel, ix, incx) = scale!(v, view(x, ix:incx:n))
calc(::BM_BLAS, ::FUNT, x, v, n, nel, ix, incx) = BLAS.scal!(nel, v, pointer(x,ix), incx)
calc(::BM_Broadcast, ::FUNT, x, v, n, nel, ix, incx) = begin
    vx = view(x, ix:incx:n)
    broadcast!(*,vx,vx,v)
end
calc(::BM_Forloop, ::FUNT, x, v, n, nel, ix, incx) = scaleloopgen(x, v, n, ix, incx)
function scaleloopgen(x::Array, v, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] *= v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, v, n, nel, ix, incx) = unsafe_fast_scale!(x, ix, incx, v, nel)


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), float64(1.000001), n, nmax2nel(2, 40, n), 2, 40)
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
	s::(Vector{Float64}, Float64, Int, Int, Int, Int)
    x, v, n, nel, ix, incx = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, v, n, nel, ix, incx)
    return nothing
end

