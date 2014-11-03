
# assume ix = 1 and incx = 1

calc(::BM_TEST, ::FUNT, x, y, z, v) = calc(BM_JBase1(), FUNT(), x, y, z, v)
calc(::BM_JBase1, ::FUNT, x, y, z, v) = begin
    copy!(x, y)
end
calc(::BM_JBase2, ::FUNT, x, y, z, v) = begin
    copy!(x, y)
end
calc(::BM_BLAS, ::FUNT, x, y, z, v) = begin
    n = length(x)
    BLAS.blascopy!(n, y, 1, x, 1)
    BLAS.axpy!(n, v, z, 1, x, 1)
end
calc(::BM_Broadcast, ::FUNT, x, y, z, v) = broadcast!(xpyv, x, y, z, v)
xpyv(x,y,v) = x + v*y
calc(::BM_Forloop, ::FUNT, x, y, z, v) = addarrloop1(x, y, z, v)
function addarrloop1(x::Array, y::Array, z::Array, v)
    @inbounds for i=1:length(x)
        x[i] = y[i] + v*z[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, z, v) = unsafe_fast_addscal!(x, 1, y, z, v, length(x))


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), rand(n), convert(Float64,1+1e-8))

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

