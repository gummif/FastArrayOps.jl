
# assume ix = 1 and incx = 1

calc(::BM_TEST, ::FUNT, x, y, z) = calc(BM_JBase1(), FUNT(), x, y, z)
calc(::BM_JBase1, ::FUNT, x, y, z) = begin
    copy!(x, y)
    #scale!(x, v)
end
calc(::BM_JBase2, ::FUNT, x, y, z) = begin
    copy!(x, y)
    #scale!(v, x)
end
calc(::BM_BLAS, ::FUNT, x, y, z) = begin
    n = length(x)
    BLAS.blascopy!(n, y, 1, x, 1)
    BLAS.axpy!(n, 1.0, z, 1, x, 1)
end
calc(::BM_Broadcast, ::FUNT, x, y, z) = broadcast!(+, x, y, z)
calc(::BM_Forloop, ::FUNT, x, y, z) = addarrloop1(x, y, z)
function addarrloop1(x::Array, y::Array, z::Array)
    @inbounds for i=1:length(x)
        x[i] = y[i] + z[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, z) = unsafe_fast_add!(x, 1, y, z, length(x))


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), rand(n))

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

