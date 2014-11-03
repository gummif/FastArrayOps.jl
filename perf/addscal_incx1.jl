
# assume ix = 1 and incx = 1

calc(::BM_TEST, ::FUNT, x, y, v) = calc(BM_JBase1(), FUNT(), x, y, v)
calc(::BM_JBase1, ::FUNT, x, y, v) = begin
    #scale!(reshape(x,(1,length(x))),y)
    copy!(x, y)
end
calc(::BM_JBase2, ::FUNT, x, y, v) = begin
    copy!(x, y)
    #scale!(v, x)
end
calc(::BM_BLAS, ::FUNT, x, y, v) = begin
    n = length(x)
    BLAS.axpy!(n, v, y, 1, x, 1)
end
calc(::BM_Broadcast, ::FUNT, x, y, v) = broadcast!(xpyv, x, x, y, v)
xpyv(x,y,v) = x + v*y
calc(::BM_Forloop, ::FUNT, x, y, v) = addarrloop1(x, y, v)
function addarrloop1(x::Array, y::Array, v)
    @inbounds for i=1:length(x)
        x[i] += v*y[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y, v) = unsafe_fast_addscal!(x, 1, y, v, length(x))


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), rand(n), convert(Float64,1+1e-8))

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

