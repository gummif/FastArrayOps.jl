
# assume ix = 1 and incx = 1, fill 1

calc(::BM_TEST, ::FUNT, x, v) = calc(BM_JBase1(), FUNT(), x, v)
calc(::BM_JBase1, ::FUNT, x, v) = fill!(x, v)
calc(::BM_JBase2, ::FUNT, x, v) = calc(BM_JBase1(), FUNT(), x, v)
calc(::BM_BLAS, ::FUNT, x, v) = calc(BM_JBase1(), FUNT(), x, v)
calc(::BM_Broadcast, ::FUNT, x, v) = broadcast!(+, x, v)
calc(::BM_Forloop, ::FUNT, x, v) = fillloop1(x, v)
function fillloop1(x::Array, v)
    @inbounds for i=1:length(x)
        x[i] = v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, v) = unsafe_fast_fill!(x, 1, v, length(x))


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = 
    (rand(n), convert(Float64,0))

function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    calc(Op(), FUNT(), s...)
    return nothing
end

