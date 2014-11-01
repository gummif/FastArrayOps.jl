
# assume ix = 1 and incx = 1

calc(::BM_TEST, ::FUNT, x, y) = calc(BM_JBase1(), FUNT(), x, y)
calc(::BM_JBase1, ::FUNT, x, y) = begin
    scale!(reshape(x,(1,length(x))),y)
end
calc(::BM_JBase2, ::FUNT, x, y) = begin
    copy!(x, y)
    #scale!(v, x)
end
calc(::BM_BLAS, ::FUNT, x, y) = begin
    n = length(x)
    BLAS.blascopy!(n, y, 1, x, 1)
    #BLAS.scal!(n, v, x, 1)
end
calc(::BM_Broadcast, ::FUNT, x, y) = broadcast!(*, x, x, y)
calc(::BM_Forloop, ::FUNT, x, y) = scalearrloop1(x, y)
function scalearrloop1(x::Array, y::Array)
    @inbounds for i=1:length(x)
        x[i] *= y[i]
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, y) = unsafe_fast_scale!(x, 1, y, length(x))


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), rand(n))
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    s::(Vector{Float64}, Float64)
    x, y = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, y)
    return nothing
end

