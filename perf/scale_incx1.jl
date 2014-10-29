
# assume incx == 1

calc(::BM_TEST, ::FUNT, x, v, incx) = calc(BM_JBase1(), FUNT(), x, v, incx)
calc(::BM_JBase1, ::FUNT, x, v, incx) = scale!(x, v)
calc(::BM_JBase2, ::FUNT, x, v, incx) = scale!(v, x)
calc(::BM_BLAS, ::FUNT, x, v, incx) = BLAS.scal!(length(x), v, x, incx)
calc(::BM_Broadcast, ::FUNT, x, v, incx) = broadcast!(.*, x, x, v)
calc(::BM_Forloop, ::FUNT, x, v, incx) = scaleloop1(x, v, length(x))
function scaleloop1(x, v, n)
    @inbounds for i=1:n
        x[i] *= v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, v, incx) = unsafe_fast_scale!(x, 1, v, length(x), incx)
#calc(::BM_FAO, x, v, incx) = perfscale(length(x), v, x, incx)

function perfscale(n, v, x, incx)
    if n<13 && incx == 1
        @inbounds for i = 1:n
            x[i] *= v
        end
    else
        unsafe_fast_scale!(x, 1, v, n, incx)
    end
    return x
end

#typealias BenchStart{FUNT} (Vector{Float64}, Float64, Int)
Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), convert(Float64,1+1e-8), 1) # x, v, incx
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    s::(Vector{Float64}, Float64, Int)
    x, v, incx = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, v, incx)
    return nothing
end

