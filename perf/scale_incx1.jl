
# assume incx == 1

calc(::BM_TEST, x, v, incx) = calc(BM_JBase1(), x, v, incx)
calc(::BM_JBase1, x, v, incx) = scale!(x, v)
calc(::BM_JBase2, x, v, incx) = scale!(v, x)
calc(::BM_BLAS, x, v, incx) = BLAS.scal!(length(x), v, x, incx)
calc(::BM_Broadcast, x, v, incx) = broadcast!(.*, x, x, v)
calc(::BM_Forloop, x, v, incx) = begin
    @inbounds for i=1:length(x)
        x[i] *= v
    end
    return x
end
calc(::BM_FAO, x, v, incx) = unsafe_fast_scale!(x, 1, v, length(x), incx)
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

typealias BenchStart (Vector{Float64}, Float64, Int)
Base.start(p::BenchCase, n::Int) = (rand(n), convert(Float64,1+1e-8), 1) # x, v, incx
function Base.run{Op}(p::BenchCase{Op}, n::Int, s::BenchStart)  # bench type, config, start value
    x, v, incx = s
    op = Op()
    calc(op, x, v, incx)
    return nothing
end

