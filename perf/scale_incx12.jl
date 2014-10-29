
# 

calc(::BM_TEST, x, v, ix, incx) = calc(BM_JBase1(), x, v, ix, incx)
calc(::BM_JBase1, x, v, ix, incx) = scale!(view(x, ix:incx:length(x)),v)
calc(::BM_JBase2, x, v, ix, incx) = scale!(v, view(x, ix:incx:length(x)))
calc(::BM_BLAS, x, v, ix, incx) = begin
    n = length(x)
    BLAS.scal!(div(n,incx), v, pointer_to_array(pointer(x,ix), n-ix, false), incx)
end
calc(::BM_Broadcast, x, v, ix, incx) = begin
    vx = view(x, ix:incx:length(x))
    broadcast!(.*,vx,vx,v)
end
calc(::BM_Forloop, x, v, ix, incx) = begin
    @inbounds for i=ix:incx:length(x)
        x[i] *= v
    end
end
calc(::BM_FAO, x, v, ix, incx) = unsafe_fast_scale!(x, ix, v, div(length(x),incx), incx)


typealias BenchStart (Vector{Float64}, Float64, Int, Int)
Base.start(p::BenchCase, n::Int) = (rand(n), float64(1.000001), 2, 12)
function Base.run{Op}(p::BenchCase{Op}, n::Int, s::BenchStart)  # bench type, config, start value
    x, v, ix, incx = s
    op = Op()
    calc(op, x, v, ix, incx)
    return nothing
end

