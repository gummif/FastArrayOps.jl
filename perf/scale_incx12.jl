
# 

calc(::BM_TEST, ::FUNT, x, v, n, nel, ix, incx) = calc(BM_JBase1(), FUNT(), x, v, n, nel, ix, incx)
calc(::BM_JBase1, ::FUNT, x, v, n, nel, ix, incx) = scale!(view(x, ix:incx:n),v)
calc(::BM_JBase2, ::FUNT, x, v, n, nel, ix, incx) = scale!(v, view(x, ix:incx:n))
calc(::BM_BLAS, ::FUNT, x, v, n, nel, ix, incx) = BLAS.scal!(nel, v, pointer(x,ix), incx)
    #n = length(x)
    #BLAS.scal!(div(n-ix,incx)+1, v, pointer_to_array(pointer(x,ix), n-ix, false), incx)
    #BLAS.scal!(div(length(x)-ix,incx)+1, v, pointer(x,ix), incx)
#end
calc(::BM_Broadcast, ::FUNT, x, v, n, nel, ix, incx) = begin
    vx = view(x, ix:incx:n)
    broadcast!(.*,vx,vx,v)
end
calc(::BM_Forloop, ::FUNT, x, v, n, nel, ix, incx) = scaleloopgen(x, v, n, ix, incx)
function scaleloopgen(x, v, n, ix, incx)
    @inbounds for i=ix:incx:n
        x[i] *= v
    end
    return x
end
calc(::BM_FAO, ::FUNT, x, v, n, nel, ix, incx) = unsafe_fast_scale!(x, ix, v, nel, incx)
#calc(::BM_FAO, ::FUNT, x, v, ix, incx) = FastArrayOps.fast_gen_scal(x, v, n, ix, incx)


#typealias BenchStart{FUNT} (Vector{Float64}, Float64, Int, Int)
Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n), float64(1.000001), n, div(n-2,40)+1, 2, 40)
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
	s::(Vector{Float64}, Float64, Int, Int, Int, Int)
    x, v, n, nel, ix, incx = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x, v, n, nel, ix, incx)
    return nothing
end

