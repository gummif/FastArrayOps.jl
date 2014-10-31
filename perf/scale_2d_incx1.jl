
# scale columns of 2d Array by first row elements

calc(::BM_TEST, ::FUNT, x) = calc(BM_JBase1(), FUNT(), x)
calc(::BM_JBase1, ::FUNT, x) = begin
    ix = 1
    nx = size(x,1)
    @inbounds for j = 1:size(x,2)
        vx = view(x, ix:ix+nx-1)
        scale!(vx, x[ix])
        ix += nx
    end
    return x
end
calc(::BM_JBase2, ::FUNT, x) = begin
    ix = 1
    nx = size(x,1)
    @inbounds for j = 1:size(x,2)
        vx = view(x, ix:ix+nx-1)
        scale!(x[ix], vx)
        ix += nx
    end
    return x
end
calc(::BM_BLAS, ::FUNT, x) = begin
    ix = 1
    nx = size(x,1)
    @inbounds for j = 1:size(x,2)
        BLAS.scal!(nx, x[ix], pointer(x,ix), 1)
        ix += nx
    end
    return x
end
calc(::BM_Broadcast, ::FUNT, x) = begin
    ix = 1
    nx = size(x,1)
    @inbounds for j = 1:size(x,2)
        vx = view(x, ix:ix+nx-1)
        broadcast!(.*, vx, vx, x[ix])
        ix += nx
    end
    return x
end
calc(::BM_Forloop, ::FUNT, x) = scaleloop1_2d(x)
function scaleloop1_2d(x::Array)
    ix = 1
    nx = size(x,1)
    @inbounds for j = 1:size(x,2)
        v = x[ix]
        for i = ix:ix+nx-1
            x[i] *= v
        end
        ix += nx
    end
    return x
end
calc(::BM_FAO, ::FUNT, x) = begin
    ix = 1
    nx = size(x,1)
    @inbounds for i = 1:size(x,2)
        unsafe_fast_scale!(x, ix, x[ix], nx)
        ix += nx
    end
    return x
end


Base.start{T}(p::BenchCase{T,FUNT}, n::Int) = (rand(n,16), )
function Base.run{Op,FUNT}(p::BenchCase{Op,FUNT}, n::Int, s)  # bench type, config, start value
    s::(Vector{Float64}, Float64, Int)
    x, = s
    op = Op()
    funt = FUNT()
    calc(op, funt, x)
    return nothing
end

