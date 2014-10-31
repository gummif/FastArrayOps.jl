

## SAFE ARGUMENT CHECKS

macro fast_check1(x, ix, incx, n)
    quote
        0 < $(esc(incx)) || throw(ArgumentError("non-positive increment"))
        0 < $(esc(ix)) || throw(BoundsError())
        $(esc(ix))-1+$(esc(n))*$(esc(incx)) <= length($(esc(x))) || throw(BoundsError())
    end
end
macro fast_check2(x, ix, incx, y, iy, incy, n)
    quote
        (0 != $(esc(incx)) && 0 != $(esc(incy))) || throw(ArgumentError("zero increment"))
        (0 < $(esc(ix)) && 0 < $(esc(iy))) || throw(BoundsError())
        $(esc(ix))-1+$(esc(n))*abs($(esc(incx))) <= length($(esc(x))) || throw(BoundsError())
        $(esc(iy))-1+$(esc(n))*abs($(esc(incy))) <= length($(esc(y))) || throw(BoundsError())
    end
end


## FOR-LOOP MACROS

# x = a*x
macro scale_for(x, ix, incx, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))-$(esc(one))+$(esc(n))*$(esc(incx))
            $(esc(x))[i] *= $(esc(a))
        end
    end
end
macro scale_for_inc1(x, ix, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] *= $(esc(a))
        end
    end
end

# x = a*y
macro scale_foroop(x, ix, incx, y, iy, incy, a, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))]*$(esc(a))
        end
    end
end
macro scale_foroop_inceq(x, ix, incx, y, iy, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))-$(esc(one))+$(esc(n))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[d+i]*$(esc(a))
        end
    end
end
macro scale_foroop_inc1(x, ix, y, iy, a, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[d+i]*$(esc(a))
        end
    end
end
macro scale_foroop_inc1ieq(x, ix, y, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i]*$(esc(a))
        end
    end
end


# x = y  (same as x = a*y without the a*)
macro copy_foroop(x, ix, incx, y, iy, incy, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))]
        end
    end
end
macro copy_foroop_inceq(x, ix, incx, y, iy, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))-$(esc(one))+$(esc(n))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[d+i]
        end
    end
end
macro copy_foroop_inc1(x, ix, y, iy, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[d+i]
        end
    end
end
macro copy_foroop_inc1ieq(x, ix, y, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i]
        end
    end
end



## BLAS MACROS

# x = a*x
macro scale_blas(f, elty, x, ix, incx, a, n)
    quote
        px = pointer($(esc(x)), $(esc(ix)))
        ccall(($(esc(f)),libblas), Void,
              (Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
              &($(esc(n))), &($(esc(a))), px, &($(esc(incx))))
    end
end

# x = y     
macro copy_blas(f, elty, x, ix, incx, y, iy, incy, n)
    quote
        px = pointer($(esc(x)), $(esc(ix)))
        py = pointer($(esc(y)), $(esc(iy)))
        ccall(($(esc(f)),libblas), Void,
              (Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
              &($(esc(n))), py, &($(esc(incy))), px, &($(esc(incx))))
    end
end






