

## SAFE ARGUMENT CHECKS

macro fast_check1(x, ix, incx, n)
    quote
        0 < $(esc(incx)) || throw(ArgumentError("non-positive increment"))
        0 < $(esc(ix)) || throw(BoundsError())
        $(esc(ix))+($(esc(n))-1)*$(esc(incx)) <= length($(esc(x))) || throw(BoundsError())
    end
end
macro fast_check2(x, ix, incx, y, iy, incy, n)
    quote
        (0 != $(esc(incx)) && 0 != $(esc(incy))) || throw(ArgumentError("zero increment"))
        (0 < $(esc(ix)) && 0 < $(esc(iy))) || throw(BoundsError())
        $(esc(ix))-($(esc(n))-1)*abs($(esc(incx))) <= length($(esc(x))) || throw(BoundsError())
        $(esc(iy))-($(esc(n))-1)*abs($(esc(incy))) <= length($(esc(y))) || throw(BoundsError())
    end
end
macro fast_check3(x, ix, incx, y, iy, incy, z, iz, incz, n)
    quote
        (0 != $(esc(incx)) && 0 != $(esc(incy)) && 0 != $(esc(incz))) || throw(ArgumentError("zero increment"))
        (0 < $(esc(ix)) && 0 < $(esc(iy)) && 0 < $(esc(iz))) || throw(BoundsError())
        $(esc(ix))-($(esc(n))-1)*abs($(esc(incx))) <= length($(esc(x))) || throw(BoundsError())
        $(esc(iy))-($(esc(n))-1)*abs($(esc(incy))) <= length($(esc(y))) || throw(BoundsError())
        $(esc(iz))-($(esc(n))-1)*abs($(esc(incz))) <= length($(esc(z))) || throw(BoundsError())
    end
end

## FOR-LOOP MACROS

# scale by scalar
# x = a*x
# x = a*y
# scale by array
# x = x.*y
# x = y.*z
# add scalar
# x = x + a
# x = y + a
# add array
# x = x + y
# x = y + z
# add array times scalar
# x = x + a*y
# x = y + a*z
# copy array
# x = y


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
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
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

# x = x.*y  (same as x = y with extra x*)
macro scalearr_for(x, ix, incx, y, iy, incy, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(x))[$(esc(ix))+i*$(esc(incx))]*$(esc(y))[$(esc(iy))+i*$(esc(incy))]
        end
    end
end
macro scalearr_for_inceq(x, ix, incx, y, iy, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(x))[i]*$(esc(y))[d+i]
        end
    end
end
macro scalearr_for_inc1(x, ix, y, iy, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i]*$(esc(y))[d+i]
        end
    end
end
macro scalearr_for_inc1ieq(x, ix, y, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i]*$(esc(y))[i]
        end
    end
end

# x = y.*z 
macro scalearr_foroop(x, ix, incx, y, iy, incy, z, iz, incz, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        $(esc(incz)) < 0 && ($(esc(iz)) = $(esc(iz))+($(esc(n))-$(esc(one)))*abs($(esc(incz))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))]*$(esc(z))[$(esc(iz))+i*$(esc(incz))]
        end
    end
end
macro scalearr_foroop_inceq(x, ix, incx, y, iy, z, iz, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[dy+i]*$(esc(z))[dz+i]
        end
    end
end
macro scalearr_foroop_inc1(x, ix, y, iy, z, iz, n, one)
    quote
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[dy+i]*$(esc(z))[dz+i]
        end
    end
end
macro scalearr_foroop_inc1ieq(x, ix, y, z, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i]*$(esc(z))[i]
        end
    end
end

# x = x + a (same as x = a*x with +)
macro add_for(x, ix, incx, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] += $(esc(a))
        end
    end
end
macro add_for_inc1(x, ix, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] += $(esc(a))
        end
    end
end

# x = y + a (same as x = a*y with +)
macro add_foroop(x, ix, incx, y, iy, incy, a, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))] + $(esc(a))
        end
    end
end
macro add_foroop_inceq(x, ix, incx, y, iy, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[d+i] + $(esc(a))
        end
    end
end
macro add_foroop_inc1(x, ix, y, iy, a, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[d+i] + $(esc(a))
        end
    end
end
macro add_foroop_inc1ieq(x, ix, y, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i] + $(esc(a))
        end
    end
end

# x = x + y (same as x = x.*y with +)
macro addarr_for(x, ix, incx, y, iy, incy, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(x))[$(esc(ix))+i*$(esc(incx))] + $(esc(y))[$(esc(iy))+i*$(esc(incy))]
        end
    end
end
macro addarr_for_inceq(x, ix, incx, y, iy, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[d+i]
        end
    end
end
macro addarr_for_inc1(x, ix, y, iy, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[d+i]
        end
    end
end
macro addarr_for_inc1ieq(x, ix, y, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[i]
        end
    end
end

# x = y + z (same as x = y.*z with +)
macro addarr_foroop(x, ix, incx, y, iy, incy, z, iz, incz, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        $(esc(incz)) < 0 && ($(esc(iz)) = $(esc(iz))+($(esc(n))-$(esc(one)))*abs($(esc(incz))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))] + $(esc(z))[$(esc(iz))+i*$(esc(incz))]
        end
    end
end
macro addarr_foroop_inceq(x, ix, incx, y, iy, z, iz, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[dy+i] + $(esc(z))[dz+i]
        end
    end
end
macro addarr_foroop_inc1(x, ix, y, iy, z, iz, n, one)
    quote
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[dy+i] + $(esc(z))[dz+i]
        end
    end
end
macro addarr_foroop_inc1ieq(x, ix, y, z, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i] + $(esc(z))[i]
        end
    end
end

# x = x + a*y (same as x = a*y with extra x+)
macro addarrscal_for(x, ix, incx, y, iy, incy, a, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(x))[$(esc(ix))+i*$(esc(incx))] + $(esc(y))[$(esc(iy))+i*$(esc(incy))]*$(esc(a))
        end
    end
end
macro addarrscal_for_inceq(x, ix, incx, y, iy, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[d+i]*$(esc(a))
        end
    end
end
macro addarrscal_for_inc1(x, ix, y, iy, a, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[d+i]*$(esc(a))
        end
    end
end
macro addarrscal_for_inc1ieq(x, ix, y, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(x))[i] + $(esc(y))[i]*$(esc(a))
        end
    end
end

# x = y + a*z (same as x = y + z with extra a*)
macro addarrscal_foroop(x, ix, incx, y, iy, incy, z, iz, incz, a, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        $(esc(incz)) < 0 && ($(esc(iz)) = $(esc(iz))+($(esc(n))-$(esc(one)))*abs($(esc(incz))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(y))[$(esc(iy))+i*$(esc(incy))] + $(esc(z))[$(esc(iz))+i*$(esc(incz))]*$(esc(a))
        end
    end
end
macro addarrscal_foroop_inceq(x, ix, incx, y, iy, z, iz, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(y))[dy+i] + $(esc(z))[dz+i]*$(esc(a))
        end
    end
end
macro addarrscal_foroop_inc1(x, ix, y, iy, z, iz, a, n, one)
    quote
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[dy+i] + $(esc(z))[dz+i]*$(esc(a))
        end
    end
end
macro addarrscal_foroop_inc1ieq(x, ix, y, z, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(y))[i] + $(esc(z))[i]*$(esc(a))
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
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
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

# scale by scalar
# x = a*x       # done: scale_blas
# x = a*y       # done: copy_blas, scale_blas
# scale by array
# x = x.*y      # done, slow: vecmult_blas
# x = y.*z      # done, slow: vecmultoop_blas
# add scalar
# x = x + a     # ?
# x = y + a     # ? probably very slow: fill(x,a) axpy: x += 1*y
# add array
# x = x + y     # done: axpy_blas with a=1
# x = y + z     # done: copy_blas, axpy_blas with a=1
# add array times scalar
# x = x + a*y   # done: axpy_blas
# x = y + a*z   # done: copy_blas, axpy_blas
# copy array
# x = y         # done: copy_blas



# x = a*x
macro scale_blas(f, elty, x, ix, incx, a, n)
    quote
        #px = pointer($(esc(x)), $(esc(ix)))
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*sizeof($(esc(elty)))
        ccall(($(esc(f)),libblas), Void,
              (Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
              &($(esc(n))), &($(esc(a))), px, &($(esc(incx))))
    end
end

# x = y     
macro copy_blas(f, elty, x, ix, incx, y, iy, incy, n)
    quote
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*sizeof($(esc(elty)))
        py = convert(Ptr{$(esc(elty))},$(esc(y))) + ($(esc(iy))-1)*sizeof($(esc(elty)))
        ccall(($(esc(f)),libblas), Void,
              (Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
              &($(esc(n))), py, &($(esc(incy))), px, &($(esc(incx))))
    end
end  

# x = x + a*y
macro axpy_blas(f, elty, x, ix, incx, y, iy, incy, a, n)
    quote
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*sizeof($(esc(elty)))
        py = convert(Ptr{$(esc(elty))},$(esc(y))) + ($(esc(iy))-1)*sizeof($(esc(elty)))
        ccall(($(esc(f)),libblas), Void,
                (Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{$(esc(elty))}, Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
                 &($(esc(n))), &($(esc(a))), py, &($(esc(incy))), px, &($(esc(incx))))
    end
end

# x = y.*z
macro vecmultoop_blas(f, elty, x, ix, incx, y, iy, incy, z, iz, incz, n)
    quote
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*sizeof($(esc(elty)))
        py = convert(Ptr{$(esc(elty))},$(esc(y))) + ($(esc(iy))-1)*sizeof($(esc(elty)))
        pz = convert(Ptr{$(esc(elty))},$(esc(z))) + ($(esc(iz))-1)*sizeof($(esc(elty)))
        uplo::BlasChar = 'U'
        k::BlasInt = 0
        a::$(esc(elty)) = 1.0
        b::$(esc(elty)) = 0.0
        lda::BlasInt = $(esc(incy))
        ccall(($(esc(f)),libblas), Void,
                (Ptr{Uint8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$(esc(elty))}, 
                 Ptr{$(esc(elty))}, Ptr{BlasInt}, Ptr{$(esc(elty))}, Ptr{BlasInt}, 
                 Ptr{$(esc(elty))}, Ptr{$(esc(elty))}, Ptr{BlasInt}),
                 &uplo, &($(esc(n))), &k, &a, 
                 py, &lda, pz, &($(esc(incz))),
                 &b, px, &($(esc(incx))))
    end
end

# x = x.*y
macro vecmult_blas(f, elty, x, ix, incx, y, iy, incy, n)
    quote
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*sizeof($(esc(elty)))
        py = convert(Ptr{$(esc(elty))},$(esc(y))) + ($(esc(iy))-1)*sizeof($(esc(elty)))
        uplo::BlasChar = 'U'
        trans::BlasChar = 'N'
        diag::BlasChar = 'N'
        k::BlasInt = 0
        lda::BlasInt = $(esc(incy)) #*$(esc(n))
        ccall(($(esc(f)),libblas), Void,
                (Ptr{Uint8}, Ptr{Uint8}, Ptr{Uint8}, 
                 Ptr{BlasInt}, Ptr{BlasInt},
                 Ptr{$(esc(elty))}, Ptr{BlasInt},
                 Ptr{$(esc(elty))}, Ptr{BlasInt}),
                 &uplo, &trans, &diag, 
                 &($(esc(n))), &k, 
                 py, &lda,
                 px, &($(esc(incx))))
    end
end

