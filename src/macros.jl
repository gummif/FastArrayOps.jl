

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

# name convention, binary operations: 
# op(a, x) = scalarr1, op(x, a) = arr1scal
# op(x, y) = arr2xy, op(y, x) = arr2yx
# foroop for same variable on both sides of =

# x = op(a, x)
macro scalarr1_for(op, x, ix, incx, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))-$(esc(one))+$(esc(n))*$(esc(incx))
            $(esc(x))[i] = $(esc(op))( $(esc(a)), $(esc(x))[i] )
        end
    end
end
macro scalarr1_for_inc1(op, x, ix, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(a)), $(esc(x))[i] )
        end
    end
end

# x = op(a, y)
macro scalarr1_foroop(op, x, ix, incx, y, iy, incy, a, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(op))( $(esc(a)), $(esc(y))[$(esc(iy))+i*$(esc(incy))] )
        end 
    end
end
macro scalarr1_foroop_inceq(op, x, ix, incx, y, iy, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(op))( $(esc(a)), $(esc(y))[d+i] )
        end
    end
end
macro scalarr1_foroop_inc1(op, x, ix, y, iy, a, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(a)), $(esc(y))[d+i] )
        end
    end
end
macro scalarr1_foroop_inc1ieq(op, x, ix, y, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(a)), $(esc(y))[i] )
        end
    end
end

# x = op(x, y)
macro arr2xy_for(op, x, ix, incx, y, iy, incy, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(op))( $(esc(x))[$(esc(ix))+i*$(esc(incx))], $(esc(y))[$(esc(iy))+i*$(esc(incy))] )
        end 
    end
end
macro arr2xy_for_inceq(op, x, ix, incx, y, iy, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(op))( $(esc(x))[i], $(esc(y))[d+i] )
        end
    end
end
macro arr2xy_for_inc1(op, x, ix, y, iy, n, one)
    quote
        d = $(esc(iy)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(x))[i], $(esc(y))[d+i] )
        end
    end
end
macro arr2xy_for_inc1ieq(op, x, ix, y, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(x))[i], $(esc(y))[i] )
        end
    end
end

# x = op(y, z)
macro arr2yz_foroop(op, x, ix, incx, y, iy, incy, z, iz, incz, n, zero, one)
    quote
        $(esc(incx)) < 0 && ($(esc(ix)) = $(esc(ix))+($(esc(n))-$(esc(one)))*abs($(esc(incx))))
        $(esc(incy)) < 0 && ($(esc(iy)) = $(esc(iy))+($(esc(n))-$(esc(one)))*abs($(esc(incy))))
        $(esc(incz)) < 0 && ($(esc(iz)) = $(esc(iz))+($(esc(n))-$(esc(one)))*abs($(esc(incz))))
        @inbounds for i = $(esc(zero)):$(esc(n))-$(esc(one))
            $(esc(x))[$(esc(ix))+i*$(esc(incx))] = $(esc(op))( $(esc(y))[$(esc(iy))+i*$(esc(incy))], $(esc(z))[$(esc(iz))+i*$(esc(incz))] )
        end
    end
end
macro arr2yz_foroop_inceq(op, x, ix, incx, y, iy, z, iz, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))+($(esc(n))-$(esc(one)))*$(esc(incx))
            $(esc(x))[i] = $(esc(op))( $(esc(y))[dy+i], $(esc(z))[dz+i] )
        end
    end
end
macro arr2yz_foroop_inc1(op, x, ix, y, iy, z, iz, n, one)
    quote
        dy = $(esc(iy)) - $(esc(ix))
        dz = $(esc(iz)) - $(esc(ix))
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(y))[dy+i], $(esc(z))[dz+i] )
        end
    end
end
macro arr2yz_foroop_inc1ieq(op, x, ix, y, z, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(op))( $(esc(y))[i], $(esc(z))[i] )
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


# x = y  (same as x = op(a, y) without the op(a, ))
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
# for inc1, inc1ieq use memcpy_c
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

# x = a
macro fill_for(x, ix, incx, a, n, one)
    quote
        $(esc(incx)) = abs($(esc(incx)))
        @inbounds for i = $(esc(ix)):$(esc(incx)):$(esc(ix))-$(esc(one))+$(esc(n))*$(esc(incx))
            $(esc(x))[i] = $(esc(a))
        end
    end
end
# for inc1 use memset_c
macro fill_for_inc1(x, ix, a, n, one)
    quote
        @inbounds for i = $(esc(ix)):$(esc(ix))-$(esc(one))+$(esc(n))
            $(esc(x))[i] = $(esc(a))
        end
    end
end



## C LIBRARY MACROS

# x = y
macro memcpy_c(elty, x, ix, y, iy, n)
    quote
        selty = sizeof($(esc(elty)))
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*selty
        py = convert(Ptr{$(esc(elty))},$(esc(y))) + ($(esc(iy))-1)*selty
        ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Uint),
          px, py, $(esc(n))*selty)
    end
end

# x = a
macro memset_c(elty, x, ix, a, n)
    quote
        a::Int32 = $(esc(a))
        selty = sizeof($(esc(elty)))
        px = convert(Ptr{$(esc(elty))},$(esc(x))) + ($(esc(ix))-1)*selty
        ccall(:memset, Ptr{Void}, (Ptr{Void}, Int32, Csize_t),
          px, a, $(esc(n))*selty)
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

