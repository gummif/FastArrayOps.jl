
FastArrayOps.jl
---------

A [Julia](https://github.com/JuliaLang/julia) package for fast inplace Array operations.

See license (MIT) in LICENSE.md.


Usage
---------

Install via the package manager and load with `using`

```julia
julia> Pkg.add("FastArrayOps")
julia> using FastArrayOps
```

API
---------

* All functions are inplace, overwriting the first argument. 
* No Array allocations
* Element types supported are real and complex floating point numbers. 
* `x,y,z` are Arrays of any size, `a` is a scalar, `incx` is the stride of `x`, `ix` is the starting index of `x`, `n` is the number of elements to use or modify. A negative stride `incx` reverses the indexing order.

```julia
# scale by scalar
fast_scale!(x, ix, a, n, incx)                          # x = a*x
fast_scale!(x, ix, y, iy, a, n, incx, incy)             # x = a*y
# scale by array
fast_scale!(x, ix, y, iy, n, incx, incy)                # x = x.*y
fast_scale!(x, ix, y, iy, z, iz, n, incx, incy, incz)   # x = y.*z
# add scalar
fast_add!(x, ix, a, n, incx)                            # x = x + a
fast_add!(x, ix, y, iy, a, n, incx, incy)               # x = y + a
# add array
fast_add!(x, ix, y, iy, n, incx, incy)                  # x = x + y
fast_add!(x, ix, y, iy, z, iz, n, incx, incy, incz)     # x = y + z
# add array times constant
fast_add!(x, ix, a, y, iy, n, incx, incy)               # x = x + a*y
fast_add!(x, ix, y, iy, a, z, iz, n, incx, incy, incz)  # x = y + a*z
# copy array
fast_copy!(x, ix, y, iy, n, incx, incy)                 # x <- y
```

Benchmarks
---------

![Scale1](/perf/scale_incx1.png)

![Scale12](/perf/scale_incx12.png)



