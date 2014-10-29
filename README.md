
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
* Element types supported are real and complex floating point. 
* `x,y,z` are Arrays of any size, `a` is a scalar, `incx` is the stride of `x`, `ix` is the starting index of `x`, `n` is the number of elements to use or modify.

```julia
# scale by scalar (x = a*x)
fast_scale!(x, ix, a, n, incx)
fast_scale!(x, ix, y, iy, a, n, incx, incy)
# scale by array (x = a*y)
fast_scale!(x, ix, y, iy, n, incx, incy)
fast_scale!(x, ix, y, iy, z, iz, n, incx, incy, incz)
# add scalar (x = x + a)
fast_add!(x, ix, a, n, incx)
fast_add!(x, ix, y, iy, a, n, incx, incy)
# add array (x = x + y)
fast_add!(x, ix, y, iy, n, incx, incy)
fast_add!(x, ix, y, iy, z, iz, n, incx, incy, incz)
# add array times constant (x = x + a*y)
fast_add!(x, ix, a, y, iy, n, incx, incy)
fast_add!(x, ix, y, iy, a, z, iz, n, incx, incy, incz)
# copy (x <- y)
fast_copy!(x, ix, y, iy, n, incx, incy)
```

Benchmarks
---------

![Scale1](/perf/scale_incx1.png)

![Scale12](/perf/scale_incx12.png)



