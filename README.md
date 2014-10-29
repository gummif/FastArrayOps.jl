
FastArrayOps.jl
===============

A [Julia](https://github.com/JuliaLang/julia) package for fast inplace Array operations.

Contains a large collection of benchmarks in `/perf`. They (should) show that the functions in this package provide the fastest way of performing the specified operation in julia, on average over a large range of parameters.

Different or new implementations or benchmark cases are welcome as issues or pull requests.

Currently implemented:

* fast_scale!(x, ix, a, n, incx)
* fast_scale!(x, ix, y, iy, a, n, incx, incy)

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

* All functions are inplace (or out-of-place), overwriting the first argument.
* No Array allocations.
* Element types supported are real and complex floating point numbers. 
* `x,y,z` are Arrays of any size, `a` is a scalar, `incx` is the stride of `x`, `ix` is the starting index of `x`, `n` is the number of elements to use or modify. A negative stride `incx` reverses the indexing order.
* All functions `f` have an `unsafe_f` version without any argument or bounds checks (e.g. `unsafe_fast_scale!`).

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
# add array times scalar
fast_add!(x, ix, a, y, iy, n, incx, incy)               # x = x + a*y
fast_add!(x, ix, y, iy, a, z, iz, n, incx, incy, incz)  # x = y + a*z
# copy array
fast_copy!(x, ix, y, iy, n, incx, incy)                 # x <- y
```

Benchmarks
---------

Benchmarks with `Float64` Arrays on Intel Core i7 3612QM 2,1-3,1GHz quad core 64bit. See directory `/perf` for details. Out-of-place functions are marked with `oop`.

![Scale1](/perf/scale_incx1.png)

![Scale12](/perf/scale_incxnu.png)

![Scale1](/perf/scale_oop_incx1.png)

![Scale12](/perf/scale_oop_incxnu.png)



