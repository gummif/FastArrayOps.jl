
Benchmarks
---------

All results can bee seen in RESULTS.md.

* Benchmarks with `Float64` Arrays.
* Out-of-place functions are marked with `oop`. 
* Benchmarks specs: Julia 0.3, Linux 3.14, Intel Core i7 3612QM 2,1-3,1GHz quad core 64bit, OpenBLAS, single thread.
* For non-unit stride or start index, the `JBase1`, `JBase2`, and `Broadcast` benchmarks use ArrayViews' `view` function.
