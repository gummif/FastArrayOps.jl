using FastArrayOps
using BenchmarkLite
using ArrayViews

if !isdefined(:BenchCase)
    type BenchCase{Op} <: Proc end
end
type BM_TEST end
type BM_JBase1 end
type BM_JBase2 end
type BM_BLAS end
type BM_Broadcast end
type BM_Forloop end
type BM_FAO end

procs = Proc[ 
    BenchCase{BM_TEST}(),
    BenchCase{BM_JBase1}(),
    BenchCase{BM_JBase2}(),
    BenchCase{BM_FAO}(),
    BenchCase{BM_BLAS}(),
    BenchCase{BM_Forloop}(),
    BenchCase{BM_Broadcast}(),
    BenchCase{BM_TEST}()
    ]
FAOproc = 4
# procedure name
Base.string{Op}(::BenchCase{Op}) = string("$Op"[4:end])
# procedure codes
Base.length(p::BenchCase, n::Int) = n
Base.isvalid(p::BenchCase, n::Int) = (n > 0)
Base.done(p::BenchCase, n, s) = nothing

