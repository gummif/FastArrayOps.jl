using FastArrayOps
using BenchmarkLite
using ArrayViews

blas_set_num_threads(1)
type BenchCase{Op, FT} <: Proc end

abstract MethodType
abstract FuncType
type BM_TEST end
type BM_JBase1 end
type BM_JBase2 end
type BM_BLAS end
type BM_Broadcast end
type BM_Forloop end
type BM_FAO end


# procedure name
Base.string{Op, FT}(::BenchCase{Op, FT}) = string("$Op"[4:end])
Base.string{FT<:FuncType}(::FT) = lowercase(string("$FT"))
# procedure codes
Base.length(p::BenchCase, n::Int) = n
Base.isvalid(p::BenchCase, n::Int) = (n > 0)
Base.done(p::BenchCase, n, s) = nothing

function runFAO(procs::Vector{Proc}, cfgs::Vector; 
             nruns::Int=0, verbose::Int=2, allowgc::Bool=true, logger::IO=STDOUT)

    bt = BenchmarkTable(cfgs, procs)
    m = length(cfgs)
    n = length(procs)
    for (j, p) in enumerate(procs)
        procname = string(p)
        verbose >= 1 && println(logger, "Benchmarking $procname ...")

        for (i, cfg) in enumerate(cfgs)
            cfgname = string(cfg)
            (nr, et) = run(p, cfg; nruns=nruns, allowgc=allowgc) #warm up
            (nr, et) = run(p, cfg; nruns=nruns, allowgc=allowgc)

            bt.nruns[i, j] = nr
            bt.etime[i, j] = et

            verbose >= 2 && 
                println(logger, "  $procname with cfg = $cfgname: nruns = $nr, elapsed = $et secs")
        end
    end
    return bt
end

