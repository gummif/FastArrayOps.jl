
procs = Proc[ 
    BenchCase{BM_TEST,FUNT}(),
    BenchCase{BM_JBase1,FUNT}(),
    BenchCase{BM_JBase2,FUNT}(),
    BenchCase{BM_FAO,FUNT}(),
    BenchCase{BM_BLAS,FUNT}(),
    BenchCase{BM_Forloop,FUNT}(),
    BenchCase{BM_Broadcast,FUNT}(),
    BenchCase{BM_TEST,FUNT}()
    ]
@assert length(procs) == PROCLEN
FAOproc = 4
## perform the benchmark
println("Running log:")
println("--------------------")
# nruns::Int = 0, duration::Float64=1.0, allowgc::Bool=true
#rtable = run(procs, cfgs; duration=0.3)
rtable = runFAO(procs, cfgs, verbose=1, nruns=Nruns)  
@assert isa(rtable, BenchmarkTable)
println()
## show results
show(rtable; unit=:mps, cfghead="len")

