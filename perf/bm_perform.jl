
## perform the benchmark
println("Running log:")
println("--------------------")
rtable = run(procs, cfgs; duration=0.2)
@assert isa(rtable, BenchmarkTable)
println()
## show results
show(rtable; unit=:mps, cfghead="len")

