
const PLOTBM = true     # save plot to img?

include("bm_setup.jl")

if PLOTBM
    include("bm_plot.jl")
end

const STDCP = 2:length(procs)-1
cfgs = [2,4,8,12,16,50,200,1000,4000,20000,120000,600000]  #int(logspace(0.5,6,10))

const method_list = (("scale_incx1", "scale!", STDCP),
                     ("scale_incx12", "scale! (incx=12, ix=2)", STDCP)
                        )

for mi = 2:length(method_list)
    bm_method, bAA, chooseproc = method_list[mi]
    println(bAA)
    include(string(bm_method,".jl"))
    include("bm_perform.jl")
    if PLOTBM
        savetableimg(rtable, bAA, bm_method, FAOproc, chooseproc)
    end
end


