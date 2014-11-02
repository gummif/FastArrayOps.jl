
const PLOTBM = true     # save plot to img?
const PROCLEN = 8
const STDNRUNS = 5000
const STDCP = 2:PROCLEN-1
const STDCFGS = [2,4,8,12,26,50,200,1000,4000,20000,120000,600000,1000000]

include("bm_setup.jl")
if PLOTBM
    include("bm_plot.jl")
end

type Scale_incx1 <: FuncType end
type Scale_incxnu <: FuncType end
type Scale_2d_incx1 <: FuncType end
type Scale_oop_incx1 <: FuncType end
type Scale_oop_incxnu <: FuncType end
type Scalearr_incx1 <: FuncType end
type Scalearr_incxnu <: FuncType end
type Scalearr_oop_incx1 <: FuncType end
type Scalearr_oop_incxnu <: FuncType end
type Addarr_incx1 <: FuncType end
type Addarr_oop_incx1 <: FuncType end
const METHODLIST = ((Scale_incx1, "scale (ix=1, incx=1)", STDNRUNS, STDCP, STDCFGS),
                     (Scale_incxnu, "scale (ix=2, incx=40)", STDNRUNS<<2, STDCP, STDCFGS[6:end]),
                     (Scale_2d_incx1, "scale 2-d (columns)", STDNRUNS>>5, STDCP, STDCFGS),
                     (Scale_oop_incx1, "scale oop (ix,iy=1, incx,incy=1)", STDNRUNS, STDCP, STDCFGS),
                     (Scale_oop_incxnu, "scale oop (ix,iy=2, incx,incy=40)", STDNRUNS<<2, STDCP, STDCFGS[6:end]),
                     (Scalearr_incx1, "scale by Array (ix=1, incx=1)", STDNRUNS>>3, [2,4,6,7], STDCFGS),
                     (Scalearr_incxnu, "scale by Array (ix=2, incx=40)", STDNRUNS, [4,6,7], STDCFGS),
                     (Scalearr_oop_incx1, "scale by Array oop (ix=1, incx=1)", STDNRUNS>>3, [2,4,6,7], STDCFGS),
                     (Scalearr_oop_incxnu, "scale by Array oop (ix=2, incx=40)", STDNRUNS, [4,6,7], STDCFGS),
                     (Addarr_incx1, "add Array (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                     (Addarr_oop_incx1, "add Array oop (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                        )
const NM = length(METHODLIST)

global bm_method, bm_title, chooseproc, cfgs, FUNT, Nruns
for mi = NM:NM   #NM
    FUNT, bm_title, Nruns, chooseproc, cfgs = METHODLIST[mi]
    bm_method = string(FUNT())
    println(bm_title)
    include(string(bm_method,".jl"))
    include("bm_perform.jl")
    if PLOTBM
        savetableimg(rtable, bm_title, bm_method, FAOproc, chooseproc)
    end
end



