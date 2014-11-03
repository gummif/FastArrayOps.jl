
const PLOTBM = true     # save plot to img?
const PROCLEN = 8
const STDNRUNS = 5000000
const STDCP = 2:PROCLEN-1
const STDCFGS = [2,4,8,12,26,50,200,1000,4000,20000,120000,400000,1000000,6000000]

include("bm_setup.jl")
if PLOTBM
    include("bm_plot.jl")
end

# scale, scalar
type Scale_incx1 <: FuncType end
type Scale_incxnu <: FuncType end
type Scale_2d_incx1 <: FuncType end
type Scale_oop_incx1 <: FuncType end
type Scale_oop_incxnu <: FuncType end
# scale, array
type Scalearr_incx1 <: FuncType end
type Scalearr_incxnu <: FuncType end
type Scalearr_oop_incx1 <: FuncType end
type Scalearr_oop_incxnu <: FuncType end
# add, array
type Addarr_incx1 <: FuncType end
type Addarr_oop_incx1 <: FuncType end
# addscal
type Addscal_incx1 <: FuncType end
type Addscal_oop_incx1 <: FuncType end
# copy

const METHODLIST = ((Scale_incx1, "scale (ix=1, incx=1)", STDNRUNS, STDCP, STDCFGS),
                     (Scale_incxnu, "scale (ix=2, incx=40)", STDNRUNS<<3, STDCP, STDCFGS[6:end]),
                     (Scale_2d_incx1, "scale 2-d (columns)", STDNRUNS>>5, STDCP, STDCFGS),
                     (Scale_oop_incx1, "scale oop (ix,iy=1, incx,incy=1)", STDNRUNS, STDCP, STDCFGS),
                     (Scale_oop_incxnu, "scale oop (ix,iy=2, incx,incy=40)", STDNRUNS<<1, STDCP, STDCFGS[6:end]),
                     (Scalearr_incx1, "scale by Array (ix=1, incx=1)", STDNRUNS>>3, [2,4,6,7], STDCFGS),
                     (Scalearr_incxnu, "scale by Array (ix=2, incx=40)", STDNRUNS, [4,6,7], STDCFGS[6:end]),
                     (Scalearr_oop_incx1, "scale by Array oop (ix=1, incx=1)", STDNRUNS>>3, [2,4,6,7], STDCFGS),
                     (Scalearr_oop_incxnu, "scale by Array oop (ix=2, incx=40)", STDNRUNS, [4,6,7], STDCFGS[6:end]),
                     (Addarr_incx1, "add Array (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                     (Addarr_oop_incx1, "add Array oop (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                     (Addscal_incx1, "add Array*scalar (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                     (Addscal_oop_incx1, "add Array*scalar oop (ix=1, incx=1)", STDNRUNS>>3, [4,5,6,7], STDCFGS),
                        )
const NM = length(METHODLIST)

global bm_method, bm_title, chooseproc, cfgs, FUNT, Nruns
for mi = 1:NM   #NM
    FUNT, bm_title, Nruns, chooseproc, cfgs = METHODLIST[mi]
    bm_method = string(FUNT())
    println(bm_title)
    include(string(bm_method,".jl"))
    t_a = time()
    include("bm_perform.jl")
    t_b = time()
    println("time: ", string(round(t_b-t_a, 1)), "s\n")
    if PLOTBM
        savetableimg(rtable, bm_title, bm_method, FAOproc, chooseproc)
    end
end



