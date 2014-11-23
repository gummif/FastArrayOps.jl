
const IDLIST = (
# SAFE ARGUMENT CHECKS
    ("fast_check1", false),
    ("fast_check2", false),
    ("fast_check3", false),
    
# FOR-LOOP MACROS
    # x = op(a, x)
    ("scalarr1_for", true),
    ("scalarr1_for_inc1", true),
    # x = op(a, y)
    ("scalarr1_foroop", true),
    ("scalarr1_foroop_inceq", true),
    ("scalarr1_foroop_inc1", true),
    ("scalarr1_foroop_inc1ieq", true),
    # x = op(x, y)
    ("arr2xy_for", true),
    ("arr2xy_for_inceq", true),
    ("arr2xy_for_inc1", true),
    ("arr2xy_for_inc1ieq", true),
    # x = op(y, z)
    ("arr2yz_foroop", true),
    ("arr2yz_foroop_inceq", true),
    ("arr2yz_foroop_inc1", true),
    ("arr2yz_foroop_inc1ieq", true),
    # x = x + a*y
    ("addarrscal_for", false),
    ("addarrscal_for_inceq", false),
    ("addarrscal_for_inc1", false),
    ("addarrscal_for_inc1ieq", false),
    # x = y + a*z
    ("addarrscal_foroop", false),
    ("addarrscal_foroop_inceq", false),
    ("addarrscal_foroop_inc1", false),
    ("addarrscal_foroop_inc1ieq", false),
    # x = y
    ("copy_foroop", false),
    ("copy_foroop_inceq", false),
    ("copy_foroop_inc1", false),
    ("copy_foroop_inc1ieq", false),
    # x = a
    ("fill_foroop", false),
    ("fill_foroop_inc1", false),
    
# C LIBRARY MACROS
    # x = y
    ("c_memcpy", false),
    # x = a
    ("c_memset", false),
    
# BLAS MACROS
    # x = a*x
    ("blas_scale", false),
    # x = y
    ("blas_copy", false),
    # x = x + a*y
    ("blas_axpy", false),
    
# SET MACROS
    ("set1_inc1", false),
    ("set2_inceq", false),
    ("set2_inc1", false),
    ("set2_inc1ieq", false),
    ("set3_inceq", false),
    ("set3_inc1", false),
    ("set3_inc1ieq", false),
)

# name convention, binary operations: 
# op(a, x) = scalarr1, op(x, a) = arr1scal
# op(x, y) = arr2xy, op(y, x) = arr2yx
# foroop for same variable on both sides of =

# "compile" the source srcfile to destfile
function makeFAO()
    destfile = "FastArrayOps.jl"   # warning: gets overwritten
    srcfile = "FastArrayOps_src.jl"
    dir = "templates/"
    opid = "\$op"
    preid = "##!"
    postid = "\n"
    fext = ".jl"
    oplist = ("*","/","+","-")

    file = open(destfile, "w")
    f = readall(srcfile)

    for (id, opstr) in IDLIST
        idfs = readall(string(dir, id, fext))
        idfs = string("# src: ", string(id, fext), "\n", idfs)
        if opstr
            for op in oplist
                idfsop = replace(idfs, opid, op)
                idstr = string(preid, id, "#", op, postid)
                if searchindex(f, idstr) != 0
                    f = replace(f, idstr, idfsop)
                end
            end
        else
            f = replace(f, string(preid, id, postid), idfs)
        end
    end

    idmiss = searchindex(f, preid)
    idmiss != 0 && error("missing ids found : ", f[idmiss:idmiss+30])

    write(file, f)
    close(file)
end

makeFAO()

