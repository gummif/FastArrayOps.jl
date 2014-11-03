using DataFrames
using Gadfly

function savetableimg(A::BenchmarkTable, title::String, filename::String, 
                        relproc::Int, indp = 1:size(A, 2))
    #require("DataFrames")
    #require("Gadfly")
    x = A.cfgs
    nn = size(A, 1)
    np = size(A, 2)
    dfa = Array(DataFrames.DataFrame, length(indp))
    
    for i = 1:length(indp)
        proc = indp[i]
        y = Array(Float64, nn)
        for j = 1:nn
            rel = get(A[j, relproc], BenchmarkLite.Gps())
            y[j] = rel/get(A[j, proc], BenchmarkLite.Gps())
        end
        dfa[i] = DataFrames.DataFrame(x=x, y=y, method=string(rtable.procs[proc]))
    end
    df = vcat(dfa...)
    
    p = Gadfly.plot(df, x="x", y="y", color="method", 
            Scale.x_log10, Scale.y_log10, Theme(line_width=0.7mm, default_point_size=0.5mm), Geom.line, Geom.point,
            Guide.XLabel("Array length [n]"), Guide.YLabel("Rel. [s / G ops.] to FAO"), 
            Guide.title(title))
    Gadfly.draw(SVG(string(filename,".svg"), 16cm, 10cm), p)
    Gadfly.draw(PNG(string(filename,".png"), 16cm, 10cm), p)
    return nothing
end


