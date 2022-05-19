using Plots, DataFrames, Statistics, Measures, LaTeXStrings
# plotlyjs()


y = [1, 0.1, 0.01, 0.001]
x=-log.(y)

x1=x./2

y1 = [0.01,0.01,0.01,0.01]
x2 = [0,2,4, 7]

y3 = [0.001,10^-2]
x3 = [2.302585,2.302585]

x4 = [4.60517,4.60517]



plot(x,y, title="Tolerance vs susceptibility", yaxis=:log10, yticks=y, legend=false, xlabel="Time (hours)", ylabel="Fraction of survivors", linecolor=:red, 
    annotations=([3],[0.1],"Tolerant"), xticks=[0,2,4,6], xlims=(0,6), ylims=(1e-3,1), annotationfontsize=8, annotationcolor=:red)
plot!(x1,y, annotations=([2.4],[0.03],"Susceptible"), linecolor=:purple3, annotationfontsize=8, annotationcolor=:purple3)
plot!(x2,y1, linecolor=:black, linewidth=0.3)
plot!(x3,y3, linecolor=:black, linestyle=:dash, linewidth=0.3)
plot!(x4,y3, linecolor=:black, linestyle=:dash, linewidth=0.3)




