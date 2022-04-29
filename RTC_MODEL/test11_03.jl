using Plots, Dierckx

y1 = [0,5,10,15, 20]
y2 = [2,5,6.5,8,10]
x = range(0,20, length=5)

plot(x, [y1, y2])

spl = Spline1D(1:5, y1[1:5]-y2)
x0 = roots(spl)
y0 = Spline1D(1:5, y2)(x0)


spl1 = Spline1D(x, y1.-y2)
x1=roots(spl1)

