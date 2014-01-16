from math import exp,sqrt,pi

reflvl = 20
imin = -12.
imax =  12.
exact = sqrt(2*pi)

f = lambda x: exp(-x**2/2)

integral = .5*(imax-imin)*(f(imax)+f(imin))

print "%20s %20s %20s %20s"%("Ref. Step", "Integ. Val", "Exact Val","Error (%)")

for iref in range(reflvl):
    nsteps = 1 << (iref)
    dx = (imax-imin)/nsteps
    iimprove = sum( f(imin+dx/2+dx*i) for i in range(nsteps))
    integral = .5*(integral+(imax-imin)*iimprove/(nsteps))
    print "%20d %20.5f %20.5f %20.5e"%(iref,integral,exact,(exact-integral)/exact*100)
