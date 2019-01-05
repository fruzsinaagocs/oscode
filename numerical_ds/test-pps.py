#!/usr/bin/env python
from solver import Solver
import scipy.integrate
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import numpy
import scipy
import time

mp = 1
phi_p = 23.293
m = 5e-6#4.51e-6
n = 2
calls = 0
gcalls = 0

# Necessary for getting w, g

def V(phi):
    # Inflationary potential
    return m**2*phi**n

def dV(phi):
    # Gradient of inflationary potential
    return n*m**2*phi**(n-1)

def f(y, t):
    # ODE system describing inflating FRW universe
    # y = [phi, dphi, a, H]
    return numpy.array([y[1],-3.0*y[1]*numpy.sqrt(1.0/(3*mp**2)*(0.5*y[1]**2 +
    V(y[0]))) - dV(y[0]),y[2]*y[3],(-1.0/(3*mp**2))*(y[1]**2 - V(y[0])) -
    y[3]**2])

def ic(t):
    # Initial conditions in kinetic dominance
    return numpy.array([phi_p - numpy.sqrt(2.0/3.0)*mp*numpy.log(t),
    -numpy.sqrt(2.0/3.0)*mp/t, t**(1.0/3.0), 1.0/(3.0*t)])

#def horizon_exit(t, y):
#    return y[2]*y[3] - 100*k 

def solve_bg():
    # Routine to solve inflating FRW background
    t0 = 1
    tf = 1e6
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=1e4)
    sol = (
    scipy.integrate.odeint(f,y0,tevals,rtol=3e-14,atol=3e-14))
    ws = 1.0/sol[:,2]
    dy = numpy.array([f(soli, 0.0) for soli in sol])
    gs = 1.5*sol[:,3] + dy[:,1]/sol[:,1] - dy[:,3]/sol[:,3]
    return tevals, ws, gs

def main():
    
    # Solve background once over large range
    ts, ws, gs = solve_bg()
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws) 
    gfit = scipy.interpolate.interp1d(ts,gs)
    
    # Define parameters of spectrum (ic)
    start = 1e4
    finish = 8e5
    x0, dx0 = hd()
    krange = numpy.logspace(-5,0,1000)
   
    # Parameters of solver
    rk = False
    rtol = 1e-4
    atol = 0.0
    t = start
    
    for k in krange:
        
        def wnew(t):
            global calls
            calls += 1 
            return numpy.exp(logwfit(t))
    
        def gnew(t):
            global gcalls
            gcalls += 1
            return gfit(t)
    
        starttime = time.process_time()
        solver = Solver(wnew,gnew,t=start,x=x0,dx=dx0,rtol=rtol,atol=atol)
        
        for step in solver.evolve(rk):
            if t >= finish:
                x = step['x']
                dx = step['dx']
                break
        
        endtime = time.process_time()
        print('calls: ',calls,gcalls)
        print('time: ', endtime-starttime)
    
    
    
    fig, axes = plt.subplots(1,1, sharex=False)

    # Real part of analytic and RKWKB solution
    axes.semilogx(ts[wkbs==False],numpy.real(xs[wkbs==False]),'rx')
    axes.semilogx(ts[wkbs==True],numpy.real(xs[wkbs==True]),'gx')
    axes.semilogx(tevals, sol2[:,0])
    axes.set_ylabel('$\mathcal{Re}(x)$')
   
    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
