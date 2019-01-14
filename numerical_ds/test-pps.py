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
m = 4.5e-6#4.51e-6
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

def hd(k,y0,rks):
    # Initial conditions for the perturbations
    z = y0[1]*y0[2]/y0[3]
    dy0 = f(y0,0.0)
    dz_z = y0[3] + dy0[1]/y0[1] - dy0[3]/y0[3]
    a = -1/(z*(2.0*k)**0.5*100*k)
    b = (-dz_z*10/k - 1j*10/y0[2])*a
    return numpy.abs(a*rks[0] + b*rks[1])**2*k**3/(2*numpy.pi**2)

def rst(k,y0,rks):
    # Initial conditions for the perturbations
    z = y0[1]*y0[2]/y0[3]
    dy0 = f(y0,0.0)
    dz_z = y0[3] + dy0[1]/y0[1] - dy0[3]/y0[3]
    a = -1/(z*(2.0*k)**0.5*100*k)
    b = (-1j*10/y0[2])*a
    return numpy.abs(a*rks[0] + b*rks[1])**2*k**3/(2*numpy.pi**2)

#def horizon_exit(t, y):
#    return y[2]*y[3] - 100*k 

def solve_bg(t0,tf,start):
    # Routine to solve inflating FRW background
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=5e5)
    if start not in tevals:
        tevals = numpy.append(tevals,start)
        tevals.sort()
    startindex = numpy.where(tevals==start)
    sol = (
    scipy.integrate.odeint(f,y0,tevals,rtol=3e-14,atol=3e-14))
    ws = 1.0/sol[:,2]
    dy = numpy.array([f(soli, 0.0) for soli in sol])
    gs = 1.5*sol[:,3] + dy[:,1]/sol[:,1] - dy[:,3]/sol[:,3]
    y0bg = sol[startindex].flatten()
    return tevals, ws, gs, y0bg

def main():
    
    # Define parameters of spectrum (ic)
    start = 1e4
    finish = 8e5
    krange = numpy.logspace(5,8,300)

    # Solve background once over large range
    t0 = 1.0
    tf = 1e6
    ts, ws, gs, y0 = solve_bg(t0,tf,start)
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws,kind='linear') 
    gfit = scipy.interpolate.interp1d(ts,gs,kind='linear')
       
    # Parameters of solver
    rk = False
    rtol = 1e-4
    atol = 0.0
    t = start 

    # Header of output file
    outputf = input("Enter outputfile's name: ")
    comment = input("Any comments: ")
    with open(outputf,'w') as f:
        f.write("# KD i.c. set at t0={}\n".format(t0))
        f.write("# rtol={}, atol={}\n".format(rtol,atol))
        f.write("# Perturbations' i.c. set at start={}\n".format(start))
        f.write("# {}\n".format(comment))
        f.write("# k, Rk1, dRk1, Rk2, dRK2\n")

    for k in krange:
        with open(outputf,'a') as f:
            f.write("{} ".format(k))
        ics = numpy.array([[100*k,0],[0,10*k**2]]) 
        rks = numpy.zeros(ics.shape)
        stepstot = numpy.zeros(ics.shape[0])
        stepswkb = numpy.zeros(ics.shape[0])
        timetot = numpy.zeros(ics.shape[0])
        for i,ic in enumerate(ics):
            x0, dx0 = ic
   
            def wnew(t):
                global calls
                calls += 1 
                return k*numpy.exp(logwfit(t))
        
            def gnew(t):
                global gcalls
                gcalls += 1
                return gfit(t)
        
            starttime = time.process_time()
            solver = Solver(wnew,gnew,t=start,x=x0,dx=dx0,rtol=rtol,atol=atol)
            
            for step in solver.evolve(rk):
                t = step['t']
                x = step['x']
                dx = step['dx']
                wkb = step['wkb']
                stepstot[i] += 1
                if wkb:
                    stepswkb[i] += 1
                if t >= finish:
                    break
            
            rks[i] = x, dx 
            endtime = time.process_time()
            timetot[i] = endtime - starttime
            with open(outputf,'a') as f:
                f.write("{} {} {} {} {} ".format(x,dx,stepstot[i],stepswkb[i],timetot[i]))

            #print('calls: ',calls,gcalls)
            #print('time: ', endtime-starttime)
        power1 = hd(k,y0,rks[:,0])
        power2 = rst(k,y0,rks[:,0])
        print(k, power1, power2)
        with open(outputf,'a') as f:
            f.write("{} {}\n".format(power1,power2))

if __name__=="__main__":
    main()
