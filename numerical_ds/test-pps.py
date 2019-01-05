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

def hd(k, y0):
    # Initial conditions for the perturbations Rk, dRk/dt
    # Hamiltonian Diagonalisation
    z = y0[2]*y0[1]/y0[3]
    dy = f(y0, 0.0)
    dz_z = y0[3] + dy[1]/y0[1] - dy[3]/y0[3]
    Rk = 1.0/(z*(2*k)**0.5)
    dRk = Rk*(-1j*k/y0[2] - dz_z)
    return (Rk, dRk)

def rst(k, y0):
    # Initial conditions for the perturbations Rk, dRk/dt
    # Renormalised Stress-Energy tensor
    z = y0[2]*y0[1]/y0[3]
    Rk = 1.0/(z*(2*k)**0.5)
    dRk = Rk*(-1j*k/y0[2])
    return (Rk, dRk)

def solve_bg(t0, tf, start):
    # Routine to solve inflating FRW background
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=1e4)
    tevals = numpy.append(tevals, start)
    tevals = numpy.sort(tevals)
    sol = (
    scipy.integrate.odeint(f,y0,tevals,rtol=3e-14,atol=3e-14))
    ws = 1.0/sol[:,2]
    dy = numpy.array([f(soli, 0.0) for soli in sol])
    gs = 1.5*sol[:,3] + dy[:,1]/sol[:,1] - dy[:,3]/sol[:,3]
    y0bg = sol[numpy.where(tevals==start)].flatten()
    print(y0bg)
    return tevals, y0bg, ws, gs

def main():
     
    # Define parameters of spectrum (ic)
    start = 1e4
    finish = 8e5
    krange = numpy.logspace(0,1,100)
   
    # Solve background once over large range
    t0 = 1
    tf = 1.1e6
    ts, y0bg, ws, gs = solve_bg(t0,tf,start)
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws) 
    gfit = scipy.interpolate.interp1d(ts,gs)
    
   
    # Parameters of solver
    rk = False
    rtol = 1e-4
    atol = 0.0
    t = start
    
    # Header of outputfile
    outputfile = input("Enter outputfile's name: ")
    with open(outputfile, 'w') as f:
        f.write("#KD initial conditions set at t={}\n".format(t0))
        f.write("#Perturbations start at t={}\n".format(start))
        f.write("#Inflaton mass, mp={}\n".format(m))
        f.write("#Initial field value, phi_p={}\n".format(phi_p))
        f.write("#Potential exponent, n={}\n".format(n))
        f.write("#Range of wavevectors: k={}--{}\n".format(krange[0],krange[-1]))
        f.write("#rtol={}\n".format(rtol))
        f.write("#atol={}\n".format(atol))
        f.write("#k, Rk, dRK\n")

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
        
        # set initial conditions
        x0, dx0 = hd(k,y0bg)
        solver = Solver(wnew,gnew,t=start,x=x0,dx=dx0,rtol=rtol,atol=atol)
        
        for step in solver.evolve(rk):
            t = step["t"]
            if t >= finish:
                x = step['x']
                dx = step['dx']
                break
        
        endtime = time.process_time()
        print('calls: ',calls,gcalls)
        print('time: ', endtime-starttime)
        #calls = 0
        #gcalls = 0
        with open(outputfile, 'a') as f:
            f.write("{} {} {}".format(k, x, dx))
            f.write('\n')
    
   
    data = numpy.genfromtxt(outputfile,dtype=complex)
    k = data[:,0]
    r = data[:,1]
    
    fig, axes = plt.subplots(1,1, sharex=False)
    axes.semilogx(k,k**3/(2*numpy.pi**2)*numpy.abs(r)**2)
    plt.show()

if __name__=="__main__":
    main()
