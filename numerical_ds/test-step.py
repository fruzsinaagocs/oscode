#!/usr/bin/env python
from solver import Solver
import scipy.integrate
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import numpy
import scipy
import time

k=10.0
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

def horizon_exit(t, y):
    return y[2]*y[3] - 100*k 

def solve_bg():
    # Routine to solve inflating FRW background
    t0 = 1
    tf = 5e5#1e6
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=1e4)
    sol = (
    scipy.integrate.odeint(f,y0,tevals,rtol=3e-14,atol=3e-14))
    ws = k/sol[:,2]
    dy = numpy.array([f(soli, 0.0) for soli in sol])
    gs = 1.5*sol[:,3] + dy[:,1]/sol[:,1] - dy[:,3]/sol[:,3]
    return tevals, ws, gs

def main():
    
    ts, ws, gs = solve_bg()
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws,kind=5) 
    gfit = scipy.interpolate.interp1d(ts,gs,kind=5)
        
    def wnew(t):
        global calls
        calls += 1 
        return numpy.exp(logwfit(t))

    def gnew(t):
        global gcalls
        gcalls += 1
        return gfit(t)

    # For brute-force solving MS
    def F(y,t):
        return numpy.array([y[1], -wnew(t)**2*y[0]-2*gnew(t)*y[1]])

    rk = False
    start = 1.0
    t0=1e4
    finish = 3.5*1e5
    x0 = 100*k
    dx0 = 0.0
    rtol = 1e-3
    atol = 0.0
   
    if t0 > start:
        tevals = numpy.array([start, t0])
        sol1 = scipy.integrate.odeint(F,numpy.array([x0,dx0]),tevals,rtol=1e-6,atol=1e-6)
        x0,dx0 = sol1[-1,:]

    starttime = time.process_time()
    ts, xs, dxs, wkbs, hs, oscs = [], [], [], [], [], []
    solver = Solver(wnew,gnew,t=t0,x=x0,dx=dx0,rtol=rtol,atol=atol)
       
    hs = numpy.logspace(numpy.log10(0.01),numpy.log10(3e4),1000)
    rk_steps = numpy.zeros((hs.size,2),dtype=complex)
    wkb_steps = numpy.zeros((hs.size,4),dtype=complex)
    rk_errors = numpy.zeros((hs.size,2),dtype=complex) 
    wkb_errors = numpy.zeros((hs.size,4),dtype=complex)
    wkb_rerrors = numpy.zeros((hs.size,4),dtype=complex)
    sserrors = numpy.zeros((hs.size,4),dtype=complex)
    ss = numpy.zeros((hs.size,4),dtype=complex)
    
    for i,h in enumerate(hs):
        # Take a RK and WKB step
        solver.h = h 
        x_rk, dx_rk, err_rk, ws, ws5, gs, gs5 = solver.RK_step()
        solver.ws = ws
        solver.ws5 = ws5
        solver.gs = gs
        solver.gs5 = gs5
        x_wkb, dx_wkb, err_wkb, truncerr = solver.RKWKB_step()
        
        wkb_steps[i,:] = x_wkb, dx_wkb, x_wkb, dx_wkb
        wkb_errors[i,:] = truncerr[0], truncerr[1], err_wkb[0], err_wkb[1]
        wkb_rerrors[i,:] = (numpy.abs(truncerr[0])/numpy.abs(x_wkb),
        numpy.abs(truncerr[1])/numpy.abs(dx_wkb), numpy.abs(err_wkb[0])/numpy.abs(x_wkb),
        numpy.abs(err_wkb[1])/numpy.abs(dx_wkb))
        rk_steps[i,:] = x_rk, dx_rk
        rk_errors[i,:] = err_rk[:2]
        ss[i,:] = solver.rkwkbsolver4.S0(t0,t0+h)[0], solver.rkwkbsolver4.S1(h), solver.rkwkbsolver4.S2(h), solver.rkwkbsolver4.S3(h)
        sserrors[i,:] = solver.rkwkbsolver4.S0(t0,t0+h)[1], solver.rkwkbsolver4.S1(h), solver.rkwkbsolver4.Serror[2], solver.rkwkbsolver4.S3(h)

    # scipy numerical solution
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(t0+hs[-1]),num=1e4)
    sol2 =(
    scipy.integrate.odeint(F,numpy.array([x0,dx0]),tevals,rtol=1e-6,atol=1e-6))

    
    fig, axes = plt.subplots(2,2, sharex=False)

    axes[0,0].set_title('Solution')
    axes[0,0].plot(t0+hs,numpy.real(rk_steps[:,0]),'r.',alpha=0.5)
    axes[0,0].plot(t0+hs,numpy.real(wkb_steps[:,0]),'g.',alpha=0.5)
    axes[0,0].plot(tevals,sol2[:,0],'b')
    axes[0,0].set_ylim((-100.0, 100.0))

    axes[0,1].set_title('All absolute errors')
    axes[0,1].loglog(hs, numpy.abs(wkb_rerrors[:,0]),label='WKB truncerr x')
    axes[0,1].loglog(hs, numpy.abs(wkb_rerrors[:,1]),label='WKB truncerr dx')
    axes[0,1].loglog(hs, numpy.abs(wkb_errors[:,2]),label='WKB err estimate x')
    axes[0,1].loglog(hs, numpy.abs(wkb_errors[:,3]),label='WKB err estimate dx')
    axes[0,1].loglog(hs, numpy.abs(rk_errors[:,0]),label='RK error x')
    axes[0,1].loglog(hs, numpy.abs(rk_errors[:,1]),label='RK error dx')
    axes[0,1].loglog(hs, numpy.abs(sserrors[:,0]),label='S0error')
    axes[0,1].loglog(hs, numpy.abs(sserrors[:,2]),label='S2error')
    axes[0,1].legend()

    #axes[0,1].plot(hs, numpy.ones(hs.size)*rtol, label='rtol')
    #axes[0,1].loglog(hs, abs(errs_rk1/xs_rk), 'r', alpha=0.5, label='RK error')
    #axes[0,1].loglog(hs, abs(errs_wkb/xs_wkb), 'g',alpha=0.5, label='WKB error')
    
    #axes[1,0].loglog(hs, abs(truncerrs_wkb1/dxs_wkb), label='trunc. error dx')
    #axes[1,0].loglog(hs, abs(truncerrs_wkb2/xs_wkb), label='trunc error x')
    #axes[1,0].loglog(hs, abs(numpy.sqrt((truncerrs_wkb1/dxs_wkb)**2 + (truncerrs_wkb2/xs_wkb)**2)), label='trunc error x,dx')
    axes[1,0].set_title('Relative errors and rtol')
    axes[1,0].loglog(hs, numpy.ones(hs.size)*rtol,label='rtol')
    (axes[1,0].loglog(hs, numpy.max(numpy.abs(wkb_rerrors[:,:]),axis=1),
    label='max WKBerror'))
    (axes[1,0].loglog(hs, numpy.max(numpy.abs(wkb_rerrors[:,2:]),axis=1),
    label='est WKBerror'))
    (axes[1,0].loglog(hs, numpy.max(numpy.abs(rk_errors)/numpy.abs(rk_steps),
    axis=1),label='max RKerror'))
    #axes[1,0].loglog(hs, numpy.linalg.norm(wkb_errors, axis=1)/2.0, label='WKB error')
    #axes[1,0].loglog(hs, numpy.linalg.norm(rk_errors, axis=1), label='RK error dx')
    axes[1,0].legend()

    axes[1,1].set_title('$S_i(t)$')
    axes[1,1].loglog(hs, numpy.abs(ss[:,0]),label='S0')
    axes[1,1].loglog(hs, numpy.abs(ss[:,1]),label='S1')
    axes[1,1].loglog(hs, numpy.abs(ss[:,2]),label='S2')
    axes[1,1].loglog(hs, numpy.abs(ss[:,3]),label='S3')
    axes[1,1].loglog(hs, numpy.ones(hs.size)*rtol, label='rtol')
    axes[1,1].legend()
    plt.show()

if __name__=="__main__":
    main()
