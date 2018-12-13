#!/usr/bin/env python
from solver import Solver
import scipy.integrate
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import numpy
import scipy
import time

k=0.3
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

def f(t, y):
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
    tf = 1e6
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=1e4)
    sol = (
    scipy.integrate.solve_ivp(f,(t0,tf),y0,events=horizon_exit,t_eval=tevals,rtol=1e-14,atol=1e-14))
    ws = k/sol.y[2]
    dy = f(0.0, sol.y)
    gs = 1.5*sol.y[3] + dy[1]/sol.y[1] - dy[3]/sol.y[3]
    return sol.t, ws, gs
   
# Plotting interpolating functions for w, g

def plot_w_g(ts, ws, gs, logwfit, gfit):
    fig, axes = plt.subplots(1,2,sharex=False)
    axes[0].loglog(ts,ws)
    axes[0].set_title('omega')
    axes[1].semilogx(ts,gs)
    axes[1].set_title('gamma')  
    plt.show()

    logws = numpy.log(ws)
    samples = numpy.logspace(numpy.log10(start),numpy.log10(finish),num=1e3)
    logwsfit = logwfit(samples)
    plt.plot(samples, logwsfit, '.')
    plt.title('omega fit')
    plt.show()
    
    gsfit = gfit(samples)
    plt.plot(samples, gsfit, '.')
    plt.title('gamma fit')
    plt.show()

def time_w_g(wnew,gnew):

    samples = numpy.random.rand(10000)*(finish-start)+start
    starttime = time.process_time()
    for i in range(10000):
       wnew(samples[i])
    endtime = time.process_time()
    print('calling w once takes {}s'.format((endtime - starttime)/1e4))
    starttime = time.process_time()
    for i in range(10000):
        gnew(samples[i])
    endtime = time.process_time()
    print('calling g once takes {} s'.format((endtime - starttime)/1e4))
    global calls, gcalls
    calls -= 10000
    gcalls -= 10000


def main():
   
    
    ts, ws, gs = solve_bg()
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws,kind=3) 
    gfit = scipy.interpolate.interp1d(ts,gs,kind=3)
        
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
    t=start
    finish = 3.5*1e5
    x0 = 100*k
    dx0 = 0.0
    rtol = 1e-4
    atol = 0.0

    starttime = time.process_time()
    ts, xs, dxs, wkbs, hs, oscs = [], [], [], [], [], []
    solver = Solver(wnew,gnew,t=start,x=x0,dx=dx0,rtol=rtol,atol=atol)
    
    for step in solver.evolve(rk):
        wkb = step['wkb']
        t = step['t']
        x = step['x']
        e = step['err']
        h = step['h']
        dx = step['dx']
        if wkb:
            print('wkb',t,x,e,h)
        else:
            print('rk',t,x,e,h)
    
        if t < finish:
            ts.append(t)
            xs.append(x)
            wkbs.append(wkb)
            dxs.append(dx)
            hs.append(h)
            oscs.append((n*numpy.arctan(t)/(2*numpy.pi))-(n*numpy.arctan(t-h)/(2*numpy.pi)))
        else:
            break
    
    endtime = time.process_time()
    ts = numpy.array(ts)
    xs = numpy.array(xs)
    dxs = numpy.array(dxs)
    wkbs = numpy.array(wkbs)
    hs = numpy.array(hs)
    oscs = numpy.array(oscs)
    
    print('\n number of WKB steps taken: ', ts[wkbs==True].size, '\n')
    print('total steps', ts.size)
    print('calls: ',calls,gcalls)
    print('time: ', endtime-starttime)

    # Solving brute force
    tevals = numpy.logspace(numpy.log10(start),numpy.log10(finish),num=1e4)
    sol2 =(
    scipy.integrate.odeint(F,numpy.array([x0,dx0]),tevals))
    
    fig, axes = plt.subplots(1,1, sharex=False)

    # Real part of analytic and RKWKB solution
    axes.semilogx(ts[wkbs==False],numpy.real(xs[wkbs==False]),'rx')
    axes.semilogx(ts[wkbs==True],numpy.real(xs[wkbs==True]),'gx')
    axes.semilogx(tevals, sol2[:,0])
    axes.set_ylabel('$\mathcal{Re}(x)$')
    #axes[0,1].plot(ts[wkbs==False],numpy.imag(xs[wkbs==False]),'rx')
    #axes[0,1].plot(ts[wkbs==True],numpy.imag(xs[wkbs==True]),'gx')
    #axes[0,1].set_ylabel('$\mathcal{Im}(x)$')
    # Relative error on amplitude
    #axes[1,0].semilogy(ts[wkbs==False] ,abs((sol(ts[wkbs==False])-xs[wkbs==False]))/abs(sol(ts[wkbs==False])),'rx')
    #axes[1,0].semilogy(ts[wkbs==True] ,abs((sol(ts[wkbs==True])-xs[wkbs==True]))/abs(sol(ts[wkbs==True])),'gx')
    #axes[1,0].set_ylabel('$|\Delta r|/|r|$')
    #axes[1,0].semilogy(ts, abs((dsol(ts)-dxs))/abs(dsol(ts)),'x')
    # Relative error on phase
    #axes[1,1].semilogy(ts[wkbs==False],abs((numpy.angle(sol(ts[wkbs==False]))-numpy.angle(xs[wkbs==False]))/numpy.angle(sol(ts[wkbs==False]))),'rx')
    #axes[1,1].semilogy(ts[wkbs==True],abs((numpy.angle(sol(ts[wkbs==True]))-numpy.angle(xs[wkbs==True]))/numpy.angle(sol(ts[wkbs==True]))),'gx')
    #axes[1,1].set_ylabel('$|\Delta \phi / |\phi|$')

    #axes[1,1].plot(ts, oscs,'k-');
    #axes[1,1].plot(ts[wkbs==False], oscs[wkbs==False],'rx');
    #axes[1,1].plot(ts[wkbs==True], oscs[wkbs==True],'gx');

    #ts = numpy.linspace(ts[0],ts[-1],100000)
    #axes[0,0].plot(ts,numpy.real(sol(ts)),'k-')
    #axes[0,1].plot(ts, numpy.imag(sol(ts)),'k-')

    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
