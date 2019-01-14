#!/usr/bin/env python
from solver import Solver
import scipy.integrate
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import numpy
import scipy
import time
import cProfile, pstats, io

#k=0.02
#k = 0.020327352016717128
#0.020805675382171717
#k =13530.477745798076
k = 1e5
#k = 1e7
#k = 15.3

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

def horizon_exit(t, y):
    return y[2]*y[3] - 100*k 

def solve_bg(t0,tf,start,num):
    # Routine to solve inflating FRW background
    y0 = ic(t0)
    tevals = numpy.logspace(numpy.log10(t0),numpy.log10(tf),num=num)
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
    kmax = sol[-1,2]*sol[-1,3]/100.0
    return tevals, ws, gs, y0bg, kmax
  
# Plotting interpolating functions for w, g
def plot_w_g(ts, ws, gs, wfit, gfit, start, finish):
    fig, axes = plt.subplots(1,2,sharex=False)
    axes[0].loglog(ts,ws)
    axes[0].set_title('omega')
    axes[1].semilogx(ts,gs)
    axes[1].set_title('gamma')  
    plt.show()

    logws = numpy.log(ws)
    samples = numpy.logspace(numpy.log10(2.0),numpy.log10(1e6),num=1e5)
    wsfit = wfit(samples)
    plt.semilogx(samples, wsfit,'gx')
    #plt.semilogx(ts,ws,'rx')
    plt.semilogx(ts,ws,color='red')
    plt.title('omega fit')
    plt.show()
    
    gsfit = gfit(samples)
    plt.semilogx(samples, gsfit,'gx')
    #plt.semilogx(ts,gs,'rx')
    plt.semilogx(ts,gs,color='red')
    plt.title('gamma fit')
    plt.show()

def time_w_g(wnew, gnew, start, finish):

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
  
    start = 1e4
    finish = 8e5
    t0 = 1.0
    tf = 1e6
    num = 1e5
    ts, ws, gs, y0, kmax = solve_bg(t0,tf,start,num)
    print('kmax',kmax)
    logws = numpy.log(ws)
    logwfit = scipy.interpolate.interp1d(ts,logws,kind='linear')
    gfit = scipy.interpolate.interp1d(ts,gs,kind='linear')

    def wnew(t):
        global calls
        calls += 1 
        return k*numpy.exp(logwfit(t))

    def gnew(t):
        global gcalls
        gcalls += 1
        return gfit(t)

    
    # For brute-force solving MS
    def F(y,t):
        # y = [phi, dphi, a, H]
        ybg = y[2:] 
        
        dybg = numpy.array([ybg[1],-3.0*ybg[1]*numpy.sqrt(1.0/(3*mp**2)*(0.5*ybg[1]**2 +
        V(ybg[0]))) - dV(ybg[0]),ybg[2]*ybg[3],(-1.0/(3*mp**2))*(ybg[1]**2 -
        V(ybg[0])) -
        ybg[3]**2])
        
        g = 1.5*ybg[3] + dybg[1]/ybg[1] - dybg[3]/ybg[3]
        w = k/ybg[2]
        return numpy.concatenate((numpy.array([y[1], -w**2*y[0]-2*g*y[1]]),
        dybg))

    rk = False
    t=start
    x0 = 100*k
    dx0 = 0.0
    rtol = 1e-4
    atol = 0.0
    #plot_w_g(ts, ws, gs, wfit, gfit, start, finish)
    #time_w_g(wnew, gnew, start, finish)

    # For profiling purposes
    starttime = time.process_time()
    #pr = cProfile.Profile()
    #pr.enable()


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
            #dxs.append(dx)
            #hs.append(h)
            #oscs.append((n*numpy.arctan(t)/(2*numpy.pi))-(n*numpy.arctan(t-h)/(2*numpy.pi)))
            pass
        else:
            break
    
    #pr.disable()
    endtime = time.process_time()
    ts = numpy.array(ts)
    xs = numpy.array(xs)
    #dxs = numpy.array(dxs)
    wkbs = numpy.array(wkbs)
    #hs = numpy.array(hs)
    #oscs = numpy.array(oscs)
    
    #print('\n number of WKB steps taken: ', ts[wkbs==True].size, '\n')
    #print('total steps', ts.size)
    print('calls: ',calls,gcalls)
    print('time: ', endtime-starttime)
    #s = io.StringIO()
    #sortby = 'cumulative'
    #ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
    #ps.print_stats()
    #print(s.getvalue())

    # Solving brute force
    #tevals = numpy.logspace(numpy.log10(start),numpy.log10(finish),num=1e4)
    #sol2 =(
    #scipy.integrate.odeint(F,numpy.concatenate((numpy.array([x0,dx0]),y0)),tevals,rtol=1e-8,atol=1e-10))
    
    fig, axes = plt.subplots(1,1, sharex=False)

    # Real part of analytic and RKWKB solution
    axes.semilogx(ts[wkbs==False],numpy.real(xs[wkbs==False]),'rx')
    axes.semilogx(ts[wkbs==True],numpy.real(xs[wkbs==True]),'gx')
    #axes.semilogx(tevals, sol2[:,0])
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
