#!/usr/bin/env python
from solver import Solver
import scipy.special
import scipy.interpolate
import sys
import matplotlib.pyplot as plt
import numpy
import scipy
import time

n=1e9
calls = 0

def g(t):
    return 0.0

def w(t):
    return numpy.sqrt(n**2-1) * 1 / (1 + t**2)

def sol(t):
    return 100*numpy.sqrt(1+t**2)/n * (1j*numpy.sin(n * numpy.arctan(t)) + numpy.cos(n * numpy.arctan(t)))

def dsol(t):
    return 100.0/numpy.sqrt(t**2+1)/n*( (t + 1j*n ) * numpy.cos(n*numpy.arctan(t)) + (1j*t - n ) * numpy.sin(n*numpy.arctan(t)))

def main():
   
    start = -1.1*n
    finish = 1.1*n
    rtol = 1e-4
    atol = 0.0

    # Simulate noisy data of log[w(t)]
    halftimes = numpy.logspace(-6, numpy.log10(finish), num=5*1e4/2)
    times = numpy.append(-1*numpy.flip(halftimes), [0])
    times = numpy.append(times, halftimes)
    logws = numpy.log(numpy.sqrt(n**2-1)*1/(1+times**2))
    logws = logws*(1+numpy.random.normal()*1e-14)
    logwfit = scipy.interpolate.interp1d(times, logws, kind=3) 
    #halftimesfit = numpy.logspace(-7, numpy.log10(finish), num=1e6)
    #timesfit = numpy.append(-1*numpy.flip(halftimesfit), halftimesfit)
    #logwsfit = logwfit(timesfit)
    #plt.plot(timesfit, logwsfit, '.')
    #plt.plot(timesfit, numpy.log(w(timesfit)), '-')
    #plt.show()
    #plt.semilogy(timesfit, abs((numpy.exp(logwsfit) - w(timesfit))/w(timesfit)))
    #plt.show()
    def wnew(t):
        global calls
        calls += 1 
        return numpy.exp(logwfit(t))
    
    starttime = time.process_time()
    #samples = numpy.random.rand(10000)*n
    #for i in range(10000):
    #   wnew(samples[i])
    #endtime = time.process_time()
    #print(wnew(samples[3000]), w(samples[3000]))
    #print('time: ', (endtime - starttime)/1e4)


    start = start/2
    finish = finish/2
    rk = False
    t = start
    x = sol(t)
    dx = dsol(t)
    
    ts, xs, dxs, wkbs, hs, oscs = [], [], [], [], [], []
    solver = Solver(wnew,g,t=t,x=x,dx=dx,rtol=rtol,atol=atol)
    
    for step in solver.evolve(rk):
        wkb = step['wkb']
        t = step['t']
        x = step['x']
        e = step['err']
        h = step['h']
        dx = step['dx']
        #if wkb:
        #    print('wkb',t,x,e,h)
        #else:
        #    print('rk',t,x,e,h)
    
        if t < finish:
            ts.append(t)
            xs.append(x)
            wkbs.append(wkb)
            dxs.append(dx)
            hs.append(h)
            oscs.append((n*numpy.arctan(t)/(2*numpy.pi))-(n*numpy.arctan(t-h)/(2*numpy.pi)))
        else:
            break
    
    ts = numpy.array(ts)
    xs = numpy.array(xs)
    dxs = numpy.array(dxs)
    wkbs = numpy.array(wkbs)
    hs = numpy.array(hs)
    oscs = numpy.array(oscs)
    endtime = time.process_time()
    print('calls: ',calls)
    print('time: ', endtime-starttime)
    
    fig, axes = plt.subplots(2,2, sharex=False)

    # Real part of analytic and RKWKB solution
    axes[0,0].plot(ts[wkbs==False],numpy.real(xs[wkbs==False]),'rx')
    axes[0,0].plot(ts[wkbs==True],numpy.real(xs[wkbs==True]),'gx')
    axes[0,0].set_ylabel('$\mathcal{Re}(x)$')
    axes[0,1].plot(ts[wkbs==False],numpy.imag(xs[wkbs==False]),'rx')
    axes[0,1].plot(ts[wkbs==True],numpy.imag(xs[wkbs==True]),'gx')
    axes[0,1].set_ylabel('$\mathcal{Im}(x)$')
    # Relative error on amplitude
    axes[1,0].semilogy(ts[wkbs==False] ,abs((sol(ts[wkbs==False])-xs[wkbs==False]))/abs(sol(ts[wkbs==False])),'rx')
    axes[1,0].semilogy(ts[wkbs==True] ,abs((sol(ts[wkbs==True])-xs[wkbs==True]))/abs(sol(ts[wkbs==True])),'gx')
    axes[1,0].set_ylabel('$|\Delta r|/|r|$')
    #axes[1,0].semilogy(ts, abs((dsol(ts)-dxs))/abs(dsol(ts)),'x')
    # Relative error on phase
    #axes[1,1].semilogy(ts[wkbs==False],abs((numpy.angle(sol(ts[wkbs==False]))-numpy.angle(xs[wkbs==False]))/numpy.angle(sol(ts[wkbs==False]))),'rx')
    #axes[1,1].semilogy(ts[wkbs==True],abs((numpy.angle(sol(ts[wkbs==True]))-numpy.angle(xs[wkbs==True]))/numpy.angle(sol(ts[wkbs==True]))),'gx')
    #axes[1,1].set_ylabel('$|\Delta \phi / |\phi|$')

    axes[1,1].plot(ts, oscs,'k-');
    axes[1,1].plot(ts[wkbs==False], oscs[wkbs==False],'rx');
    axes[1,1].plot(ts[wkbs==True], oscs[wkbs==True],'gx');

    print('\n number of WKB steps taken: ', ts[wkbs==True].size, '\n')
    print('total steps', ts.size)

    ts = numpy.linspace(ts[0],ts[-1],100000)
    axes[0,0].plot(ts,numpy.real(sol(ts)),'k-')
    axes[0,1].plot(ts, numpy.imag(sol(ts)),'k-')


    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
