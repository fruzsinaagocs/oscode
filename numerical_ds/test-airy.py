#!/usr/bin/env python
from solver import Solver
import scipy.special
import sys
import matplotlib.pyplot as plt
import numpy
import scipy


def w(t):
    return numpy.sqrt(t)

def sol(t):
    return scipy.special.airy(-t)[0] 

def dsol(t):
    return -scipy.special.airy(-t)[1]

def main():
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    start = 0.0
    finish = 2e3
    err = 1e-4
    
    rk = False
    t = start
    x = sol(t)
    dx = dsol(t)
    
    ts, xs, wkbs = [], [], []
    solver = Solver(w,t=t,x=x,dx=dx,err=err)
    
    for step in solver.evolve(rk):
        wkb = step['wkb']
        t = step['t']
        x = step['x']
        e = step['err']
        h = step['h']
        if wkb:
            print('wkb',t,x,e,h)
        else:
            print('rk',t,x,e,h)
    
        if t < finish:
            ts.append(t)
            xs.append(x)
            wkbs.append(wkb)
        else:
            break
    
    ts = numpy.array(ts)
    xs = numpy.array(xs)
    wkbs = numpy.array(wkbs)
    RKWKB_line, = ax2.plot(ts[wkbs==False],xs[wkbs==False],'ro')
    RKWKB_line, = ax2.plot(ts[wkbs==True],xs[wkbs==True],'go')
    ax1.semilogy(ts,abs((sol(ts)-xs)/sol(ts)),'k-')
    ax1.set_ylabel('$|\Delta x|/|x|$')

    ts = numpy.linspace(ts[0],ts[-1],100000)
    true_line, = ax2.plot(ts,sol(ts),'k-')
    #ax2.set_xlim(-0.5,0.5)
    #ax2.set_ylim(-0.5,0.5)
    #ax2.set_xlabel('$t/n$')
    
    
    #rk = True
    #t = start
    #x = sol(t)
    #dx = dsol(t)
    #
    #ts, xs = [], []
    #solver = Solver(w,t=t,x=x,dx=dx,err=err)
    #
    #for step in solver.evolve(rk):
    #    wkb = step['wkb']
    #    t = step['t']
    #    x = step['x']
    #    e = step['err']
    #    h = step['h']
    #    if wkb:
    #        print('wkb',t,x,e,h)
    #    else:
    #        print('rk',t,x,e,h)
    #
    #    if t < finish:
    #        ts.append(t)
    #        xs.append(x)
    #    else:
    #        break
    
    #RK_line, = ax1.plot(numpy.array(ts)/n,xs,'rx')
        
    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
