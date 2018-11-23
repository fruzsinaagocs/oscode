#!/usr/bin/env python
from solver import Solver
import scipy.special
import sys
import matplotlib.pyplot as plt
import numpy
import scipy

#scipy.special.seterr(all='warn')

def w(t):
    return numpy.sqrt(t)

def sol(t):
    return scipy.special.airy(-t)[0] + 1j*scipy.special.airy(-t)[2] 

def dsol(t):
    return -scipy.special.airy(-t)[1] - 1j*scipy.special.airy(-t)[3]

def main():
    
    
    start = 1.0
    finish = 1e5
    rtol = 1e-6
    atol = 1e-5

    rk = False
    t = start
    x = sol(t)
    dx = dsol(t)
    
    ts, xs, wkbs = [], [], []
    solver = Solver(w,t=t,x=x,dx=dx,rtol=rtol,atol=atol)
    
    for step in solver.evolve(rk):
        wkb = step['wkb']
        t = step['t']
        x = step['x']
        e = step['err']
        h = step['h']
        if wkb:
            print('wkb',t,x,e,h)
            #print('analytic S2:', 5/48.0*(1/(t+h)**(3/2.0) - 1/t**(3/2.0)))
        else:
            print('rk',t,x,e,h)
            #print('analytic S2:', 5/48.0*(1/(t+h)**(3/2.0) - 1/t**(3/2.0)))
 
        if t < finish:
            ts.append(t)
            xs.append(x)
            wkbs.append(wkb)
        else:
            break
    
    ts = numpy.array(ts)
    xs = numpy.array(xs)
    wkbs = numpy.array(wkbs)
    
    fig, axes = plt.subplots(2,2,sharex=True)

    axes[0,0].semilogx(ts[wkbs==False],numpy.real(xs[wkbs==False]),'rx')
    axes[0,0].semilogx(ts[wkbs==True],numpy.real(xs[wkbs==True]),'gx')
    axes[0,0].set_ylabel('$\mathcal{Re}(x)$')
    axes[0,1].semilogx(ts[wkbs==False],numpy.imag(xs[wkbs==False]),'rx')
    axes[0,1].semilogx(ts[wkbs==True],numpy.imag(xs[wkbs==True]),'gx')
    axes[0,1].set_ylabel('$\mathcal{Im}(x)$')
    # Relative error on amplitude
    axes[1,0].loglog(ts[wkbs==False] ,abs((sol(ts[wkbs==False])-xs[wkbs==False]))/abs(sol(ts[wkbs==False])),'rx')
    axes[1,0].loglog(ts[wkbs==True] ,abs((sol(ts[wkbs==True])-xs[wkbs==True]))/abs(sol(ts[wkbs==True])),'gx')
    axes[1,0].set_ylabel('$|\Delta r|/|r|$')
    # Relative error on phase
    axes[1,1].loglog(ts[wkbs==False],abs((numpy.angle(sol(ts[wkbs==False]))-numpy.angle(xs[wkbs==False]))/numpy.angle(sol(ts[wkbs==False]))),'rx')
    axes[1,1].loglog(ts[wkbs==True],abs((numpy.angle(sol(ts[wkbs==True]))-numpy.angle(xs[wkbs==True]))/numpy.angle(sol(ts[wkbs==True]))),'gx')
    axes[1,1].set_ylabel('$|\Delta \phi / |\phi|$')

    ts = numpy.linspace(ts[0],ts[-1],100000)
    axes[0,0].plot(ts,numpy.real(sol(ts)),'k-')
    axes[0,1].plot(ts, numpy.imag(sol(ts)),'k-')

    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
