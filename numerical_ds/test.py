#!/usr/bin/env python
from solver import Solver
import scipy.special
import sys
import matplotlib.pyplot as plt
import numpy
import scipy

n=80

def w(t):
    return numpy.sqrt(n**2-1) * 1 / (1 + t**2)

def dw(t):
    H = 1.0
    d1w = (-43/4.0*w(t) + 1536/35.0*w(t+1/4.0*H) - 16384/285.0*w(t+3/8.0*H) + 288/11.0*w(t+1/2.0*H) - 371293/87780.0*w(t+12/13.0*H) + 12/5.0*w(t+H))/H
    return d1w

def dwb(t):
    H = 1.0
    d1w = (-5/12.0*w(t-H) + 128/21.0*w(t-3/4.0*H) - 4096/285.0*w(t-5/8.0*H) + 120/11.0*w(t-1/2.0*H) - 371293/17556.0*w(t-1/13.0*H) + 284/15.0*w(t) )/H
    return d1w

def ddw(t):
    H = 1.0
    d2w = (1553/18.0*w(t) - 20736/35.0*w(t+1/4.0*H) + 794624/855.0*w(t+3/8.0*H) - 5040/11.0*w(t+1/2.0*H) + 10767497/131670.0*w(t+12/13.0*H) - 234/5.0*w(t+H))/H**2
    return d2w

def ddwb(t):
    H = 1.0
    d2w = (-269/18.0*w(t-H) + 22528/105.0*w(t-3/4.0*H) - 425984/855.0*w(t-5/8.0*H) + 4064/11.0*w(t-1/2.0*H) - 33045077/131670.0*w(t-1/13.0*H) + 2702/15.0*w(t) )/H**2
    return d2w

def dddw(t):
    H=1.0
    d3w = (- 1453/3.0*w(t) + 21248/5.0*w(t+1/4.0*H) - 2121728/285.0*w(t+3/8.0*H) + 44304/11.0*w(t+1/2.0*H) - 2599051/3135.0*w(t+12/13.0*H) + 2404/5.0*w(t+H))/H**3
    return d3w

def dddwb(t):
    H=1.0
    d3w = (-541/3.0*w(t-H) + 85248/35.0*w(t-3/4.0*H) - 1531904/285.0*w(t-5/8.0*H) + 40464/11.0*w(t-1/2.0*H) - 36015421/21945.0*w(t-1/13.0*H) + 5412/5.0*w(t) )/H**2
    return d3w

def sol(t):
    if n%2 ==0:
        return (-1)**(n/2) * numpy.sqrt(1+t**2)/n * numpy.sin(n * numpy.arctan(t))
    else:
        return (-1)**((n-1)/2) * numpy.sqrt(1+t**2)/n * numpy.cos(n * numpy.arctan(t))

def dsol(t):
    if n%2 ==0:
        return (-1)**(n/2) / numpy.sqrt(1+t**2) / n *( n * numpy.cos(n * numpy.arctan(t)) + t * numpy.sin(n * numpy.arctan(t)) )
    else:
        return (-1)**((n-1)/2) / numpy.sqrt(1+t**2) / n *( n * numpy.sin(n * numpy.arctan(t)) - t * numpy.cos(n * numpy.arctan(t)) )


def main():
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    start = -2*n
    finish = 2*n
    err = 1e-3
    
    rk = False
    t = start
    x = sol(t)
    dx = dsol(t)
    
    ts, xs, wkbs = [], [], []
    solver = Solver(w,dw,ddw,dddw,t=t,x=x,dx=dx,err=err)
    
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
    #solver = Solver(w,dw,ddw,dddw,t=t,x=x,dx=dx,err=err)
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
    ts = numpy.linspace(ts[0],ts[-1],100000)
    ax1.plot(ts,sol(ts),'k-')
    
    plt.show()
    #fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')

if __name__=="__main__":
    main()
