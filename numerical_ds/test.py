#!/usr/bin/env python
from solver import Solver
import scipy.special
import sys
import matplotlib.pyplot as plt
import numpy

n=80

def w(t):
    return numpy.sqrt(n**2-1) * 1 / (1 + t**2)

def dw(t):
    H = 10

    # order 4 estimate
    d1w4 = (-1/12.0*w(t+2*H) + 2/3.0*w(t+H) - 2/3.0*w(t-H) + 1/12.0*w(t-2*H) ) / H 

    # order 2 estimate
    d1w2 = (w(t+H/2) - w(t-H/2)) / H

    # analytic
    d1wa = numpy.sqrt(n**2-1) * -2 * t / (1 + t**2)**2
    
    return d1w4 

def ddw(t):
    H = 10
    
    # order 4 estimate
    d2w4 = (-1/12.0*w(t+2*H) + 4/3.0*w(t+H) - 5/2.0*w(t) + 4/3.0*w(t-H) - 1/12.0*w(t-2*H) ) / H**2

    # order 2 estimate
    d2w2 = (w(t+H) - 2*w(t) + w(t-H))/H**2

    # analytic
    d2wa = numpy.sqrt(n**2-1) * (6 * t**2 - 2 ) / (1 + t**2)**3

    return d2w4


def dddw(t):
    H=10

    # order 6 estimate
    d3w6 = (7/240.0*w(t+4*H) - 3/10.0*w(t+3*H) + 169/120.0*w(t+2*H) -
    61/30.0*w(t+H) + 61/30.0*w(t-H) - 169/120.0*w(t-2*H) + 3/10.0*w(t-3*H) -
    7/240.0*w(t-4*H) )/H**3

    # order 4 estimate
    d3w4 = (-1/8.0*w(t+3*H) + w(t+2*H) - 13/8.0*w(t+H) + 13/8.0*w(t-H) -
    w(t-2*H) + 1/8.0*w(t-3*H) ) / H**3

    # order 2 estimate
    d3w2 = (1/2.0*w(t+2*H) - w(t+H) + w(t-H) - 1/2.0*w(t-2*H))/H**3

    # analytic
    d3wa = numpy.sqrt(n**2-1) * 24 * t * (1 - t**2) / (1 + t**2)**4 
    
    #print(d3w1, d3w2, d3wa)
    return d3w6 

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
print('--------------------------------------------------------------------------------')

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
#print('--------------------------------------------------------------------------------')

#RK_line, = ax1.plot(numpy.array(ts)/n,xs,'rx')
ts = numpy.linspace(ts[0],ts[-1],100000)
ax1.plot(ts,sol(ts),'k-')
#ax1.set_xlim(-0.5,0.5)
#ax1.set_ylim(-0.5,0.5)

#fig.legend((true_line,RK_line,RKWKB_line),bbox_to_anchor=[0.5,0.5],labels=('$\\frac{(-1)^{n/2}}{n}\\sqrt{1+t^2}\\sin( n \\tan^{-1} t)$','pure RK method','RKWKB method'),loc='center')

plt.show()
#fig.savefig('/home/will/Documents/Papers/RKWKB/figures/burst_compare.pdf')
