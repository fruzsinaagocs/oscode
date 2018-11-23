#!/usr/bin/env python
from solver import Solver
import scipy.special
import sys
import matplotlib.pyplot as plt
import numpy
import scipy

n=115

def w(t):
    return (n**2-1)**0.5 * 1.0 / (1 + t**2)

def sol(t):
    return 100*numpy.sqrt(1+t**2)/n * (1j*numpy.sin(n * numpy.arctan(t)) + numpy.cos(n * numpy.arctan(t)))

def dsol(t):
    return 100.0/numpy.sqrt(t**2+1)/n*( (t + 1j*n ) * numpy.cos(n*numpy.arctan(t)) + (1j*t - n ) * numpy.sin(n*numpy.arctan(t)))

def main():
    
        
    start = -2*n
    finish = 2*n
    rtol = 1e-4
    atol = 0.0
    
    rk = False
    t = 6.0#-19.125
    x = sol(t)
    dx = dsol(t)
   
    
    hs = numpy.logspace(numpy.log10(0.0001),numpy.log10(10.0),1000)
    rk_steps = numpy.zeros((hs.size,2),dtype=complex)
    wkb_steps = numpy.zeros((hs.size,4),dtype=complex)
    rk_errors = numpy.zeros((hs.size,2),dtype=complex) 
    wkb_errors = numpy.zeros((hs.size,4),dtype=complex)
    wkb_rerrors = numpy.zeros((hs.size,4),dtype=complex)
    sserrors = numpy.zeros((hs.size,4),dtype=complex)
    ss = numpy.zeros((hs.size,4),dtype=complex)
    
    solver = Solver(w,t=t,x=x,dx=dx,rtol=rtol,atol=atol)
    
    for i,h in enumerate(hs):
        # Take a RK and WKB step
        solver.h = h 
        x_rk, dx_rk, err_rk, ws, S0 = solver.RK_step()
        solver.ws = ws
        solver.S0 = S0
        solver.S0error = err_rk[2]
        x_wkb, dx_wkb, err_wkb, truncerr = solver.RKWKB_step()
        
        wkb_steps[i,:] = x_wkb, dx_wkb, x_wkb, dx_wkb
        wkb_errors[i,:] = truncerr[0], truncerr[1], err_wkb[0], err_wkb[1]
        wkb_rerrors[i,:] = (numpy.abs(truncerr[0])/numpy.abs(x_wkb),
        numpy.abs(truncerr[1])/numpy.abs(dx_wkb), numpy.abs(err_wkb[0])/numpy.abs(x_wkb),
        numpy.abs(err_wkb[1])/numpy.abs(dx_wkb))
        rk_steps[i,:] = x_rk, dx_rk
        rk_errors[i,:] = err_rk[:2]
        ss[i,:] = S0, solver.rkwkbsolver4.S1(h), solver.rkwkbsolver4.S2(h), solver.rkwkbsolver4.S3(h)
        sserrors[i,:] = solver.S0error, solver.rkwkbsolver4.S1(h), solver.rkwkbsolver4.Serror[2], solver.rkwkbsolver4.S3(h)
        #print('\n\n',h,'\n\n')
    
    fig, axes = plt.subplots(2,2, sharex=False)

    axes[0,0].set_title('Solution')
    axes[0,0].plot(t+hs,numpy.real(rk_steps[:,0]),'r.',alpha=0.5)
    axes[0,0].plot(t+hs,numpy.real(wkb_steps[:,0]),'g.',alpha=0.5)
    axes[0,0].plot(t+hs,sol(t+hs),'b')
    axes[0,0].set_ylim((-100.0, 100.0))

    axes[0,1].set_title('All absolute errors')
    axes[0,1].loglog(hs, numpy.abs(wkb_errors[:,0]),label='WKB truncerr x')
    axes[0,1].loglog(hs, numpy.abs(wkb_errors[:,1]),label='WKB truncerr dx')
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
    (axes[1,0].loglog(hs, numpy.max(numpy.abs(wkb_rerrors),axis=1),
    label='max WKBerror'))
    #(axes[1,0].loglog(hs, numpy.max(numpy.abs(wkb_rerrors[:,2:]),axis=1),
    #label='est WKBerror'))
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
