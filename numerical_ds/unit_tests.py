#!/usr/bin/env python
import test
import nose
import numpy
import random

def assert_eq(obt, exp):
    assert obt < 10*exp, "\n Expected:\n%s\n*Obtained:\n%s" % (exp, obt)

def test_derivatives():
 
    stencils5 =(
    numpy.array([-1.0,-0.65465367070797714380,0.0,0.65465367070797714380,1.0]))
    stencils6 =(
    numpy.array([-1.0,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.0]))

    t = random.random()*100 + 1e5
    h = random.random()
    w = lambda t: t**(1/2.0)
    g = lambda t: t**(1/2.0)
    d1wa = lambda t: 1/2.0*t**(-1/2.0)
    d2wa = lambda t: -1/4.0*t**(-3/2.0)
    d3wa = lambda t: 3/8.0*t**(-5/2.0)
    d4wa = lambda t: -15/16.0*t**(-7/2.0)
    d1ga = d1wa
    d2ga = d2wa
    d3ga = d3wa
    ts6 = h/2*(1+stencils6) + t
    ts5 = h/2*(1+stencils5) + t
    ws = w(ts6)
    gs = g(ts6)
    ws5 = w(ts5)

    print(t, h)
    print(test.d4w1(ws, h), d4wa(ts6[0]))

    assert_eq(abs((test.d1w1(ws,h) - d1wa(ts6[0]))/d1wa(ts6[0])), 1e-4)
    assert_eq(abs((test.d1w2(ws,h) - d1wa(ts6[1]))/d1wa(ts6[1])), 1e-4)
    assert_eq(abs((test.d1w3(ws,h) - d1wa(ts6[2]))/d1wa(ts6[2])) , 1e-4)
    assert_eq(abs((test.d1w4(ws,h) - d1wa(ts6[3]))/d1wa(ts6[3])) , 1e-4)
    assert_eq(abs((test.d1w5(ws,h) - d1wa(ts6[4]))/d1wa(ts6[4])) , 1e-4)
    assert_eq(abs((test.d1w6(ws,h) - d1wa(ts6[5]))/d1wa(ts6[5])) , 1e-4)
    assert_eq(abs((test.d2w1(ws,h) - d2wa(ts6[0]))/d2wa(ts6[0])) , 1e-4)
    assert_eq(abs((test.d2w6(ws,h) - d2wa(ts6[5]))/d2wa(ts6[5])) , 1e-4)
    assert_eq(abs((test.d3w1(ws,h) - d3wa(ts6[0]))/d3wa(ts6[0])) , 1e-4)
    assert_eq(abs((test.d3w6(ws,h) - d3wa(ts6[5]))/d3wa(ts6[5])) , 1e-4)
#    assert_eq(abs((test.d4w1(ws,h) - d4wa(ts6[0]))/d4wa(ts6[0])) , 1e-4)
    assert_eq(abs((test.d1g1(gs,h) - d1ga(ts6[0]))/d1ga(ts6[0])) , 1e-4)
    assert_eq(abs((test.d1g6(gs,h) - d1ga(ts6[5]))/d1ga(ts6[5])) , 1e-4)
    assert_eq(abs((test.d2g1(gs,h) - d2ga(ts6[0]))/d2ga(ts6[0])) , 1e-4)
    assert_eq(abs((test.d2g6(gs,h) - d2ga(ts6[5]))/d2ga(ts6[5])) , 1e-4)
    assert_eq(abs((test.d3g1(gs,h) - d3ga(ts6[0]))/d3ga(ts6[0])) , 1e-4)
    assert_eq(abs((test.d1w2_5(ws5,h) - d1wa(ts5[1]))/d1wa(ts5[1])) , 1e-4)
    assert_eq(abs((test.d1w3_5(ws5,h) - d1wa(ts5[2]))/d1wa(ts5[2])) , 1e-4)
    assert_eq(abs((test.d1w4_5(ws5,h) - d1wa(ts5[3]))/d1wa(ts5[3])) , 1e-4)

        #print('S1: ', solver.rkwkbsolver4.S1(h))
        #'S3: ', solver.rkwkbsolver4.S3(h))
        #print('analytic S1: ', -1/4*numpy.log((t+h)/t) -
        #2/3*((t+h)**(3/2)-t**(3/2)))
        #print('error on S1: ', solver.rkwkbsolver4.Serror[1])
        #print('S2: ', solver.rkwkbsolver4.S2(h))
        #print('analytic S2: ', -1j*5/48*((t+h)**(-3/2)-t**(-3/2)) +
        #1/2*1j*((t+h)**(-1/4)-t**(-1/4)) - 1/2*1j*h)
        #print('S3: ', solver.rkwkbsolver4.S3(h))
        #print('analytic S3: ', -5/64*((t+h)**(-3) - t**(-3)) +
        #1/16*((t+h)**(-7/4) - t**(-7/4)) + 1/4*((t+h)**(-1/2)-t**(-1/2)))
        #print('g(t): ', solver.rkwkbsolver2.gs[0], numpy.sqrt(t))
        #print('dfp order 1:', solver.rkwkbsolver2.dfp(h))
        
        #print('analytic dfp order 1:', 1j*t**(1/2) -1/4*t**(-1) -t**(1/2))
        
        #print('dfp: ', solver.rkwkbsolver4.dfp(h))
        #print('analytic dfp: ', 1/2*1j*t**(1/2) -
        #1/4*(t**(-1)) - (t**(1/2)) +
        #5/32*1j*(t**(-5/2)) - 1/4*1j*(t**(-1)) +
        #15/64*(t**(-4)) - 3/16*(t**(-5/2)))
        #print('ddfp order 2: ', solver.rkwkbsolver3.ddfp(h))
        #print('analytic: ', -1j*t + 3/32*t**(-2) + 1/4*t**(-1/2) +
        #1/2*1j*t**(-1/2) - 15/32*1j*t**(-7/2) + 1/16*1j*t**(-2) -
        #25/1024*t**(-5) + 5/64*t**(-7/2) + 3/4*t)
        #print('ddfp: ', solver.rkwkbsolver4.ddfp(h))
        #print('analytic ddfp: ', 225/4096*(t**(-8)) -
        #45/512*(t**(-13/2)) - 1069/1024*(t**(-5))
        #+ 11/64*(t**(-7/2)) + 15/32*(t**(-2)) +
        #1/4*(t**(-1/2)) + 1j*1/2*(t**(-1/2)) +
        #75/1024*1j*(t**(-13/2)) - 1/8*1j*(
        #t**(-2)) - 1j*t - 45/256*1j*(t**(-5)) -
        #9/64*1j*(t**(-7/2)) + 3/4*t)



