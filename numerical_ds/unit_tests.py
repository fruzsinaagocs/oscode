#!/usr/bin/env python
import test
import nose
import numpy
import random

def test_dws():
    
    n = 80
    t = random.random()*100
    d1wa = numpy.sqrt(n**2-1) * -2 * t / (1 + t**2)**2
    d2wa = numpy.sqrt(n**2-1) * (6 * t**2 - 2 ) / (1 + t**2)**3
    d3wa = numpy.sqrt(n**2-1) * 24 * t * (1 - t**2) / (1 + t**2)**4

    print(n, t, test.dw(t), d1wa)
    print(n, t, test.ddw(t), d2wa) 
    print(n, t, test.dddw(t), d3wa)
    assert abs((test.dw(t) - d1wa)/d1wa) < 1e-3 
    assert abs((test.ddw(t) - d2wa)/d2wa) < 1e-3 
    assert abs((test.dddw(t) - d3wa)/d3wa) < 1e-3 

def test_dws_back():

    n = 80
    t = random.random()*100
    d1wa = numpy.sqrt(n**2-1) * -2 * t / (1 + t**2)**2
    d2wa = numpy.sqrt(n**2-1) * (6 * t**2 - 2 ) / (1 + t**2)**3
    d3wa = numpy.sqrt(n**2-1) * 24 * t * (1 - t**2) / (1 + t**2)**4

    print(n, t, test.dwb(t), d1wa)
    print(n, t, test.ddwb(t), d2wa) 
    print(n, t, test.dddwb(t), d3wa)
    assert abs((test.dwb(t) - d1wa)/d1wa) < 1e-3 
    assert abs((test.ddwb(t) - d2wa)/d2wa) < 1e-3 
    assert abs((test.dddwb(t) - d3wa)/d3wa) < 1e-3 


