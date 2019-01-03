#!/usr/bin/env python
import numpy

def d1w1(ws, h):
    weights = numpy.array([-15.0000000048537, 20.2828318761850,
    -8.07237453994912, 4.48936929577350, -2.69982662677053,
    0.999999999614819])
    return numpy.sum(weights*ws)/h
                                                                                
def d1w2(ws, h):
    weights = numpy.array([-3.57272991033049, 0.298532922755350e-7  ,
    5.04685352597795, -2.30565629452303, 1.30709499910514,
    -0.475562350082855]) 
    return numpy.sum(weights*ws)/h
                                                                                
def d1w3(ws, h):
    weights = numpy.array([0.969902096162109, -3.44251390568294,
    -0.781532641131861e-10, 3.50592393061265, -1.57271334190619,
    0.539401220892526]) 
    return numpy.sum(weights*ws)/h

def d1w4(ws, h):
    weights = numpy.array([-0.539401220892533, 1.57271334190621,
    -3.50592393061268, 0.782075077478921e-10, 3.44251390568290,
    -0.969902096162095])
    return numpy.sum(weights*ws)/h
                                                                                
def d1w5(ws, h):
    weights = numpy.array([0.475562350082834, -1.30709499910509,
    2.30565629452296, -5.04685352597787, -0.298533681980831e-7,
    3.57272991033053]) 
    return numpy.sum(weights*ws)/h
                                                                                
def d1w6(ws, h):
    weights = numpy.array([-0.999999999614890, 2.69982662677075,
    -4.48936929577383, 8.07237453994954, -20.2828318761854, 15.0000000048538]) 
    return numpy.sum(weights*ws)/h
                                                                                
def d2w1(ws, h):
    weights = numpy.array([140.000000016641, -263.163968874741,
    196.996471291466, -120.708905753218, 74.8764032980854, -27.9999999782328]) 
    return numpy.sum(weights*ws)/h**2  
                                                                                
def d2w6(ws, h):
    weights = numpy.array([-27.9999999782335, 74.8764032980873,
    -120.708905753221, 196.996471291469, -263.163968874744, 140.000000016642]) 
    return numpy.sum(weights*ws)/h**2 
                                                                                
def d3w1(ws, h):
    weights = numpy.array([-840.000000234078, 1798.12714381468,
    -1736.74461287884, 1322.01528240287, -879.397812956524, 335.999999851893]) 
    return numpy.sum(weights*ws)/h**3
                                                                                
def d3w6(ws, h):
    weights = numpy.array([-335.999999851897, 879.397812956534,
    -1322.01528240289, 1736.74461287886, -1798.12714381470, 840.000000234086]) 
    return numpy.sum(weights*ws)/h**3
                                                                                
def d4w1(ws, h):
    weights = numpy.array([3024.00000383582, -6923.06197480357,
    7684.77676018742, -6855.31809730784, 5085.60330881706, -2016.00000072890]) 
    return numpy.sum(weights*ws)/h**4   
                                                                                
def d1g1(gs, h):
    weights = numpy.array([-15.0000000048537, 20.2828318761850,
    -8.07237453994912, 4.48936929577350, -2.69982662677053,
    0.999999999614819])
    return numpy.sum(weights*gs)/h
                                                                                
def d1g6(gs, h):
    weights = numpy.array([-0.999999999614890, 2.69982662677075,
    -4.48936929577383, 8.07237453994954, -20.2828318761854, 15.0000000048538]) 
    return numpy.sum(weights*gs)/h
                                                                                
def d2g1(gs, h):
    weights = numpy.array([140.000000016641, -263.163968874741,
    196.996471291466, -120.708905753218, 74.8764032980854, -27.9999999782328]) 
    return numpy.sum(weights*gs)/h**2  
                                                                                
def d2g6(gs, h):
    weights = numpy.array([-27.9999999782335, 74.8764032980873,
    -120.708905753221, 196.996471291469, -263.163968874744, 140.000000016642]) 
    return numpy.sum(weights*gs)/h**2 
                                                                                
def d3g1(gs, h):
    weights = numpy.array([-840.000000234078, 1798.12714381468,
    -1736.74461287884, 1322.01528240287, -879.397812956524, 335.999999851893]) 
    return numpy.sum(weights*gs)/h**3

def d1w2_5(ws5, h):
    weights = numpy.array([-2.48198050935042, 0.560400997591235e-8,
    3.49148624058567, -1.52752523062733, 0.518019493788063])
    return numpy.sum(weights*ws5)/h
                                                                                
def d1w3_5(ws5, h):
    weights = numpy.array([0.750000000213852, -2.67316915534181,
    0.360673032443906e-10, 2.67316915534181, -0.750000000213853])
    return numpy.sum(weights*ws5)/h
                                                                                
def d1w4_5(ws5, h):
    weights = numpy.array([-0.518019493788065, 1.52752523062733,
    -3.49148624058568, -0.560400043118500e-8, 2.48198050935041])
    return numpy.sum(weights*ws5)/h

def S0(ws,ws5,t0,t1):
    return integrate(ws,ws5,t1-t0)

def S1(ws,gs,ws5,gs5,h):
    integral, error = integrate(gs,gs5,h)
    Serror[1] = error
    return (numpy.log(numpy.sqrt(ws[0]/ws[5])) -
    integral) 

def S2(ws,gs,ws5,gs5,Dws,Dws5,h):
    integrands6 = (Dws**2/ws**3 + 4*Dws*gs/ws**2 +
    4*gs**2/ws)
    integrands5 = (Dws5**2/ws5**3 +
    4*Dws5*gs5/ws5**2 + 4*gs5**2/ws5)
    integral, error = integrate(integrands6,integrands5,h)
    Serror[2] = -1/8.0*1j*error
    return (-1/4.0*(Dws[5]/ws[5]**2  + 2*gs[5]/ws[5]-
    Dws[0]/ws[0]**2 - 2*gs[0]/ws[0])
    -1/8.0*integral)

def S3(ws,gs,Dws,Dgs,DDws,h):
    S3 = (1/4.0*(gs[5]**2/ws[5]**2 - gs[0]**2/ws[0]**2) +
    1/4.0*(Dgs[-1]/ws[5]**2 -
    Dgs[0]/ws[0]**2)-3/16.0*(Dws[5]**2/ws[5]**4 -
    Dws[0]**2/ws[0]**4) + 1/8.0*(DDws[-1]/ws[5]**3 -
    DDws[0]/ws[0]**3))
    return S3

def integrate(integrand6,integrand5,h):
    legws6=(
    numpy.array([1.0/15.0,0.3784749562978469803166,0.5548583770354863530167,0.5548583770354863530167,0.3784749562978469803166,1.0/15.0]))
    legws5 =(
    numpy.array([1.0/10.0,49.0/90.0,32.0/45.0,49.0/90.0,1.0/10.0]))

    x6 = h/2.0*numpy.sum(legws6*integrand6)
    x5 = h/2.0*numpy.sum(legws5*integrand5)
    return x6, x6-x5

