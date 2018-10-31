import numpy
import scipy.integrate
import sys

def integrate(f,a,b):
    return scipy.integrate.quad(lambda x: scipy.real(f(x)), a, b)[0] + 1j * scipy.integrate.quad(lambda x: scipy.imag(f(x)), a, b)[0] 


class Solver(object):
    def __init__(self,w,dw,ddw,dddw,**kwargs):

        self.rkwkbsolver1 = RKWKBSolver1(w,dw)
        self.rkwkbsolver2 = RKWKBSolver2(w,dw,ddw)
        self.rkwkbsolver3 = RKWKBSolver3(w,dw,ddw,dddw)
        self.rksolver = RKSolver(w)

        self.t = kwargs.pop('t', 0)
        self.x = kwargs.pop('x', 1)
        self.dx = kwargs.pop('dx', 0)
        self.err = kwargs.pop('err', 1e-4)

        self.h = 1

    def evolve(self, rk=False):
        while True:
            x_rk, dx_rk, err_rk = self.RK_step()
            if not rk:
                x_wkb, dx_wkb, err_wkb = self.RKWKB_step()
                wkb = err_wkb < err_rk

            if rk:
                wkb = False

            if wkb:
                x = x_wkb
                dx = dx_wkb
                err = err_wkb
                #err = self.wkberror
            else:
                x = x_rk
                dx = dx_rk
                err = err_rk


            if err < self.err:
                self.t += self.h
                self.x = x
                self.dx = dx
                if wkb:
                    self.h *= 0.95 * (self.err/err)**(1/5.0)
                else:
                    self.h *= 0.95 * (self.err/err)**(1/5.0)
                yield {'t':self.t, 'x':self.x, 'dx':self.dx, 'h':self.h, 'err':err, 'wkb':wkb}
            else:
                if wkb:
                    self.h *= 0.95 * (self.err/err)**(1/4.0)
                else:
                    self.h *= 0.95 * (self.err/err)**(1/4.0)

    def RKWKB_step(self):
        try:
            x1, dx1 = self.rkwkbsolver1.step(self.x,self.dx,self.t,self.h)
            x2, dx2 = self.rkwkbsolver2.step(self.x,self.dx,self.t,self.h)
            x3, dx3 = self.rkwkbsolver3.step(self.x,self.dx,self.t,self.h)
        except ZeroDivisionError:
            return numpy.inf, numpy.inf, numpy.inf
        #self.wkberror = abs(x3 - x2)
        return x3, dx3, abs(x3 - x2) #abs(x3*integrate(lambda t: 3/8 * self.rkwkbsolver3.dw(t)**2/self.rkwkbsolver3.w(t)**3 - 1/4 * self.rkwkbsolver3.ddw(t)/self.rkwkbsolver3.w(t)**2,self.t,self.t+self.h))

    def RK_step(self):
        x, dx, err = self.rksolver.step(self.x,self.dx,self.t,self.h)
        return x, dx, err

class RKSolver(object):
    def __init__(self,w):
        self.w = w
        self.c = [0, 1/4, 3/8, 12/13, 1, 1/2]
        self.b5 = numpy.array([16/135, 0, 6656/12825, 28561/56430, -9/50,  2/55])
        self.b4 = numpy.array([25/216, 0, 1408/2565,  2197/4104,   -1/5,   0   ])
         
        self.a = [
                [],
                [1/4],
                [3/32,    9/32],
                [1932/2197,   -7200/2197,  7296/2197],
                [439/216, -8,  3680/513,    -845/4104],
                [-8/27,   2,   -3544/2565,  1859/4104,   -11/40]
                ]

    def f(self,t,y):
        return numpy.array([y[1],-y[0] * self.w(t)**2])

    def step(self,x0,dx0,t0,h):

        y0 = numpy.array([x0,dx0])
        k = []
        for c_s, a_s in zip(self.c, self.a):
            S = sum(a_si * k_i for a_si, k_i in zip(a_s, k))
            k_i = h * self.f(t0 + c_s * h, y0 + S)
            k.append(k_i)

        y4 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b4,k)])
        y5 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b5,k)])

        return y5[0], y5[1], abs(y5[0]-y4[0])
    
        
class RKWKBSolver(object):
    def A(self,t0,x0,dx0):
        Ap = (dx0 - x0 * self.dfm(t0)) / (self.dfp(t0) - self.dfm(t0))
        Am = (dx0 - x0 * self.dfp(t0)) / (self.dfm(t0) - self.dfp(t0))
        return Ap, Am

    def B(self,t0,dx0,ddx0):
        Bp = (ddx0 * self.dfm(t0) - dx0 * self.ddfm(t0)) / (self.ddfp(t0) * self.dfm(t0) - self.ddfm(t0) * self.dfp(t0))
        Bm = (ddx0 * self.dfp(t0) - dx0 * self.ddfp(t0)) / (self.ddfm(t0) * self.dfp(t0) - self.ddfp(t0) * self.dfm(t0))
        return Bp, Bm

    def step(self, x0, dx0, t0, h):
        ddx0 = -self.w(t0)**2 * x0
        t1 = t0 + h

        Ap, Am = self.A(t0,x0,dx0)
        Bp, Bm = self.B(t0,dx0,ddx0)

        x1 =  Ap * self.fp(t0,t1) + Am * self.fm(t0,t1)
        dx1 = Bp * self.dfp(t1) * self.fp(t0,t1) + Bm * self.dfm(t1) * self.fm(t0,t1) 

        return x1, dx1

class RKWKBSolver1(RKWKBSolver):
    def __init__(self,w,dw):
        self.w = w
        self.dw = dw

    def fp(self,t0,t1):
        return numpy.exp(1j * integrate(self.w,t0,t1))

    def dfp(self,t):
        return 1j * self.w(t) 

    def ddfp(self,t):
        return -self.w(t)**2  + 1j * self.dw(t) 

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver2(RKWKBSolver1):
    def __init__(self,w,dw,ddw):
        self.w = w
        self.dw = dw
        self.ddw = ddw

    def fp(self,t0,t1):
        return (self.w(t0)/self.w(t1))**(1/2) * super().fp(t0,t1)

    def dfp(self,t):
        return super().dfp(t) - self.dw(t)/self.w(t)/2

    def ddfp(self,t):
        return -self.w(t)**2 + 3/4 * (self.dw(t)/self.w(t))**2 - 1/2 * self.ddw(t)/self.w(t) 

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver3(RKWKBSolver2):
    def __init__(self,w,dw,ddw,dddw):
        self.w = w
        self.dw = dw
        self.ddw = ddw
        self.dddw = dddw

    def fp(self,t0,t1):
        return super().fp(t0,t1) * numpy.exp(1j * integrate(
            lambda t: 3/8 * self.dw(t)**2/self.w(t)**3 - 1/4 * self.ddw(t)/self.w(t)**2,t0,t1))

    def dfp(self,t):
        return super().dfp(t) + 1j * ( 3/8 * self.dw(t)**2/self.w(t)**3 - 1/4  * self.ddw(t)/self.w(t)**2)

    def ddfp(self,t):
        return -self.w(t)**2  - 1/4 * 1j * self.dddw(t)/self.w(t)**2 - 3/2 * 1j * self.dw(t)**3/self.w(t)**4 + 3/2 * 1j * self.dw(t) * self.ddw(t)/self.w(t)**3 - 9/64 * self.dw(t)**4/self.w(t)**6 + 3/16 * self.dw(t)**2 * self.ddw(t) / self.w(t)**5 - 1/16 * self.ddw(t)**2/self.w(t)**4

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))
