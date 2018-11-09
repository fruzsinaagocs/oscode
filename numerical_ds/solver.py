import numpy
import scipy.integrate
import sys

def integrate(f,a,b):
    return scipy.integrate.quad(lambda x: scipy.real(f(x)), a, b)[0] + 1j * scipy.integrate.quad(lambda x: scipy.imag(f(x)), a, b)[0] 


class Solver(object):
    def __init__(self,w,**kwargs):

        self.rkwkbsolver1 = RKWKBSolver1(w)
        self.rkwkbsolver2 = RKWKBSolver2(w)
        self.rkwkbsolver3 = RKWKBSolver3(w)
        self.rksolver = RKSolver(w)

        self.t = kwargs.pop('t', 0)
        self.x = kwargs.pop('x', 1)
        self.dx = kwargs.pop('dx', 0)
        self.err = kwargs.pop('err', 1e-4)

        self.h = 1

    def evolve(self, rk=False):
        while True:
            x_rk, dx_rk, err_rk, ws, S0 = self.RK_step()
            self.ws = ws
            self.S0 = S0
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
        #print(self.ws)
        self.rkwkbsolver1.ws = self.ws
        self.rkwkbsolver2.ws = self.ws       
        self.rkwkbsolver3.ws = self.ws
        self.rkwkbsolver1.S0 = self.S0
        self.rkwkbsolver2.S0 = self.S0       
        self.rkwkbsolver3.S0 = self.S0

        try:
            x1, dx1 = self.rkwkbsolver1.step(self.x,self.dx,self.t,self.h)
            x2, dx2 = self.rkwkbsolver2.step(self.x,self.dx,self.t,self.h)
            x3, dx3 = self.rkwkbsolver3.step(self.x,self.dx,self.t,self.h)
        except ZeroDivisionError:
            return numpy.inf, numpy.inf, numpy.inf
        #print(self.rkwkbsolver1.dw(self.t),self.rkwkbsolver1.Dw)
        #self.wkberror = abs(x3 - x2)
        return x3, dx3, abs(x3 - x2) #abs(x3*integrate(lambda t: 3/8 * self.rkwkbsolver3.dw(t)**2/self.rkwkbsolver3.w(t)**3 - 1/4 * self.rkwkbsolver3.ddw(t)/self.rkwkbsolver3.w(t)**2,self.t,self.t+self.h))

    def RK_step(self):
        x, dx, err, ws, S0 = self.rksolver.step(self.x,self.dx,self.t,self.h)
        return x, dx, err, ws, S0

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
        return numpy.array([y[1],-y[0] * self.w(t)**2, self.w(t)])

    def step(self,x0,dx0,t0,h):

        y0 = numpy.array([x0,dx0,0])
        k = []
        ws = []
        for c_s, a_s in zip(self.c, self.a):
            S = sum(a_si * k_i for a_si, k_i in zip(a_s, k))
            k_i = h * self.f(t0 + c_s * h, y0 + S)
            k.append(k_i)
            ws.append(k_i[-1]/h)

        y4 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b4,k)])
        y5 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b5,k)])
        w1 = ws[0]
        w2 = ws[1]
        w3 = ws[2]
        w4 = ws[5]
        w5 = ws[3]
        w6 = ws[4]
        ws = numpy.array([w1, w2, w3, w4, w5, w6])

        return y5[0], y5[1], abs(y5[0]-y4[0]), ws, y5[-1]
    
        
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
        #print(self.ws)
        self.Dw = self.d1w(h)
        self.Dw2 = self.d1w2(h)
        self.Dw3 = self.d1w3(h)
        self.Dw4 = self.d1w4(h)
        self.Dw5 = self.d1w5(h)
        self.DDw = self.d2w(h)
        self.DDDw = self.d3w(h)
        self.Dbw = self.d1wb(h)
        self.DDbw = self.d2wb(h)
        self.DDbw = self.d3wb(h) 
        #print(self.Dw, self.dw(t0))
        #print(self.DDw, self.ddw(t0))
        #print(self.DDDw, self.dddw(t0))
        
        ddx0 = -self.w(t0)**2 * x0
        t1 = t0 + h

        Ap, Am = self.A(t0,x0,dx0)
        Bp, Bm = self.B(t0,dx0,ddx0)

        x1 =  Ap * self.fp(t0,t1) + Am * self.fm(t0,t1)
        dx1 = Bp * self.dfpb(t1) * self.fp(t0,t1) + Bm * self.dfmb(t1) * self.fm(t0,t1) 

        return x1, dx1

    def d1w(self, h):
        d1w = (-43/4.0*self.ws[0] + 1536/35.0*self.ws[1] - 16384/285.0*self.ws[2] + 288/11.0*self.ws[3] - 371293/87780.0*self.ws[4] + 12/5.0*self.ws[5])/h
        return d1w

    def d1w2(self, h):
        d1w = (-35/96.0*self.ws[0] - 1136/105.0*self.ws[1] + 896/57.0*self.ws[2] - 105/22.0*self.ws[3] + 371293/702240.0*self.ws[4] - 7/24.0*self.ws[5])/h
        return d1w

    def d1w3(self, h):
        d1w = (95/768.0*self.ws[0] - 57/14.0*self.ws[1] - 72/95.0*self.ws[2] + 855/176.0*self.ws[3] - 371293/1123584.0*self.ws[4] + 57/320.0*self.ws[5])/h
        return d1w

    def d1w4(self, h):
        d1w = (-11/72.0*self.ws[0] + 352/105.0*self.ws[1] - 11264/855.0*self.ws[2] + 106/11.0*self.ws[3] + 371293/526680.0*self.ws[4] - 11/30.0*self.ws[5])/h
        return d1w

    def d1w5(self, h):
        d1w = (7315/26364.0*self.ws[0] -321024/76895.0*self.ws[1] + 1261568/125229.0*self.ws[2] - 191520/24167.0*self.ws[3] - 182663/29260.0*self.ws[4] + 17556/2197.0*self.ws[5])/h
        return d1w

    def d1wb(self, h):
        d1w = (-5/12.0*self.ws[0] + 128/21.0*self.ws[1] - 4096/285.0*self.ws[2] + 120/11.0*self.ws[3] - 371293/17556.0*self.ws[4] + 284/15.0*self.ws[5] )/h
        return d1w
    
    def d2w(self, h):
        d2w = (1553/18.0*self.ws[0] - 20736/35.0*self.ws[1] + 794624/855.0*self.ws[2] - 5040/11.0*self.ws[3] + 10767497/131670.0*self.ws[4] - 234/5.0*self.ws[5])/h**2
        return d2w
    
    def d2wb(self, h):
        d2w = (-269/18.0*self.ws[0] + 22528/105.0*self.ws[1] - 425984/855.0*self.ws[2] + 4064/11.0*self.ws[3] - 33045077/131670.0*self.ws[4] + 2702/15.0*self.ws[5] )/h**2
        return d2w
    
    def d3w(self, h):
        d3w = (- 1453/3.0*self.ws[0] + 21248/5.0*self.ws[1] - 2121728/285.0*self.ws[2] + 44304/11.0*self.ws[3] - 2599051/3135.0*self.ws[4] + 2404/5.0*self.ws[5])/h**3
        return d3w
    
    def d3wb(self, h):
        d3w = (-541/3.0*self.ws[0] + 85248/35.0*self.ws[1] - 1531904/285.0*self.ws[2] + 40464/11.0*self.ws[3] - 36015421/21945.0*self.ws[4] + 5412/5.0*self.ws[5] )/h**2
        return d3w

    #def S0(self, h):
    #    return h*self.ws[0] + h**2/2.0*self.Dw + h**3/6.0*self.DDw + h**4/24.0*self.DDDw

    def S1(self, h):
        return numpy.log(numpy.sqrt(self.ws[0]/self.ws[5])) 

    def S2(self, h):
        return  -1/4.0*(self.Dbw/self.ws[5]**2 - self.Dw/self.ws[0]**2) - 1/8.0*(16/135.0*self.Dw**2/self.ws[0]**3 + 6656/12825.0*self.Dw3**2/self.ws[2]**3 + 28561/56430.0*self.Dw4**2/self.ws[3]**3 - 9/50.0*self.Dw5**2/self.ws[4]**3 + 2/55.0*self.Dbw**2/self.ws[5]**3)*h   

class RKWKBSolver1(RKWKBSolver):
    def __init__(self,w):
        self.w = w

    def fp(self,t0,t1):
        ana = numpy.exp(1j * integrate(self.w,t0,t1))
        #print('S0',abs((ana -  numpy.exp(1j * self.S0))/ana)*100)
        return numpy.exp(1j * self.S0)

    def dfp(self,t):
        return 1j * self.w(t)
    
    def dfpb(self,t):
        return 1j * self.w(t)

    def ddfp(self,t):
        return -self.w(t)**2  + 1j * self.Dw 

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver2(RKWKBSolver1):
    def __init__(self,w):
        self.w = w

    def fp(self,t0,t1):
        ana = (self.w(t0)/self.w(t1))**(1/2)
        #print('S1',abs((ana - numpy.exp(self.S1(t1-t0)))/ana)*100)
        return numpy.exp(self.S1(t1-t0)) * super().fp(t0,t1)

    def dfp(self,t):
        return super().dfp(t) - self.Dw/self.w(t)/2

    def dfpb(self,t):
        return super().dfpb(t) - self.Dbw/self.w(t)/2

    def ddfp(self,t):
        return -self.w(t)**2 + 3/4 * (self.Dw/self.w(t))**2 - 1/2 * self.DDw/self.w(t) 

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver3(RKWKBSolver2):
    def __init__(self,w):
        self.w = w

    def fp(self,t0,t1):
        return super().fp(t0,t1) * numpy.exp(1j * self.S2(t1-t0))

    def dfp(self,t):
        return super().dfp(t) + 1j * ( 3/8 * self.Dw**2/self.w(t)**3 - 1/4  * self.DDw/self.w(t)**2)

    def dfpb(self,t):
        return super().dfpb(t) + 1j * (3/8 * self.Dbw**2/self.w(t)**3 - 1/4 * self.DDbw/self.w(t)**2)

    def ddfp(self,t):
        return -self.w(t)**2  - 1/4 * 1j * self.DDDw/self.w(t)**2 - 3/2 * 1j * self.Dw**3/self.w(t)**4 + 3/2 * 1j * self.Dw * self.DDw/self.w(t)**3 - 9/64 * self.Dw**4/self.w(t)**6 + 3/16 * self.Dw**2 * self.DDw / self.w(t)**5 - 1/16 * self.DDw**2/self.w(t)**4

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))
