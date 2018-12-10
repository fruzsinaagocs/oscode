import numpy
import sys

class Solver(object):
    def __init__(self,w,**kwargs):

        self.rkwkbsolver1 = RKWKBSolver1(w)
        self.rkwkbsolver2 = RKWKBSolver2(w)
        self.rkwkbsolver3 = RKWKBSolver3(w)
        self.rkwkbsolver4 = RKWKBSolver4(w)
        self.rksolver = RKSolver(w)

        self.t = kwargs.pop('t', 0)
        self.x = kwargs.pop('x', 1)
        self.dx = kwargs.pop('dx', 0)
        self.rtol = kwargs.pop('rtol', 1e-4)
        self.atol = kwargs.pop('atol', 0.0)

        self.h = 1.0

    def evolve(self, rk=False):
        while True:
            # Take RKF and WKB steps
            x_rk, dx_rk, err_rk, ws, ws5 = self.RK_step()
            self.ws = ws
            self.ws5 = ws5
                
            x_wkb, dx_wkb, err_wkb, truncerr = self.RKWKB_step()
            
            # Construct error estimates on steps
            deltas_wkb = (numpy.array([numpy.abs(truncerr[0])/numpy.abs(x_wkb),
            numpy.abs(truncerr[1])/numpy.abs(dx_wkb),
            numpy.abs(err_wkb[0])/numpy.abs(x_wkb),
            numpy.abs(err_wkb[1])/numpy.abs(dx_wkb)]))
            delta_rk = (numpy.max([1e-10,numpy.max(numpy.array([numpy.abs(err_rk[0])/numpy.abs(x_rk),
            numpy.abs(err_rk[1])/numpy.abs(dx_rk)]))]))
            delta_wkb = numpy.max([1e-10,numpy.max(deltas_wkb)])
            maxplace = numpy.argmax(deltas_wkb)
            
            # Predict next stepsize for each 
            h_rk = self.h*(self.rtol/delta_rk)**(1/5.0)
            #d_wkb = numpy.max(deltas_wkb[:2])
            #h_wkb = self.h*(self.rtol/d_wkb)**(1/4.0) 
            if maxplace <= 1:
                h_wkb = self.h*(self.rtol/delta_wkb)**(1/1.0)
            else:
                h_wkb = self.h*(self.rtol/delta_wkb)**(1/8.0)
            # Choose the one with larger predicted stepsize
            wkb = h_wkb > h_rk

            if wkb:
                x = x_wkb
                dx = dx_wkb
                
                # To have symmetric stepsizes but not quite symmetric switching,
                # comment these
                #if maxplace <= 1:
                delta_wkb = numpy.max(deltas_wkb[2:])
                h_next = self.h*(self.rtol/delta_wkb)**(1/8.0)
                #else:
                #    h_next = h_wkb
                
                # To have slightly asymm (<~3 steps) switching but large and
                # symmetric stepsizes in WKB, comment the next line.
                #h_next = h_wkb
                
                err = delta_wkb
                errortypes=["truncation", "S integral"]
                #print("{} error dominates".format(errortypes[maxplace//2]))
                #print("S errors: ", self.rkwkbsolver4.Serror)

            else:
                x = x_rk
                dx = dx_rk
                h_next = h_rk
                err = delta_rk

            # Check if chosen step is successful
            if h_next >= self.h:
                self.t += self.h
                self.x = x
                self.dx = dx
                h_prev = self.h
                self.h = h_next
                yield {'t':self.t, 'x':self.x, 'dx':self.dx, 'h':h_prev, 'err':err, 'wkb':wkb}
            else:
                if wkb:
                    if maxplace <=1:
                        self.h *= 0.95*(self.rtol/delta_wkb)**(1/1.0)
                    else:
                        self.h *= (self.rtol/delta_wkb)**(1/7.0)
                else:
                    self.h *= (self.rtol/delta_rk)**(1/4.0)

    def RKWKB_step(self):
        self.rkwkbsolver1.ws = self.ws
        self.rkwkbsolver2.ws = self.ws       
        self.rkwkbsolver3.ws = self.ws
        self.rkwkbsolver4.ws = self.ws
        self.rkwkbsolver1.ws5 = self.ws5
        self.rkwkbsolver2.ws5 = self.ws5       
        self.rkwkbsolver3.ws5 = self.ws5
        self.rkwkbsolver4.ws5 = self.ws5

        try:
            x1, dx1, x1_err, dx1_err = self.rkwkbsolver1.step(self.x,self.dx,self.t,self.h)
            x2, dx2, x2_err, dx2_err = self.rkwkbsolver2.step(self.x,self.dx,self.t,self.h)
            x3, dx3, x3_err, dx3_err = self.rkwkbsolver3.step(self.x,self.dx,self.t,self.h)
            x4, dx4, x4_err, dx4_err = self.rkwkbsolver4.step(self.x,self.dx,self.t,self.h)
        except ZeroDivisionError:
            return numpy.inf, numpy.inf, numpy.inf
        return x4, dx4, numpy.array([x4_err, dx4_err]), numpy.array([x4-x3,dx4-dx3])

    def RK_step(self):
        x, dx, err, ws, ws5 = self.rksolver.step(self.x,self.dx,self.t,self.h)
        return x, dx, err, ws, ws5

class RKSolver(object):
    def __init__(self,w):
        self.w = w
        
        # 6-stage, 5th order step:
        self.c5 = [0,0.117472338035267653574,0.357384241759677451843,0.642615758240322548157,0.882527661964732346426,1] 
        self.b5 = numpy.array([0.1127557227351729739820,0,0.5065579732655351768382,0.04830040376995117617928,0.3784749562978469803166,-0.04608905606850630731611])
        self.a5 = [
                [],
                [0.1174723380352676535740],
                [-0.1862479800651504276304,0.5436322218248278794734],
                [-0.6064303885508280518989,1,0.2490461467911506000559],
                [2.899356540015731406420,-4.368525611566240669139,2.133806714786316899991,0.2178900187289247091542],
                [18.67996349995727204273,-28.85057783973131956546,10.72053408420926869789,1.414741756508049078612,-0.9646615009432702537787],
                ]

        # 4-stage, 4th order step:
        self.c4 = [0,0.172673164646011428100,0.827326835353988571900,1]
        self.b4 = numpy.array([-0.08333333333333333333558,0.5833333333333333333357,0.5833333333333333333356,-0.08333333333333333333558])
        self.a4 = [
                [],
                [0.172673164646011428100],
                [-1.568317088384971429762,2.395643923738960001662],
                [-8.769507466172720011410,10.97821961869480000808,-1.208712152522079996671]
                ]

    def f(self,t,y):
        w = self.w(t)
        return numpy.array([y[1],-y[0]*w**2, w])

    def step(self,x0,dx0,t0,h):

        y0 = numpy.array([x0,dx0,0])
        ws, k5 = numpy.zeros(6,dtype=complex), numpy.zeros((6,3),dtype=complex)
        ws5, k4 = numpy.zeros(4,dtype=complex), numpy.zeros((4,3),dtype=complex)
        # 6-stage, 5th order step
        for i, c_s, a_s in zip(range(6),self.c5, self.a5):
            S = sum(a_si * k_i for a_si, k_i in zip(a_s, k5))
            k_i = h * self.f(t0 + c_s * h, y0 + S)
            k5[i,:] = k_i
            ws[i] = k_i[-1]/h
        # 4-stage, 4th order step
        k4[0] = k5[0]
        ws5[0] = ws[0]
        for i, c_s, a_s in zip(range(1,4),self.c4[1:], self.a4[1:]):
            S = sum(a_si * k_i for a_si, k_i in zip(a_s, k4))
            k_i = h * self.f(t0 + c_s * h, y0 + S)
            k4[i,:] = k_i
            ws5[i] = k_i[-1]/h


        y4 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b4,k4)])
        y5 = y0 + sum([b_i * k_i for b_i, k_i in zip(self.b5,k5)])
        delta = sum([b_i * k_i for b_i, k_i in zip(self.b5,k5)])-sum([b_i * k_i for b_i, k_i in zip(self.b4,k4)])
        #print(ws5)
        ws5 = numpy.insert(ws5,2,self.w(t0+h/2))

        return y5[0], y5[1], delta, ws, ws5
        
class RKWKBSolver(object):
    def __init__(self,w):
        self.legxs6 =(
        numpy.array([-1.0,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.0]))
        self.legws6 =(
        numpy.array([1.0/15.0,0.3784749562978469803166,0.5548583770354863530167,0.5548583770354863530167,0.3784749562978469803166,1.0/15.0]))
        self.legxs5 =(
        numpy.array([-1.0,-0.65465367070797714380,0.0,0.65465367070797714380,1.0]))
        self.legws5 =(
        numpy.array([1.0/10.0,49.0/90.0,32.0/45.0,49.0/90.0,1.0/10.0]))
        self.w = w
    
    def A(self,t0,x0,dx0):
        Ap = (dx0 - x0 * self.dfm(t0)) / (self.dfp(t0) - self.dfm(t0))
        Am = (dx0 - x0 * self.dfp(t0)) / (self.dfm(t0) - self.dfp(t0))
        return Ap, Am

    def B(self,t0,dx0,ddx0):
        Bp = (ddx0 * self.dfm(t0) - dx0 * self.ddfm(t0)) / (self.ddfp(t0) * self.dfm(t0) - self.ddfm(t0) * self.dfp(t0))
        Bm = (ddx0 * self.dfp(t0) - dx0 * self.ddfp(t0)) / (self.ddfm(t0) * self.dfp(t0) - self.ddfp(t0) * self.dfm(t0))
        return Bp, Bm

    def step(self, x0, dx0, t0, h):
        self.Dws = numpy.array([self.d1w1(h), self.d1w2(h), self.d1w3(h), self.d1w4(h), self.d1w5(h), self.d1w6(h)])
        self.Dws5 = numpy.array([self.Dws[0], self.d1w2_5(h), self.d1w3_5(h), self.d1w4_5(h), self.Dws[5]])
        self.DDws = numpy.array([self.d2w1(h), self.d2w6(h)])
        self.DDDws = numpy.array([self.d3w1(h), self.d3w6(h)])
        self.DDDDws = numpy.array([self.d4w1(h)])
        self.Serror = numpy.zeros(4, dtype=complex)
        self.Serror[0] = 1j*self.S0(t0,t0+h)[1]
        
        ddx0 = -self.ws[0]**2 * x0
        t1 = t0 + h

        Ap, Am = self.A(t0,x0,dx0)
        Bp, Bm = self.B(t0,dx0,ddx0)

        x1 =  Ap * self.fp(t0,t1) + Am * self.fm(t0,t1)
        dx1 = Bp * self.dfpb(t1) * self.fp(t0,t1) + Bm * self.dfmb(t1) * self.fm(t0,t1) 

        # Error estimate on answer based on error on S_i integrals
        error_fp = numpy.sum(numpy.abs(self.Serror))*self.fp(t0,t1)
        error_fm = numpy.conj(numpy.sum(numpy.abs(self.Serror)))*self.fm(t0,t1)
        error_dfp = self.dfpb(t1)/self.fp(t0,t1)*error_fp
        error_dfm = self.dfmb(t1)/self.fm(t0,t1)*error_fm
        error_x = Ap*error_fp + Am*error_fm
        error_dx = Bp*error_dfp + Bm*error_dfm

        return x1, dx1, error_x, error_dx

# Derivatives at Gauss-Lobatto (n=6) abscissas

    def d1w1(self, h):
        weights = numpy.array([-15.0000000048537, 20.2828318761850,
        -8.07237453994912, 4.48936929577350, -2.69982662677053,
        0.999999999614819])
        return numpy.sum(weights*self.ws)/h

    def d1w2(self, h):
        weights = numpy.array([-3.57272991033049, 0.298532922755350e-7  ,
        5.04685352597795, -2.30565629452303, 1.30709499910514,
        -0.475562350082855]) 
        return numpy.sum(weights*self.ws)/h

    def d1w3(self, h):
        weights = numpy.array([0.969902096162109, -3.44251390568294,
        -0.781532641131861e-10, 3.50592393061265, -1.57271334190619,
        0.539401220892526]) 
        return numpy.sum(weights*self.ws)/h
    
    def d1w4(self, h):
        weights = numpy.array([-0.539401220892533, 1.57271334190621,
        -3.50592393061268, 0.782075077478921e-10, 3.44251390568290,
        -0.969902096162095])
        return numpy.sum(weights*self.ws)/h

    def d1w5(self, h):
        weights = numpy.array([0.475562350082834, -1.30709499910509,
        2.30565629452296, -5.04685352597787, -0.298533681980831e-7,
        3.57272991033053]) 
        return numpy.sum(weights*self.ws)/h

    def d1w6(self, h):
        weights = numpy.array([-0.999999999614890, 2.69982662677075,
        -4.48936929577383, 8.07237453994954, -20.2828318761854, 15.0000000048538]) 
        return numpy.sum(weights*self.ws)/h

    def d2w1(self, h):
        weights = numpy.array([140.000000016641, -263.163968874741,
        196.996471291466, -120.708905753218, 74.8764032980854, -27.9999999782328]) 
        return numpy.sum(weights*self.ws)/h**2  

    def d2w6(self, h):
        weights = numpy.array([-27.9999999782335, 74.8764032980873,
        -120.708905753221, 196.996471291469, -263.163968874744, 140.000000016642]) 
        return numpy.sum(weights*self.ws)/h**2 

    def d3w1(self, h):
        weights = numpy.array([-840.000000234078, 1798.12714381468,
        -1736.74461287884, 1322.01528240287, -879.397812956524, 335.999999851893]) 
        return numpy.sum(weights*self.ws)/h**3

    def d3w6(self, h):
        weights = numpy.array([-335.999999851897, 879.397812956534,
        -1322.01528240289, 1736.74461287886, -1798.12714381470, 840.000000234086]) 
        return numpy.sum(weights*self.ws)/h**3

    def d4w1(self, h):
        weights = numpy.array([3024.00000383582, -6923.06197480357,
        7684.77676018742, -6855.31809730784, 5085.60330881706, -2016.00000072890]) 
        return numpy.sum(weights*self.ws)/h**4   

# First derivatives at Gauss-Lobatto (n=3) abscissas, excluding endpoints
    
    def d1w2_5(self,h):
        weights = numpy.array([-2.48198050935042, 0.560400997591235e-8,
        3.49148624058567, -1.52752523062733, 0.518019493788063])
        return numpy.sum(weights*self.ws5)/h

    def d1w3_5(self,h):
        weights = numpy.array([0.750000000213852, -2.67316915534181,
        0.360673032443906e-10, 2.67316915534181, -0.750000000213853])
        return numpy.sum(weights*self.ws5)/h

    def d1w4_5(self,h):
        weights = numpy.array([-0.518019493788065, 1.52752523062733,
        -3.49148624058568, -0.560400043118500e-8, 2.48198050935041])
        return numpy.sum(weights*self.ws5)/h

##########################################################################

    def S0(self,t0,t1):
        return self.integrate(self.ws,self.ws5,t1-t0)

    def S1(self, h):
        self.Serror[1] = 0.0
        return numpy.log(numpy.sqrt(self.ws[0]/self.ws[5])) 

    def S2(self, h):
        integrands6 = self.Dws**2/self.ws**3 
        integrands5 = self.Dws5**2/self.ws5**3
        integral, error = self.integrate(integrands6,integrands5,h)
        self.Serror[2] = -1/8.0*1j*error
        return -1/4.0*(self.Dws[5]/self.ws[5]**2 - self.Dws[0]/self.ws[0]**2) -1/8.0*integral

    def S3(self, h):
        S3 = -3/16.0*(self.Dws[5]**2/self.ws[5]**4 - self.Dws[0]**2/self.ws[0]**4) + 1/8.0*(self.DDws[-1]/self.ws[5]**3 - self.DDws[0]/self.ws[0]**3)
        #print(S3)
        #self.Serror[3] = S3*0.1
        return S3

    def integrate(self,integrand6,integrand5,h):
        x6 = h/2.0*numpy.sum(self.legws6*integrand6)
        x5 = h/2.0*numpy.sum(self.legws5*integrand5)
        return x6, x6-x5

class RKWKBSolver1(RKWKBSolver):

    def fp(self,t0,t1):
        return numpy.exp(1j * self.S0(t0,t1)[0])

    def dfp(self,t):
        return 1j * self.ws[0]
    
    def dfpb(self,t):
        return 1j * self.ws[-1]

    def ddfp(self,t):
        return -self.ws[0]**2  + 1j * self.Dws[0]

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver2(RKWKBSolver1):

    def fp(self,t0,t1):
        ana = (self.ws[-1]/self.ws[0])**(1/2)
        return numpy.exp(self.S1(t1-t0)) * super().fp(t0,t1)

    def dfp(self,t):
        return super().dfp(t) - self.Dws[0]/self.ws[0]/2

    def dfpb(self,t):
        return super().dfpb(t) - self.Dws[5]/self.ws[-1]/2

    def ddfp(self,t):
        return (-self.ws[0]**2 + 3/4 * (self.Dws[0]/self.ws[0])**2 - 1/2 *
        self.DDws[0]/self.ws[0])

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver3(RKWKBSolver2):

    def fp(self,t0,t1):
        return super().fp(t0,t1) * numpy.exp(1j * self.S2(t1-t0))

    def dfp(self,t):
        return (super().dfp(t) + 1j * ( 3/8 * self.Dws[0]**2/self.ws[0]**3 - 1/4
        * self.DDws[0]/self.ws[0]**2))

    def dfpb(self,t):
        return (super().dfpb(t) + 1j * (3/8 * self.Dws[5]**2/self.ws[-1]**3 - 1/4
        * self.DDws[-1]/self.ws[-1]**2))

    def ddfp(self,t):
        return (-self.ws[0]**2  - 1/4 * 1j * self.DDDws[0]/self.ws[0]**2 - 3/2 *
        1j * self.Dws[0]**3/self.ws[0]**4 + 3/2 * 1j * self.Dws[0] *
        self.DDws[0]/self.ws[0]**3 - 9/64 * self.Dws[0]**4/self.ws[0]**6 + 3/16
        * self.Dws[0]**2 * self.DDws[0] / self.ws[0]**5 - 1/16 *
        self.DDws[0]**2/self.ws[0]**4)

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))

class RKWKBSolver4(RKWKBSolver3):

    def fp(self,t0,t1):
        return super().fp(t0,t1) * numpy.exp(self.S3(t1-t0)) 

    def dfp(self,t):
        return super().dfp(t) + 3/4.0*self.Dws[0]**3/self.ws[0]**5 - 3/4.0*self.Dws[0]*self.DDws[0]/self.ws[0]**4 + 1/8.0*self.DDDws[0]/self.ws[0]**3 

    def dfpb(self,t):
        return super().dfpb(t) + 3/4.0*self.Dws[-1]**3/self.ws[5]**5 - 3/4.0*self.Dws[-1]*self.DDws[-1]/self.ws[5]**4 + 1/8.0*self.DDDws[-1]/self.ws[5]**3

    def ddfp(self,t):
        return -3/16.0*self.Dws[0]*self.DDws[0]*self.DDDws[0]/self.ws[0]**7 - 297/64.0*self.Dws[0]**4/self.ws[0]**6 + 1/8.0*self.DDDDws[0]/self.ws[0]**3 - 3/16.0*self.DDws[0]**2/self.ws[0]**4 + 9/16.0*self.Dws[0]**6/self.ws[0]**10 + 1/64.0*self.DDDws[0]**2/self.ws[0]**6 + 99/16.0*self.Dws[0]**2*self.DDws[0]/self.ws[0]**5 - 5/4.0*self.Dws[0]*self.DDDws[0]/self.ws[0]**4 - self.ws[0]**2 + 9/16.0*self.Dws[0]**2*self.DDws[0]**2/self.ws[0]**8 - 9/8.0*self.Dws[0]**4*self.DDws[0]/self.ws[0]**9 + 3/16.0*self.Dws[0]**3*self.DDDws[0]/self.ws[0]**8 - 1/16.0*1j*self.DDws[0]*self.DDDws[0]/self.ws[0]**5 + 3/32.0*1j*self.Dws[0]**2*self.DDDws[0]/self.ws[0]**6 + 9/16.0*1j*self.Dws[0]**5/self.ws[0]**8 + 3/8.0*1j*self.DDws[0]**2*self.Dws[0]/self.ws[0]**6 - 15/16.0*1j*self.DDws[0]*self.Dws[0]**3/self.ws[0]**7

    def fm(self,t0,t1):
        return numpy.conj(self.fp(t0,t1))

    def dfm(self,t):
        return numpy.conj(self.dfp(t))

    def dfmb(self,t):
        return numpy.conj(self.dfpb(t))

    def ddfm(self,t):
        return numpy.conj(self.ddfp(t))
   
