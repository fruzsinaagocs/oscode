import _solver
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt

l=1.0
m=0.5
#E=1.065286

def Ei(n):
    """ 
    Energy eigenvalues (if analytic solution available)
    """
    return np.sqrt(2.0)*(n-0.5)

def V(t):
    """ 
    Potential well 
    """
    return t**2 + l*t**4

def w(t,E):
    """
    Frequency term in the Schrodinger equation
    """
    return np.sqrt(2*m*(complex(E)-V(t)));

def f(E):
    """
    Function to minimize wrt E to give the energy eigenvalues
    """

    # Boundaries of integration
    tl = -((E)**(0.25))-2.0
    tr = -tl
    tm = 0.5
    
    # Grid of w, g
    t = np.linspace(tl.real,tr.real,30000)
    ws = np.log(w(t,E))
    g = np.zeros(t.shape)
    sol_l = _solver.solve(t,ws,g,tl,tm,0,1e-3,logw=True,rtol=1e-5)
    sol_r = _solver.solve(t,ws,g,tr,tm,0,1e-3,h=-1,logw=True,rtol=1e-5)
    psi_l = sol_l["sol"][-1] 
    psi_r = sol_r["sol"][-1] 
    dpsi_l = sol_l["dsol"][-1]
    dpsi_r = sol_r["dsol"][-1]
#    print(sol_l["types"])
#    print(sol_l["t"])
    try:
        return abs(dpsi_l/psi_l - dpsi_r/psi_r)
    except ZeroDivisionError:
        return 1000.0

bounds = [
(0,2),(3,5),(8,9),(12,14),(17,19),(88,89),(95,97),(103,105),(110,112),(118,120),(416.5,417.5),(1035,1037),(21930,21940),(471100,471110)]
#print("Energy levels:")
#ns = range(1,10)
#Es = [Ei(n) for n in ns]
#for e,n in zip(Es,ns):
#    print("n = "+str(n)+", E="+ str(e))
#print("\n")
#res = minimize_scalar(f,bounds=(1.5*800,1.5*801),method='bounded')
for bound in bounds:
    res = minimize_scalar(f,bounds=bound,method='bounded')
    print(res.x,f(res.x),res.success,res.message)
#res = minimize_scalar(f,bounds=bounds[-1],method='bounded')
#print(res.x,f(res.x))
#print(f(471103.70))
#print(f(471103.76))
#print(f(471103.77))
#print(f(471103.78))
#print(f(471103.79))
#print(f(471103.80))
#print(f(471103.90))
#print(f(471103.85))
#print(f(471104.00))
evalues = np.linspace(471103.96,471104.04,100)
fvalues = [f(evalue) for evalue in evalues]
plt.plot(evalues,fvalues)
plt.show()
#print(f(471103.99))
#print(f(471103.985))
#print(f(471103.98))
#print(f(471103.975))
#print(f(471103.97))

