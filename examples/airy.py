import pyoscode
import numpy as np
from scipy.special import airy
from matplotlib import pyplot as plt

# Define the frequency and friction term over the range of integration
#ts = np.logspace(0,2,5000)
ts = np.linspace(1,100,5000)
ws = np.sqrt(ts)
gs = np.zeros_like(ws)
# Define the range of integration and the initial conditions
ti = 1.0
tf = 100.0
x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
t_eval = np.linspace(ti,ti+1.0,10)
# Solve the system
# Test forward integration
#sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
# Test backward integration 
#x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # X
t = np.asarray(sol['t'])
x = np.asarray(sol['sol'])
dense = np.asarray(sol['x_eval'])
print("dense output: ",dense[0])
types = np.asarray(sol['types'])
# Plot the solution
ana_t = np.linspace(ti,tf,5000)
plt.plot(ana_t,[1j*airy(-T)[0] - airy(-T)[2] for T in ana_t],color='black',lw=0.5,label='true solution')
plt.plot(t[types==0],1j*x[types==0],'.',color='C0',label='RK steps')
plt.plot(t[types==1],1j*x[types==1],'.',color='C1',label='WKB steps')
plt.plot(t_eval, 1j*dense, color='C1', label='dense output')
plt.legend()
#plt.xlim((1.0,35.0))
plt.ylim((-1.0,1.0))
plt.xlabel('t')
plt.ylabel('Ai(-t)')
plt.show()
#plt.savefig('airy-example.png')

