import pyoscode
import numpy
from scipy.special import airy
from matplotlib import pyplot as plt

# Define the frequency and friction term over the range of integration
ts = numpy.linspace(1,1000,5000)
ws = numpy.sqrt(ts)
gs = numpy.zeros_like(ws)
# Define the range of integration and the initial conditions
ti = 1.0
tf = 1000.0
x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
# Solve the system
sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0)
t = numpy.asarray(sol['t'])
x = numpy.asarray(sol['sol'])
types = numpy.asarray(sol['types'])
# Plot the solution
ana_t = numpy.linspace(1,35.0,1000)
plt.plot(ana_t,[airy(-T)[0] for T in ana_t],label='true solution')
plt.plot(t[types==0],x[types==0],'.',color='red',label='RK steps')
plt.plot(t[types==1],x[types==1],'.',color='green',label='WKB steps')
plt.legend()
plt.xlim((1.0,35.0))
plt.ylim((-1.0,1.0))
plt.xlabel('t')
plt.ylabel('Ai(-t)')
plt.savefig('airy-example.png')

