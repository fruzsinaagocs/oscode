import pyoscode
import numpy as np
from scipy.special import airy
from matplotlib import pyplot as plt

# Define the frequency and friction term over the range of integration
ts = np.linspace(1,100,5000)
ws = np.sqrt(ts)
gs = np.zeros_like(ws)

# Define the frequency and friction term as functions
def w(t):
    return np.sqrt(t)

def g(t):
    return 0.0

# Helper function for plotting
def plot_airy(sol):
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    types = np.asarray(sol['types'])
    # Plot the solution
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.array([airy(-T)[0] + 1j*airy(-T)[2] for T in ana_t])
    plt.plot(ana_t,ana_x,color='black',lw=0.5,label='true solution')
    plt.plot(t[types==0],x[types==0],'.',color='C0',label='RK steps')
    plt.plot(t[types==1],x[types==1],'.',color='C1',label='WKB steps')
    plt.plot(t_eval, dense, color='C1', label='dense output')
    plt.legend()
    plt.ylim((-1.0,1.0))
    plt.xlabel('t')
    plt.ylabel('Ai(-t)')
    plt.show()
    #plt.savefig('airy-example.png')


# Define the range of integration and the initial conditions
ti = 1.0
tf = 100.0
x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
dx0 = -airy(-ti)[1] -1j*airy(-ti)[3]
t_eval = np.linspace(ti+20,ti+40,500)

# Solve the system: w, g, arrays
# Test forward integration
sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
plot_airy(sol)

#Solve the system: w, g functions
sol = pyoscode.solve_fn(w, g, ti, tf, x0, dx0, t_eval=t_eval)
plot_airy(sol)




