import matplotlib.pyplot as plt
import numpy as np

# Code to plot results from solving the burst equation with oscode in C++,
# written in "output.txt"
n = 100
f = "output.txt"
data = np.genfromtxt(f, delimiter=', ')
t = data[:,0]
x = data[:,1]
types = data[:,3]
# Analytic solution (real part)
ts = np.linspace(-200,200,40000)
xs = 100.0*(1+ts**2)**0.5/100.0*(np.cos(100.0*np.arctan(ts)))
fig, ax = plt.subplots(1,2)

ax[0].plot(ts,xs,label='true solution',lw=1)
ax[0].plot(t[types==1],x[types==1],'.',color='green',label='WKB steps')
ax[0].plot(t[types==0],x[types==0],'.',color='red',label='RK steps')
ax[0].set_xlabel('t')
ax[0].set_ylabel('Re[x(t)]')

ax[1].plot(ts,xs,label='true solution',lw=1)
ax[1].plot(t[types==1],x[types==1],'.',color='green',label='WKB steps')
ax[1].plot(t[types==0],x[types==0],'.',color='red',label='RK steps')
ax[1].set_xlim((-2,2))
ax[1].set_ylim((-10,10))
ax[1].legend()

plt.tight_layout()
plt.show()
#plt.savefig('burst-example.png')
