import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import re

# Airy equation with w(t), g(t) given analytically
f1 = 'test/test-airy-vhighprec.txt'
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair, 4:parse_pair, 5:parse_pair})
times = np.logspace(0,2,5000)
analytic = np.array([sp.airy(-ti)[0] + 1j*sp.airy(-ti)[2] for ti in times ])
danalytic = np.array([-sp.airy(-ti)[1] - 1j*sp.airy(-ti)[3] for ti in times ])
t = data[:,0]
x = data[:,1]
dx = data[:,2]
wkb = data[:,3]
sols = data[:,4]
dsols = data[:,5]
#sols = np.array([sp.airy(-ti)[0] + 1j*sp.airy(-ti)[2] for ti in t ])
errs = np.abs(sols - x)/np.abs(sols)
derrs = np.abs(dsols - dx)/np.abs(dsols)
fig,axes=plt.subplots(1,2)

axes[0].loglog(t, derrs, 'x', color='black')

axes[0].semilogx(times,analytic, color='black')
axes[0].semilogx(t[wkb==1],x[wkb==1],'x',color='green')
axes[0].semilogx(t[wkb==0],x[wkb==0],'x',color='red')

axes[0].semilogx(times,danalytic)
axes[0].semilogx(t[wkb==1],dx[wkb==1],'x',color='green')
axes[0].semilogx(t[wkb==0],dx[wkb==0],'x',color='red')

axes[0].set_title('Example solution x(t)')
axes[1].loglog(t, errs,'x',color='black')
axes[1].set_title('Relative error')
plt.show()


# To parse complex numbers of the form (a,b) from file
pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
def parse_pair(s):
    return complex(*map(float, pair.match(s).groups()))

# Burst equation with w(t), g(t) given analytically
f1 = 'test/test-burst-n1e5'
n =  1e5 
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair})
times = np.linspace(-1.5*n,1.5*n,5000)
analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])
danalytic = np.array([100.0/np.sqrt(ti**2+1)/n*( (ti + 1j*n ) * np.cos(n*np.arctan(ti)) + (1j*ti - n ) * np.sin(n*np.arctan(ti))) for ti in times ])
t = data[:,0]
x = data[:,1]
dx = data[:,2]
wkb = data[:,3]
sols = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t])
dsols = np.array([100.0/np.sqrt(ti**2+1)/n*( (ti + 1j*n ) * np.cos(n*np.arctan(ti)) + (1j*ti - n ) * np.sin(n*np.arctan(ti))) for ti in t])
errs = np.abs(sols - x)/np.abs(sols)
derrs = np.abs(dsols - dx)/np.abs(dsols)
fig,axes=plt.subplots(1,2)

axes[0].loglog(t, derrs, 'x', color='black')

axes[0].plot(times,analytic, color='black')
axes[0].plot(t[wkb==1],x[wkb==1],'x',color='green')
axes[0].plot(t[wkb==0],x[wkb==0],'x',color='red')

axes[0].semilogx(times,danalytic)
axes[0].semilogx(t[wkb==1],dx[wkb==1],'x',color='green')
axes[0].semilogx(t[wkb==0],dx[wkb==0],'x',color='red')

axes[0].set_title('Example solution x(t)')
axes[1].semilogy(t, errs,'x',color='black')
axes[1].set_title('Relative error')
plt.show()















data = np.genfromtxt('test/timepps2_cont.txt',dtype=complex, filling_values=0.0)
data2 = np.genfromtxt('test/timepps2.txt',dtype=complex,filling_values=0.0)
data3 = np.genfromtxt('test/timepps1.txt',dtype=complex,filling_values=0.0)
data4 = np.genfromtxt('test/timepps1cont.txt',dtype=complex,filling_values=0.0)
fig,axes = plt.subplots(2,2,sharex=True)
axes[0,0].loglog(data[:-110,0],data[:-110,-1],color='blue',label='rst,cubic',alpha=0.5)
axes[0,0].loglog(data2[:,0],data2[:,-1],color='blue',alpha=0.5)
axes[0,0].loglog(data[:-110,0],data[:-110,-2],color='orange',label='hd,cubic',alpha=0.7)
axes[0,0].loglog(data2[:,0],data2[:,-2],color='orange',alpha=0.7)
axes[0,0].loglog(data3[:,0],data3[:,-1],color='green',label='rst,linear',alpha=0.5)
axes[0,0].loglog(data4[:,0],data4[:,-1],color='green',alpha=0.5)
axes[0,0].loglog(data3[:,0],data3[:,-2],color='red',label='hd,linear',alpha=0.7)
axes[0,0].loglog(data4[:,0],data4[:,-2],color='red',alpha=0.7)

axes[0,1].loglog(data[:-110,0],data[:-110,3],label='total steps, cubic',alpha=0.7)
axes[0,1].loglog(data2[:,0],data2[:,3],alpha=0.7,color='orange')
axes[0,1].loglog(data[:,0],data[:,4],label='wkb steps, cubic',alpha=0.7,color='blue') 
axes[0,1].loglog(data2[:,0],data2[:,4],alpha=0.7,color='blue')
axes[0,1].loglog(data3[:,0],data3[:,3],label='total steps, linear',alpha=0.7,color='red') 
axes[0,1].loglog(data4[:,0],data4[:,3],alpha=0.7,color='red')
axes[0,1].loglog(data3[:,0],data3[:,4],label='wkb steps, linear',alpha=0.7,color='green') 
axes[0,1].loglog(data4[:,0],data4[:,4],alpha=0.7,color='green')
axes[0,1].set_title('Step breakdown')
axes[0,1].legend()

axes[1,0].loglog(data[:,0],data[:,5],color='blue',label='cubic spline',alpha=0.5)
axes[1,0].loglog(data2[:,0],data2[:,5],color='blue',alpha=0.5)
axes[1,0].loglog(data3[:,0],data3[:,5],color='orange',alpha=0.5)
axes[1,0].loglog(data4[:,0],data4[:,5],color='orange',label='linear',alpha=0.5)
axes[1,0].legend()
axes[1,0].set_title('Total time')
axes[1,0].set_ylabel('time/s')

axes[0,0].set_title('PPS with different interpolation methods')
axes[0,0].set_xlabel('k (Planck units)')
axes[0,0].set_ylabel('$P_{\mathcal{R}_k}$')
axes[0,0].legend()
#plt.savefig('test/stepstot_cubic.pdf')
plt.show()


data = np.genfromtxt('test/timeinterp_cubic.txt')
data2 = np.genfromtxt('test/timeinterp_linear.txt')
plt.xlabel('total number of points, N')
plt.ylabel('time/s')
plt.title('Timing test of interpolating methods')
plt.loglog(data[:,0],data[:,1],label='to create interp object, cubic')
plt.loglog(data[:,0],data[:,2],label='to call, cubic')
plt.loglog(data2[:,0],data2[:,1],label='to create interp object, linear')
plt.loglog(data2[:,0],data2[:,2],label='to call, linear')
plt.legend()
plt.savefig('plots/interptimetest.pdf')
plt.show()
