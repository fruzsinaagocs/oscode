import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import re

# To parse complex numbers of the form (a,b) from file
pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
def parse_pair(s):
    return complex(*map(float, pair.match(s).groups()))

# Airy equation with w(t), g(t) given analytically
f1 = 'test/airy/airycorr_o1.txt'
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
axes[1].loglog(t[wkb==1], errs[wkb==1],'x',color='green')
axes[1].loglog(t[wkb==0], errs[wkb==0],'x',color='red')
axes[1].set_title('Relative error')
plt.show()


# Burst equation with w(t), g(t) given analytically
# Single solution at low n 
f1 = 'test/burst/d4w1less_n40.txt' #'plots/burstn40.txt'
n = 40.0 
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair})
times = np.linspace(-2*n,2*n,10000)
analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])
t = data[:,0]
x = data[:,1]
dx = data[:,2]
wkb = data[:,3]
sols = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t])

plt.figure()
plt.style.use('fyr')
plt.plot(times,analytic, color='black',label='true solution')
plt.plot(t[wkb==1],x[wkb==1],'x',color='green',label='WKB step')
plt.plot(t[wkb==0],x[wkb==0],'x',color='red',label='RK step')
#plt.xlim((-40,40))
#plt.ylim((-60,60))
plt.xlabel('t')
plt.ylabel('$\Re{\{x(t)\}}$')
plt.legend()
#plt.savefig("plots/burstn40_x_fyrstyle.pdf")
plt.show()

# Error progression in this example
errs = np.abs(sols - x)/np.abs(sols)
plt.figure()
plt.style.use('fyr')
plt.semilogy(t,errs,'x-',color='black')
plt.ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
plt.xlabel('t')
#plt.savefig("plots/burstn40_err_fyrstyle.pdf")
plt.show()

# Large n example: number of oscillations traversed
f = "plots/burstn1e5_tol-4.txt"
n = 1e5
data = np.loadtxt(f,dtype=complex,converters={1:parse_pair, 2:parse_pair})
ts = data[:,0]
xs = data[:,1]
wkbs = data[:,-1]
oscs = np.zeros(ts.size) 
for i,t in enumerate(ts): 
    if i!=0:
        oscs[i] = ((n*np.arctan(t)/(2*np.pi))-(n*np.arctan(ts[i-1])/(2*np.pi)))
    else:
        oscs[i] = None

plt.figure();
plt.style.use("fyr")
plt.semilogy(ts,oscs,color='black')
plt.semilogy(ts[wkbs==1],oscs[wkbs==1],'.',label='WKB step',color='green')
plt.semilogy(ts[wkbs==0],oscs[wkbs==0],'.',label='RK step',color='red')
plt.xlim((-n,n))
plt.xlabel('t')
plt.ylim((1e-2,1e4))
plt.ylabel('oscillations traversed')
plt.legend()
plt.savefig('plots/burstn1e5_oscs.pdf')

# Changing rtol at n=1e5, effect on relative error progression 
f1 = 'plots/burstn1e5_tol-4.txt'
f2 = 'plots/burstn1e5_tol-5.txt'
f3 = 'plots/burstn1e5_tol-6.txt'
n = 1e5 
data1 = np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair})
data2 = np.loadtxt(f2,dtype=complex,converters={1:parse_pair, 2:parse_pair})
data3 = np.loadtxt(f3,dtype=complex,converters={1:parse_pair, 2:parse_pair})
t1 = data1[:,0]
x1 = data1[:,1]
wkb1 = data1[:,3]
sols1 = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t1])
errs1 = np.abs(sols1 - x1)/np.abs(sols1)
#
t2 = data2[:,0]
x2 = data2[:,1]
wkb2 = data2[:,3]
sols2 = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t2])
errs2 = np.abs(sols2 - x2)/np.abs(sols2)
#
t3 = data3[:,0]
x3 = data3[:,1]
wkb3 = data3[:,3]
sols3 = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t3])
errs3 = np.abs(sols3 - x3)/np.abs(sols3)
#
plt.figure()
plt.style.use('fyr')
plt.semilogy(t1,errs1,'-',color='black',label='rtol=$10^{-4}$')
plt.semilogy(t2,errs2,'-.',color='black',label='rtol=$10^{-5}$')
plt.semilogy(t3,errs3,'--',color='black',label='rtol=$10^{-6}$')
plt.ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
plt.xlabel('t')
plt.legend()
plt.savefig("plots/burstn1e5_rtols.pdf")

# Effect of changing rtol and n on runtime
f1 = 'plots/bursttimingtol-4.txt'
f2 = 'plots/bursttimingtol-5.txt'
f3 = 'plots/bursttimingtol-6.txt'
data1 = np.loadtxt(f1, delimiter=',')
data2 = np.loadtxt(f2, delimiter=',')
data3 = np.loadtxt(f3, delimiter=',')
ns = data1[:,0]
t1 = data1[:,-1]
t2 = data2[:,-1]
t3 = data3[:,-1]
pivot = t1[249]
t1 = t1/pivot
t2 = t2/pivot
t3 = t3/pivot
plt.figure()
plt.style.use('fyr')
plt.loglog(10**ns,t1,'-',color='black',label='rtol=$10^{-4}$')
plt.loglog(10**ns,t2,'-.',color='black',label='rtol=$10^{-5}$')
plt.loglog(10**ns,t3,'--',color='black',label='rtol=$10^{-6}$')
plt.loglog([10**ns[249],10**ns[249]],[0.3,1.0],'k:')
plt.loglog([4.0,10**ns[249]],[1.0, 1.0], 'k:')
plt.ylabel('relative runtime')
plt.xlabel('n')
plt.xlim((4.0,2.5*10**(10)))
plt.ylim((0.3,4.5))
plt.legend()
plt.savefig("plots/bursttiming.pdf")

# Effect of changing n, rtol on number of steps
f1 = 'plots/bursttimingtol-4_stepscorr.txt'
f2 = 'plots/bursttimingtol-5_stepscorr.txt'
f3 = 'plots/bursttimingtol-6_stepscorr.txt'
data1 = np.loadtxt(f1, delimiter=',')
data2 = np.loadtxt(f2, delimiter=',')
data3 = np.loadtxt(f3, delimiter=',')
ns = data1[:,0]
ts1 = data1[:,1]
#ts2 = data2[:,1]
#ts3 = data3[:,1]
ss1 = data1[:,2]
#ss2 = data2[:,2]
#ss3 = data3[:,2]
wkbs1 = data1[:,3]
#wkbs2 = data2[:,3]
#wkbs3 = data3[:,3]
plt.figure()
plt.style.use('fyr')
plt.loglog(10**ns,ts1,':',color='black',label='total steps')
plt.loglog(10**ns,ss1,'-',color='black',label='accepted steps')
plt.loglog(10**ns,wkbs1,'-',color='green',label='WKB steps')
plt.loglog(10**ns,ss1-wkbs1,'-',color='red',label='RK steps')
#plt.loglog(10**ns,t2,'-.',color='black',label='rtol=$10^{-5}$')
#plt.loglog(10**ns,t3,'--',color='black',label='rtol=$10^{-6}$')
#plt.loglog([10**ns[249],10**ns[249]],[0.3,1.0],'k:')
plt.ylabel('number of steps')
plt.xlabel('n')
plt.legend()
plt.savefig("plots/burststeps.pdf")

# Version 2 of the above plot - vary rtol rather than n
from cycler import cycler
grays = cycler(color=['0.00','0.20', '0.40', '0.50', '0.60', '0.70'])
f1 = 'plots/bursttimingn1e1.txt'
f2 = 'plots/bursttimingn1e2.txt'
f3 = 'plots/bursttimingn1e3.txt'
f4 = 'plots/bursttimingn1e4.txt'
f5 = 'plots/bursttimingn1e5.txt'
f6 = 'plots/bursttimingn1e6.txt'
data1 = np.loadtxt(f1, delimiter=',')
data2 = np.loadtxt(f2, delimiter=',')
data3 = np.loadtxt(f3, delimiter=',')
data4 = np.loadtxt(f4, delimiter=',')
data5 = np.loadtxt(f5, delimiter=',')
data6 = np.loadtxt(f6, delimiter=',')
rtols = data1[:,0]
t1 = data1[:,1]
t2 = data2[:,1]
t3 = data3[:,1]
t4 = data4[:,1]
t5 = data5[:,1]
t6 = data6[:,1]
pivot_t = t1[199]
pivot_rtol = rtols[199]
fig, ax = plt.subplots(1,1)
plt.style.use('fyr')
ax.set_prop_cycle(grays)
plt.loglog(10**rtols,t1/pivot_t,label='n=$10^{1}$')
plt.loglog(10**rtols,t2/pivot_t,label='n=$10^{2}$')
plt.loglog(10**rtols,t3/pivot_t,label='n=$10^{3}$')
plt.loglog(10**rtols,t4/pivot_t,label='n=$10^{4}$')
plt.loglog(10**rtols,t5/pivot_t,label='n=$10^{5}$')
plt.loglog(10**rtols,t6/pivot_t,label='n=$10^{6}$')
plt.loglog([10**pivot_rtol, 10**pivot_rtol],[1e-1,1.0],':',color='black')
plt.loglog([1e-7, 10**pivot_rtol],[1.0,1.0],':',color='black')
plt.ylabel('relative runtime')
plt.xlabel("relative tolerance `rtol'")
plt.legend()
plt.show()
#plt.savefig('plots/bursttiming_rtol.pdf')

# Mukhanov--Sasaki equation
# NAG files: ms1 has k=0.1, ms2 k=10.0, ms3 k=0.03

fsolver = "test/ms/ms-singlek-k5e-1.txt"
fnag = "test/ms/nag-ms-k5e-1.txt"
dsolver = np.loadtxt(fsolver,dtype=complex,converters={1:parse_pair, 2:parse_pair})
dnag = np.loadtxt(fnag,delimiter=",")
t = dsolver[:,0]
x = dsolver[:,1]
dx = dsolver[:,2]
wkb = dsolver[:,3]
tnag = dnag[:,0]
xnag = dnag[:,1]

plt.figure()
plt.style.use('fyr')
plt.semilogx(tnag,xnag, color='black',label='true solution', lw=0.5)
plt.semilogx(t[wkb==1],x[wkb==1],'x',color='green',label='WKB step')
plt.semilogx(t[wkb==0],x[wkb==0],'x',color='red',label='RK step')
plt.xlim((9900,2e5))
#plt.ylim((-60,60))
plt.xlabel('t')
plt.ylabel('$\Re{\{x(t)\}}$')
plt.legend()
plt.savefig("plots/ms-k51e-1.pdf")
#plt.show()

# PPS
f = "test/ms/pps-timed.txt"
data = np.genfromtxt(f,delimiter=",",dtype=float,missing_values='nan',filling_values=0.0)
k = data[:,0]
phd = data[:,1]
prst = data[:,2]
plt.figure()
plt.style.use("fyr")
plt.xlabel("k")
plt.ylabel("$P_{\mathcal{R}}(k)$")
plt.loglog(k, phd, '-',label='HD',color='gray')
plt.loglog(k, prst, '-',label='RST')
plt.legend()
#plt.show()
plt.savefig("plots/example-pps.pdf")

# Timing comparison with BINGO 
from cycler import cycler
grays = cycler(color=['0.00','0.40', '0.60', '0.50', '0.60', '0.70'])

fsolver = "test/ms/pps-timed.txt"
fbingo1 = "plots/bingo-pps-1e2.txt"
fbingo2 = "plots/bingo-pps-1e3.txt"
dsolver = np.genfromtxt(fsolver,delimiter=",")
dbingo1 = np.genfromtxt(fbingo1)
dbingo2 = np.genfromtxt(fbingo2)
k = dsolver[:,0]
t = dsolver[:,-1]
kbingo1 = dbingo1[:,0]
tbingo1 = dbingo1[:,-1]
kbingo2 = dbingo2[:,0]
tbingo2 = dbingo2[:,-1]

fig,ax=plt.subplots(1,1)
ax.set_prop_cycle(grays)
plt.style.use('fyr')
plt.loglog(kbingo1,tbingo1,'.',label='BINGO - $k/aH=10^2$ start')
plt.loglog(kbingo2,tbingo2,'.',label='BINGO - $k/aH=10^3$ start')
plt.loglog(k,t,'.',label='RKWKB')
plt.xlabel('k')
plt.ylabel('runtime/s')
plt.legend()
plt.show()
plt.savefig("plots/bingo-rkwkb-abst.pdf")

fig,ax=plt.subplots(1,1)
ax.set_prop_cycle(grays)
plt.style.use('fyr')
plt.loglog(kbingo1,tbingo1/t,'.',label='$k/aH=10^2$ start')
plt.loglog(kbingo1,tbingo2/t,'.',label='$k/aH=10^3$ start')
plt.xlabel('k')
plt.ylabel('$t_{\mathrm{BINGO}}/t_{\mathrm{RKWKB}}$')
plt.legend()
plt.show()
plt.savefig("plots/bing-rkwkb-relt.pdf")

# Plot of quadratic potential with a step
# Params:
# phi_0: location of step, BINGO default: 14.668
# alpha: height of step, BINGO: 0.001606
# delta_phi: width of step, BINGO: 0.031108
# m: mass of inflaton field, BINGO: 7.15e-6

m = 7.15e-6
phi_0 = 23.0
alpha = 0.002
delta_phi = 0.03
phi = np.linspace(25,0,5000)
V = 0.5*m**2*phi**2*(1.0+alpha*np.tanh((phi-phi_0)/delta_phi))
fig,ax=plt.subplots(1,1)
plt.plot(phi, V)
plt.show()

# PPS with quadratic step in potential
f1 = "test/ms/pps-quadratic-hankel-3-rtol5.txt"
f2 = "test/ms/pps-quadstep-hankel-3-rtol5_opt.txt"
data1 = np.loadtxt(f1,delimiter=", ",dtype=complex,converters={1:parse_pair, 2:parse_pair})
data2 = np.loadtxt(f2,delimiter=", ",dtype=complex,converters={1:parse_pair, 2:parse_pair})
k1 = data1[:,0]
phd1 = data1[:,3]
prst1 = data1[:,4]
pkd1 = data1[:,5]
k2 = data2[:,0]
phd2 = data2[:,3]
prst2 = data2[:,4]
pkd2 = data2[:,5]
plt.figure()
plt.style.use("default")
plt.xlabel("k")
plt.ylabel("$P_{\mathcal{R}}(k)$")
#plt.loglog(k1, phd1, '-',label='hd - quadratic')
#plt.loglog(k2, phd2, '-',label='hd - step')
#plt.loglog(k1, prst1, '-',label='rst - quadratic')
#plt.loglog(k2, prst2, '-',label='rst - step')
#plt.loglog(k1, pkd1, label='kd - quadratic')
plt.loglog(k2, pkd2, label='kd - optimised')
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig("plots/pps-kd-step-comparison.pdf")

# Background
f = "test/ms/pps-quadstep-3-bg.txt"
data = np.genfromtxt(f,delimiter=',')
t = data[:,0]
phi = data[:,1]
plt.semilogx(t,phi)
plt.show()








#axes[0].semilogx(times,danalytic)
#axes[0].semilogx(t[wkb==1],dx[wkb==1],'x',color='green')
#axes[0].semilogx(t[wkb==0],dx[wkb==0],'x',color='red')
axes[1].semilogy(t[wkb==1], errs[wkb==1],'x',color='green')
axes[1].semilogy(t[wkb==0], errs[wkb==0],'x',color='red')
axes[1].set_title('Relative error')
dsols = np.array([100.0/np.sqrt(ti**2+1)/n*( (ti + 1j*n ) * np.cos(n*np.arctan(ti)) + (1j*ti - n ) * np.sin(n*np.arctan(ti))) for ti in t])
errs = np.abs(sols - x)/np.abs(sols)
derrs = np.abs(dsols - dx)/np.abs(dsols)

