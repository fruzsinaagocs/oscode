import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import re
import math

# To parse complex numbers of the form (a,b) from file
pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
def parse_pair(s):
    return complex(*map(float, pair.match(s).groups()))

# Airy equation with w(t), g(t) given analytically
# FIGURE 1
f1 = "test/airy/airy-example.txt" #'test/airy/airycorr_o1.txt'
frk = "test/airy/airy-example-rkonly.txt"
datark = np.loadtxt(frk,dtype=complex,converters={1:parse_pair, 2:parse_pair,
4:parse_pair, 5:parse_pair})
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair,
4:parse_pair, 5:parse_pair})
times = np.logspace(0,2,5000)
analytic = np.array([sp.airy(-ti)[0] + 1j*sp.airy(-ti)[2] for ti in times ])
danalytic = np.array([-sp.airy(-ti)[1] - 1j*sp.airy(-ti)[3] for ti in times ])
t = data[:,0]
trk = datark[:,0]
x = data[:,1]
xrk = datark[:,1]
dx = data[:,2]
wkb = data[:,3]
wkbrk = datark[:,3]
sols = data[:,4]
solsrk = datark[:,4]
dsols = data[:,5]
errs = np.abs(sols - x)/np.abs(sols)
derrs = np.abs(dsols - dx)/np.abs(dsols)
errsrk = np.abs(solsrk - xrk)/np.abs(solsrk)
# Example solution
plt.style.use('paper')
fig,ax=plt.subplots(1,2,figsize=(5.95,2.375))
ax[0].semilogx(times,analytic, color='black',label='true solution',lw=0.5)
ax[0].semilogx(t[wkb==1],x[wkb==1],'^',ms=2.5,color='green',label='WKB step')
ax[0].semilogx(t[wkb==0],x[wkb==0],'.',color='red', label='RK step')
ax[0].set_xlim((0,40.0))
ax[0].set_ylim((-0.5, 0.7))
ax[0].set_xlabel("$t$")
ax[0].set_ylabel("$\Re\{x\}$")
ax[0].text(-0.2,1.1, '(a)', transform=ax[0].transAxes)
ax[0].legend()
# Error progression
ax[1].loglog(t, errs, '-', color='black')
ax[1].loglog(t[wkb==1], errs[wkb==1],'^',ms=2.5,color='green',label='WKB step')
ax[1].loglog(t[wkb==0], errs[wkb==0],'.',color='red',label='RK step')
ax[1].loglog(trk[wkbrk==0], errsrk[wkbrk==0],'.',color='red')
ax[1].loglog(trk, errsrk, '-', color='black')
ax[1].set_ylim((1e-6,1e-2))
ax[1].set_xlabel("$t$")
ax[1].set_ylabel("relative error, $\\frac{|\Delta x|}{|x|}$")
ax[1].text(-0.22,1.1, '(b)', transform=ax[1].transAxes)
#plt.axvline(x=3.79)
#plt.text(3.79/(7.0/3.79), 2e-4, 'RK', verticalalignment='center',
#horizontalalignment='right')
#plt.text(7.0, 2e-4, 'WKB', verticalalignment='center',
#horizontalalignment='left')
ax[1].legend()
#plt.show()
plt.tight_layout()
#plt.show()
plt.savefig("plots/airy-merged-dots-labelled.pdf")

# Burst equation with w(t), g(t) given analytically
# Single solution at n=40 and its error progression on one plot
# FIGURE 2  
f1 = 'plots/burstn40.txt'
f2 = 'plots/burst40_rkonly.txt' # Pure RK 
n = 40 
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair})
data2 = np.loadtxt(f2,dtype=complex,converters={1:parse_pair, 2:parse_pair})
times = np.linspace(-2*n,2*n,10000)
analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])
# RKWKB
t = data[:,0]
x = data[:,1]
dx = data[:,2]
wkb = data[:,3]
sols = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t])
errs = np.abs(sols - x)/np.abs(sols)
# Pure RK
t2 = data2[:,0]
x2 = data2[:,1]
dx2 = data2[:,2]
wkb2 = data2[:,3]
sols2 = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) +
np.cos(n * np.arctan(ti))) for ti in t2])
errs2 = np.abs(sols2 - x2)/np.abs(sols2)
# Plotting
plt.style.use('paper')
fig,ax = plt.subplots(1,2,figsize=(5.95,2.375))
ax[0].plot(times,analytic, color='black',label='true solution')
ax[0].plot(t[wkb==1],x[wkb==1],'^',color='green',ms=2.5,label='WKB step')
ax[0].plot(t[wkb==0],x[wkb==0],'.',color='red',label='RK step')
ax[0].set_xlim((-40,40))
ax[0].set_ylim((-60,60))
ax[0].set_xlabel('$t$')
ax[0].set_ylabel('$\Re{\{x(t)\}}$')
ax[0].text(-0.2,1.1, '(a)', transform=ax[0].transAxes)
ax[0].legend()
ax[1].semilogy(t2,errs2,color='black')
ax[1].semilogy(t2,errs2,'.',color='red')
ax[1].semilogy(t,errs,color='black')
ax[1].semilogy(t[wkb==1],errs[wkb==1],'^',color='green',ms=2.5,label='WKB step')
ax[1].semilogy(t[wkb==0],errs[wkb==0],'.',color='red',label='RK step')
ax[1].set_ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
ax[1].set_xlabel('$t$')
ax[1].set_ylim(1e-7,1e-1)
ax[1].set_xlim(-60,60)
ax[1].text(-0.23,1.1, '(b)', transform=ax[1].transAxes)
ax[1].legend()
#plt.show()
plt.tight_layout()
plt.savefig("plots/burst_n40merged-dots-labelled.pdf")
#plt.show()

# Large n example: number of oscillations traversed
# FIGURE 4
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

fig,ax = plt.subplots(1,1)
plt.style.use("paper")
plt.semilogy(ts,oscs,color='black')
plt.semilogy(ts[wkbs==1],oscs[wkbs==1],'^',ms=2.5,label='WKB step',color='green')
plt.semilogy(ts[wkbs==0],oscs[wkbs==0],'.',label='RK step',color='red')
plt.xlim((-n,n))
plt.xlabel('$t$')
plt.ylim((1e-2,1e4))
plt.ylabel('oscillations traversed')
ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2),useOffset=False)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('plots/burstn1e5_oscs.pdf')

# Changing rtol at n=1e5, effect on relative error progression 
# FIGURE 3
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
fig,ax = plt.subplots(1,1)
plt.style.use('paper')
plt.semilogy(t1,errs1,'-',color='black',label='rtol=$10^{-4}$')
plt.semilogy(t2,errs2,'-.',color='black',label='rtol=$10^{-5}$')
plt.semilogy(t3,errs3,'--',color='black',label='rtol=$10^{-6}$')
plt.ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
plt.xlabel('$t$')
plt.xlim((-2*n,2*n))
plt.ylim((1e-9, 1e-1))
ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2),useOffset=False)
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("plots/burstn1e5_rtols.pdf")

# Effect of changing rtol and n on runtime
# FIGURE 5
f1 = 'plots/bursttimingtol-4.txt'
f2 = 'plots/bursttimingtol-5.txt'
f3 = 'plots/bursttimingtol-6.txt'
data1 = np.loadtxt(f1, delimiter=',')
data2 = np.loadtxt(f2, delimiter=',')
data3 = np.loadtxt(f3, delimiter=',')
ns = data1[::3,0]
t1 = data1[::3,-1]
t2 = data2[::3,-1]
t3 = data3[::3,-1]
pivot = t1[84]
pivotno = 84 
t1 = t1/pivot
t2 = t2/pivot
t3 = t3/pivot
plt.figure()
plt.style.use('paper')
plt.loglog(10**ns,t1,'-',color='black',label='rtol=$10^{-4}$')
plt.loglog(10**ns,t2,'-.',color='black',label='rtol=$10^{-5}$')
plt.loglog(10**ns,t3,'--',color='black',label='rtol=$10^{-6}$')
plt.loglog([10**ns[pivotno],10**ns[pivotno]],[0.3,1.0],'k:')
plt.loglog([4.0,10**ns[pivotno]],[1.0, 1.0], 'k:')
plt.ylabel('relative runtime')
plt.xlabel('$n$')
plt.xlim((1e1,10**(10)))
plt.ylim((0.3,4.5))
plt.legend()
#plt.show()
plt.tight_layout()
plt.savefig("plots/bursttiming.pdf")

# Effect of changing n, rtol on number of steps
# FIGURE 7
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
plt.style.use('paper')
plt.loglog(10**ns,ts1,':',color='black',label='total steps')
plt.loglog(10**ns,ss1,'-',color='black',label='accepted steps')
plt.loglog(10**ns,wkbs1,'--',color='green',ms=2.5,label='WKB steps')
plt.loglog(10**ns,ss1-wkbs1,'-',color='red',label='RK steps')
#plt.loglog(10**ns,t2,'-.',color='black',label='rtol=$10^{-5}$')
#plt.loglog(10**ns,t3,'--',color='black',label='rtol=$10^{-6}$')
#plt.loglog([10**ns[249],10**ns[249]],[0.3,1.0],'k:')
plt.xlim((1e1,1e10))
plt.ylim((2e1, 1e3))
plt.ylabel('number of steps')
plt.xlabel('$n$')
plt.legend()
#plt.show()
plt.savefig("plots/burststeps.pdf")

# Version 2 of the above plot - vary rtol rather than n
# FIGURE 6
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
plt.style.use('paper')
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
plt.xlim((1e-7, 1e-3))
plt.ylim((2e-1,1e2))
plt.legend()
#plt.show()
plt.tight_layout()
plt.savefig('plots/bursttiming_rtol.pdf')

# Mukhanov--Sasaki equation
# NAG files: ms1 has k=0.1, ms2 k=10.0, ms3 k=0.03

fsolver = "test/ms/ms-singlek-diffgrid3.txt"#"test/ms/ms-singlek-k5e-1.txt"
fnag = "test/ms/nag-ms-kd1.txt" #"test/ms/nag-ms-k5e-1.txt"
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
plt.plot(tnag,xnag, color='black',label='true solution',lw=0.5)
plt.plot(t[wkb==1],x[wkb==1],'^',color='green',ms=2.5,label='WKB step')
plt.plot(t[wkb==0],x[wkb==0],'x',color='red',label='RK step')
#plt.xlim((9900,2e5))
#plt.ylim((-60,60))
plt.xlabel('t')
plt.ylabel('$\Re{\{x(t)\}}$')
plt.legend()
#plt.savefig("plots/ms-k51e-1.pdf")
plt.show()

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
plt.show()
#plt.savefig("plots/example-pps.pdf")

# Timing comparison with BINGO - old and not like-to-like, not used in paper
from cycler import cycler
grays = cycler(color=['0.00','0.40', '0.60', '0.50', '0.60', '0.70'])

fsolver = "test/ms/pps-hd_rst_opt.txt"
fbingo1 = "plots/bingo-pps-k1e3-1e2.txt"
fbingo2 = "plots/bingo-pps-k1e3-1e3.txt"
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
#plt.show()
plt.savefig("plots/bingo-rkwkb-abst_opt.pdf")

fig,ax=plt.subplots(1,1)
ax.set_prop_cycle(grays)
plt.style.use('fyr')
plt.loglog(kbingo1,tbingo1/t,'.',label='$k/aH=10^2$ start')
plt.loglog(kbingo1,tbingo2/t,'.',label='$k/aH=10^3$ start')
plt.xlabel('k')
plt.ylabel('$t_{\mathrm{BINGO}}/t_{\mathrm{RKWKB}}$')
plt.legend()
#plt.show()
plt.savefig("plots/bing-rkwkb-relt_opt.pdf")

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

# PPS with different (HD, RST, KD) initial conditions
f1 = "test/ms/pps-kd_opt.txt"
f2 = "test/ms/pps-hd_rst_opt.txt"#"test/ms/pps-logt-1"#"test/ms/pps-hd_rst_opt.txt"
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
plt.loglog(k2, phd2, '-',label='hd')
plt.loglog(k2, prst2, '-',label='rst')
plt.loglog(k1, pkd1, label='kd')
plt.loglog(k2,pkd2,label='kd - logt')
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig("plots/pps-kd-step-comparison.pdf")

# Background
f = "test/ms/pps-bg.txt"
data = np.genfromtxt(f,delimiter=',')
t = data[:,0]
aH = data[:,1]
phi = data[:,2]
fig, ax = plt.subplots(1,1)
plt.axvline(x=1e4)
plt.loglog(t,aH**(-1)/100, label='BINGO start')
plt.loglog(t,aH**(-1), label='Hubble horizon')
#ax.axhspan(1e-3, 1e3, alpha=0.5, color='green')
plt.xlabel('t')
plt.ylabel('lengthscales')
plt.legend()
plt.show()

#plt.semilogx(t,phi,label='$\phi$')
plt.show()

# Testing grid fineness
plt.style.use("default")
#f = "test/ms/pps-kd-diffgrid2-t10.txt"
f2 = "test/ms/pps-testinggrid9-scaled.txt"
data = np.genfromtxt(f,delimiter=', ',dtype=complex,converters={1:parse_pair, 2:parse_pair})
data2 = np.genfromtxt(f2,delimiter=', ',dtype=complex,converters={1:parse_pair, 2:parse_pair})
k = data[:,0]
hd = data[:,3]
rst = data[:,4]
kd = data[:,5]
k2 = data2[:,0]
hd2 = data2[:,3]
rst2 = data2[:,4]
kd2 = data2[:,5]
fig, ax = plt.subplots(1,1)
#plt.loglog(k,hd,label='hd')
#plt.loglog(k,rst,label='rst')
#plt.loglog(k,kd,label='kd-'+f)
plt.loglog(k2,kd2,label='kd-'+f2)
plt.xlabel('k')
plt.ylabel('$P_{\mathcal{R}}(k)$')
plt.legend()
plt.show()

# BINGO PPS comparison 

# RKWKB: pps-bingo-start1e2-highk-lowk.txt
#        pps-bingo-start2e2-highk-lowk.txt
# BINGO: test/ms/pps-bingo-comparison-ltl-moreac-abs0.txt
#        test/ms/pps-bingo-comparison-start2e2.txt

# FIGURE 9 
from cycler import cycler
from matplotlib.ticker import FuncFormatter
grays = cycler(color=['0.00','0.60', '0.40', '0.60', '0.80', '0.90'])

fsolver = "test/ms/pps-bingo-start1e2-highk-lowk.txt"
#"test/ms/pps-bingo-start1e2-highk-tol4.txt"
fbingo1 = "test/ms/pps-bingo-comparison-start2e2.txt"
fbingo2 = "test/ms/pps-bingo-comparison-ltl-moreacc-pt2.txt"
fbingo3 = "test/ms/pps-bingo-comparison-ltl-moreacc-pt3.txt"
fbingo4 = "test/ms/pps-bingo-comparison-ltl-moreacc-abs0.txt"

dsolver = np.loadtxt(fsolver,delimiter=", ",dtype=complex,converters={1:parse_pair, 2:parse_pair})
dbingo1 = np.genfromtxt(fbingo1)
dbingo2 = np.genfromtxt(fbingo2)
dbingo3 = np.genfromtxt(fbingo3)
dbingo4 = np.genfromtxt(fbingo4)

t = dsolver[:,0]
x = dsolver[:,3]
tbingo1 = dbingo1[:,0]
xbingo1 = dbingo1[:,1]
tbingo2 = dbingo2[:,0]
xbingo2 = dbingo2[:,1]
tbingo3 = dbingo3[:,0]
xbingo3 = dbingo3[:,1]
tbingo4 = dbingo4[:,0]
xbingo4 = dbingo4[:,1]

fig,ax=plt.subplots(1,1)
plt.style.use('paper')
ax.set_prop_cycle(grays)
plt.loglog(t,x,lw=1.0,label='$\mathrm{RKWKB}$')
plt.loglog(tbingo1,xbingo1,'--',lw=1.0,label='$\mathrm{BINGO}$')
#plt.loglog(tbingo3[:4900],xbingo3[:4900],label='BINGO, rtol$=10^{-6}$',alpha=0.8,lw=1.3)
#plt.loglog(tbingo2[:4500],xbingo2[:4500],label='BINGO, rtol$=10^{-5}$',alpha=0.8,lw=1.3)
#plt.loglog(tbingo1[:4900],xbingo1[:4900],label='BINGO, rtol$=10^{-4}$',alpha=0.8,lw=1.3)
ax.text(0.0,1.01,'$\\times 10^{-9}$',transform=plt.gca().transAxes)
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e9)))
ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e9)))

plt.xlabel('$k$ [Mpc${}^{-1}$]')
plt.ylabel('$P_{\mathcal{R}}(k)$')
plt.xlim((1e-5,1e8))
plt.ylim((6e-10,4e-9))
plt.legend()
plt.tight_layout()
plt.savefig("plots/rkwkb-v-bingo-pps-labelled.pdf")
#plt.show()

# Timing comparison with BINGO 2 - like-to-like
# FIGURES 10 - 11
from cycler import cycler
grays = cycler(color=['0.00','0.40', '0.40', '0.60', '0.80', '0.90'])

fsolver1 = "test/ms/pps-bingo-start1e2-highk-lowk.txt"
fsolver2 = "test/ms/pps-bingo-start2e2-highk-lowk.txt"
fbingo1 = "test/ms/pps-bingo-comparison-ltl-moreacc-abs0.txt"
fbingo2 = "test/ms/pps-bingo-comparison-start2e2.txt"

dsolver1 = np.genfromtxt(fsolver1,delimiter=",")
dsolver2 = np.genfromtxt(fsolver2,delimiter=",")
dbingo1 = np.genfromtxt(fbingo1)
dbingo2 = np.genfromtxt(fbingo2)


k1 = dsolver1[:,0]
t1 = dsolver1[:,-1]
k2 = dsolver2[:,0]
t2 = dsolver2[:,-1]

kbingo1 = dbingo1[:,0]
tbingo1 = dbingo1[:,-1]
kbingo2 = dbingo2[:,0]
tbingo2 = dbingo2[:,-1]

# FIGURE 10
# tbingo/trkwkb plot
fig,ax=plt.subplots(1,1)
plt.style.use('paper')
ax.set_prop_cycle(grays)
plt.semilogx(k1,tbingo1/t1)
plt.xlabel('$k$ [Mpc${}^{-1}$]')
plt.ylabel('relative runtime, $t_{\mathrm{BINGO}}/t_{\mathrm{RKWKB}}$')
plt.xlim((1e-5,1e8))
#plt.show()
plt.savefig("plots/rkwkb-v-bingo-reltimes-semilog-labelled.pdf")


# abs runtimes
plt.semilogx(k1,t1,label='rkwkb, 100')
plt.semilogx(k2,t2,label='rkwkb, 200')
plt.semilogx(kbingo1, tbingo1,label='bingo 100')
plt.semilogx(kbingo2, tbingo2, label='bingo 200')

# relative runtime within trkwkb
# FIGURE 11
tpivot = t1[2502]
plt.style.use('paper')
plt.semilogx(k1,t1/tpivot)
plt.semilogx([1e-5,k1[2502]],[1,1], 'k:')
plt.semilogx([k1[2502],k1[2502]],[0,1],'k:')
plt.xlabel('$k$ [Mpc${}^{-1}$]')
plt.ylabel('relative runtime')
plt.xlim((1e-5,1e8))
plt.ylim((0.5,2.5))
#plt.show()
plt.savefig("plots/rkwkb-reltimes-labelled.pdf")

# RK single-k solution of BINGO eqs with NAG
f1 = "test/ms/bingo-singlek-k1e-5.txt"
f1ref = "test/ms/bingo-singlek-k1e-5-ref.txt"
f1wkb = "test/ms/rkwkb-single-k1e-5.txt"
d1 = np.genfromtxt(f1)
d1ref = np.genfromtxt(f1ref)
d1wkb = np.genfromtxt(f1wkb,dtype=complex,converters={1:parse_pair, 2:parse_pair})
n1 = d1[:,0]
n1ref = d1ref[:,0]
n1wkb = d1wkb[:,0]
rk1 = d1[:,1]
rk1ref = d1ref[:,1]
rk1wkb = d1wkb[:,1]
rk1steps = d1wkb[:,3]

f2 = "test/ms/bingo-singlek-k1e8.txt"
f2ref = "test/ms/bingo-singlek-k1e8-ref.txt"
f2wkb = "test/ms/rkwkb-single-k1e8.txt"
d2 = np.genfromtxt(f2)
d2ref = np.genfromtxt(f2ref)
d2wkb = np.genfromtxt(f2wkb,dtype=complex,converters={1:parse_pair, 2:parse_pair})
n2 = d2[:,0]
rk2 = d2[:,1]
rk2 = rk2*rk1[0]/rk2[0]
n2 = n2 - (n2[0] - n1[0]) 
n2ref = d2ref[:,0]
rk2ref = d2ref[:,1]
rk2ref = rk2ref*rk1ref[0]/rk2ref[0]
n2ref = n2ref - (n2ref[0] - n1ref[0])
n2wkb = d2wkb[:,0]
rk2wkb = d2wkb[:,1].real
rk2wkb = rk2wkb*rk1wkb[0].real/rk2wkb[0].real
n2wkb = n2wkb - (n2wkb[0] - n1wkb[0])

plt.style.use('default')
#fig,ax=plt.subplots(2,1,sharex=True)
plt.plot(n1ref,rk1ref,label="$k=10^{-5}$")
plt.plot(n2ref,rk2ref,label="$k=10^{8}$")
plt.legend()
plt.xlabel('N')
plt.ylabel('$\Re{\{ \mathcal{R}_k \}} $')
plt.savefig('plots/modes-highk-lowk.pdf')
#plt.show()

# FIGURE 12
plt.style.use('paper')
fig,ax = plt.subplots(2,1,sharex=True,sharey=True)#,figsize=(5.95,2.75))
ax[0].plot(n1ref,rk1ref)
ax[0].plot(n1,rk1,'.',color='red',label='RK step')
ax[0].text(-0.08,1.1, '(a)', transform=ax[0].transAxes)

ax[0].legend()
ax[1].plot(n1ref,rk1ref)
ax[1].plot(n1wkb[rk1steps==0],rk1wkb[rk1steps==0],'.',color='red',label='RK step')
ax[1].plot(n1wkb[rk1steps==1],rk1wkb[rk1steps==1],'^',color='green',ms=2.5,label='WKB step')
ax[1].legend(loc='upper right')
#fig.text(0.5,0.02,'$N$',ha='center',va='center')
fig.text(0.02,0.5,'$\Re{\{\mathcal{R}_k\}}$',ha='center',va='center',rotation='vertical')
ax[0].ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
ax[1].ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
ax[0].set_xlim((n1ref[0]-0.1,n1ref[-1]))
ax[1].set_xlim((n1ref[0]-0.1,n1ref[-1]))
ax[1].text(-0.08,1.1, '(b)', transform=ax[1].transAxes)

plt.xlabel("$N$")
#plt.ylabel("$\Re{\{\mathcal{R}_k\}}$")
plt.tight_layout()
plt.savefig('plots/rkwkb-bingo-singlek-labelled.pdf')#-elsevier-thin.pdf')
#plt.show()

# Kinetically dominated PPS in terms of N
# background check
def V(phi):
    m = 7.147378e-6
    return 0.5*m**2*phi**2

def dV(phi):
   m = 7.147378e-6
   return m**2*phi

fbg = "test/ms/kd-bg.txt"
dbg = np.genfromtxt(fbg,delimiter=", ")
N = dbg[:,0]
phi = dbg[:,1]
dphi = dbg[:,2]
ddphi = dbg[:,3]
H = dbg[:,4]
w = 1.0/(np.exp(N)*H)
g = -0.25*dphi**2 -(3.0-(0.5*dphi*dphi)) - (6.0-(dphi*dphi))*dV(phi)/(2.0*V(phi)*dphi) + 1.5
#plt.plot(N,phi)
#plt.semilogy(N,w)
plt.plot(N,2*ddphi/dphi)
plt.plot(N,2*(-(3.0-(0.5*dphi*dphi)) -
(6.0-(dphi*dphi))*dV(phi)/(2.0*V(phi)*dphi)),color='red')
#plt.semilogy(N,2e-4*w)
#plt.semilogy(N,2e4*w)
plt.show()

# KD PPS
# FIGURE 13
from matplotlib.ticker import FuncFormatter
f = "test/ms/pps-N-kd-sparsebg2.txt" #"test/ms/pps-testcomplexfn.txt" 
d = np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
p1 = d[:,3] # Bunch-Davies
p2 = d[:,4] # attempted kd
#plt.loglog(k,p1)
plt.style.use('paper')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("$k$ [Mpc${}^{-1}$]")
ax.set_ylabel("$P_{\mathcal{R}}(k)$")
ax.loglog(k[150:],p2[150:])
ax.text(0.0,1.01,'$\\times 10^{-9}$',transform=plt.gca().transAxes)
#ax.set_xlim((1e-3,1e1))
#ax.set_ylim((6e-10,2e-9))
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e9)))
ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e9)))
plt.tight_layout()
#plt.show()
plt.savefig("plots/example-pps-kd-N-labelled.pdf")


# Single-k solutions
fs1 = "test/ms/singlek-kd-bd-2a.txt"
fs2 = "test/ms/singlek-kd-bd-2a.txt"
fref = "test/ms/singlek-kd-bg-2a-ref.txt"
d1 = np.genfromtxt(fs1,dtype=complex,converters={1:parse_pair,2:parse_pair,4:parse_pair})
d2 = np.genfromtxt(fs2,dtype=complex,converters={1:parse_pair,2:parse_pair,4:parse_pair})
dref = np.genfromtxt(fref)
n1 = d1[:,0]
n2 = d2[:,0]
nref = dref[:,0]
rk1 = d1[:,1]
rk2 = d2[:,1]
rkref = dref[:,1]
wkb1 = d1[:,3]
wkb2 = d2[:,3]
plt.plot(n1[wkb1==1],(rk1[wkb1==1]),'x',color='green')
plt.plot(n1[wkb1==0],(rk1[wkb1==0]),'x',color='red')
plt.plot(nref,rkref,color='black')
#plt.plot(n1[wkb1==1],rk1[wkb2==1],'x',color='green')
#plt.plot(n1[wkb1==0],rk1[wkb2==0],'x',color='red')
#plt.plot(nref,dref[:,1],color='black')
plt.show()

# Kinetically dominated closed universe
# background check: horizon
f = "test/ms/kd-closed-bg-tip.txt"
f2 = "test/ms/kd-closed-bg-2.txt"
data = np.genfromtxt(f,delimiter=', ',comments='0, 0, 0')
data2 = np.genfromtxt(f2,delimiter=', ')
N = data[:,0]
phi = data[:,1]
dphi = data[:,2]
logo_k = data[:,3]
dlogo_k = data[:,4]
dE_E = data[:,7]
plt.style.use('paper')
#plt.plot(N,dlogo_k,label='dlogok')
#plt.plot(N,dE_E,label='ddphi')
#plt.plot(N,logo_k)
plt.plot(N,np.exp(0.5*(logo_k + np.log(data[:,5]))),color='red' ,label='$1/(aH)^2 $,mimic flat')
#plt.plot(N,np.exp(logo_k)*data[:,5],color='red' ,label='$1/(aH)^2 $,mimic flat')

#plt.plot(N, data[:,-2],color='red',label='$\gamma$,mimic flat')
#plt.semilogy(N,np.exp(-0.5*logo_k-N))
#plt.plot(N,data[:,5],label='$\omega$')
#plt.plot(N, logo_k,color='red',label='$K=+1$, closed')
#plt.plot(data2[:,0], data2[:,3],color='blue',label='$K=+1$, closed')
plt.legend()
plt.xlabel('$N$')
#plt.ylabel('$-2 \log{(aH)}$')
#plt.plot(N, phi, '.')
#plt.savefig('plots/closed-bg.pdf')
plt.show()

# For sake of debugging closed universe code
# Kinetically dominated PPS in terms of N
# background check
def V(phi):
    m = 5e-6#7.147378e-6
    return 0.5*m**2*phi**2

def dV(phi):
   m = 5e-6#7.147378e-6
   return m**2*phi

#k=0.0248
fbg = "test/ms/kd-bg.txt"
dbg = np.genfromtxt(fbg,delimiter=", ")
N = dbg[:,0]
phi = dbg[:,1]
dphi = dbg[:,2]
ddphi = dbg[:,3]
H = dbg[:,4]
w = 1.0/(np.exp(N)*H)
g = -0.25*dphi**2 -(3.0-(0.5*dphi*dphi)) - (6.0-(dphi*dphi))*dV(phi)/(2.0*V(phi)*dphi) + 1.5
#plt.plot(N,4.0-2.0*V(phi)/H**2,color='red',label='dlogok flat2')
#plt.plot(N,2/3.0*(dphi**2 - V(phi)/H**2),color='red',label='dlogok flat')
plt.plot(N,np.log(1.0/(H*np.exp(N))**2),label='$1/(aH)^2$ flat')
#plt.plot(N,(-(3.0-(0.5*dphi*dphi))*dphi -
#(6.0-(dphi*dphi))*dV(phi)/(2.0*V(phi))),color='red',label='ddphi flat')
#plt.plot(N,phi)
plt.plot(N,g,label='$\gamma$,flat')
#plt.plot(N,np.log(w/k),color='red')
#plt.plot(N,2*ddphi/dphi,label='2$\ddot{\phi}/{\phi}$ explicit')
#plt.plot(N,2*(-(3.0-(0.5*dphi*dphi)) -
#(6.0-(dphi*dphi))*dV(phi)/(2.0*V(phi)*dphi)),color='red',label='2$\ddot{\phi}/{\phi}$ implicit')
#plt.semilogy(N,2e-4*w)
#plt.semilogy(N,2e4*w)
plt.legend()
plt.show()

# PPS
f = "test/ms/pps-test.txt"
d = np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
p1 = d[:,3] # BD 
p2 = d[:,4] # KD
#plt.loglog(k,p1)
plt.style.use('paper')
plt.xlabel("$k$/Mpc${}^{-1}$")
plt.ylabel("$P_{\mathcal{R}}(k)$")
#plt.xlim((8e-1,1e1))
#plt.ylim((1e-3,1e2))
plt.loglog(k,p1)
#fig.text(0.5,0.02,'$N$',ha='center',va='center')
#plt.show()
plt.savefig("plots/lowk-exponential.pdf")

# single-k solution
f = "test/ms/kd-mimicflat2.txt"
d = np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
rk = d[:,1]
wkb = d[:,3]
plt.style.use('default')
plt.plot(k[wkb==1],rk[wkb==1],'^',color='green',ms=2.5,label='wkb step')
plt.plot(k[wkb==0],rk[wkb==0],'x',color='red',label='rk step')
#plt.plot(k,rk)

plt.legend()
plt.show()

# Effect of decreasing curvature on a curved universe
# PPS
from cycler import cycler
#grays = cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728','#9467bd',
#'#8c564b', '#e377c2', '#7f7f7f','#bcbd22', '#17becf'])
cmap = plt.cm.tab20
grays = cycler('color',cmap(np.linspace(0,1,20)))
nospectra=6
atoday = 4.3e4
fs = ["test/ms/pps-kd-closed-{}intcorr.txt".format(i) for i in range(1,nospectra+1)]
okis = np.logspace(-3,-3+nospectra-1,nospectra) 
ds = [np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair}) for f in fs]
ks = [d[:,0]*atoday for d in ds]
ps = [d[:,3] for d in ds]  # BD 
plt.style.use('paper')
fig,ax=plt.subplots(1,1,figsize=(5.95,4.0))
ax.set_prop_cycle(grays)
ax2=ax.twiny()
#plt.xlim((1e-5,1e1))
#plt.ylim((6e-10,2e-9))
pivotamp = ps[0][-1]
for i in range(nospectra):
    if(i<nospectra/2):
        ax.loglog(ks[i],ps[i]*pivotamp/ps[i][-1],label='$\Omega_k^i={}$'.format(okis[i]),lw=1.0)
    else:
        ax.loglog(ks[i],ps[i]*pivotamp/ps[i][-1],label='$\Omega_k^i={}$'.format(okis[i]),lw=1.0)
axxs=np.logspace(0,8,9)
#print(axxs)
ax2xs=axxs/atoday
#print(ax2xs)
ax2.set_xticks(axxs)
ax2.set_xticklabels(ax2xs)
ax2.set_xlabel("$k\\big/\\big(\\big(\\frac{h}{\mathrm{0.7}}\\big) \\big(\\frac{\Omega_{k,\mathrm{0}}}{\mathrm{0.01}}\\big)$ Mpc${}^{-1}$\\big)\n ")
ax.set_ylabel("$P_{\mathcal{R}}(k)/m^2$")
ax.set_xlabel("comoving $k$")
ax2.loglog(ax2xs,np.ones_like(ax2xs),alpha=0) 
ax.legend()
plt.tight_layout()
#plt.show()
#plt.savefig("plots/closed-spectra-log-paper-colours.pdf")

# Schrodinger equation
# FIGURE 8 
import scipy.special as sp
K=2
ns=[2,20,40,60,80,100]
m=1
Es=[np.sqrt(K/m)*(n-0.5) for n in ns]
gam = np.sqrt(m*K)
fls = ["test/schrodinger/schrodinger-plot-{}l.txt".format(n) for n in ns]
frs = ["test/schrodinger/schrodinger-plot-{}r.txt".format(n) for n in ns]
dls = [np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair}) for
f in fls]
drs = [np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair}) for
f in frs]
tls = [np.linspace(np.real(d[0,0]),np.real(d[-1,0]),1000) for d in dls]
trs = [np.linspace(np.real(d[0,0]),np.real(d[-1,0]),1000) for d in drs]
As = [(gam/np.pi)**(1/4.0)*1.0/np.sqrt(2.0**(n-1)*math.factorial(n-1)) for n in ns]
scale = 10
lanalytics = [scale*A*sp.eval_hermite(n-1,np.sqrt(gam)*t)*np.exp(-0.5*gam*t**2) for A,t,n in zip(As,tls,ns)]
ranalytics = [scale*A*sp.eval_hermite(n-1,np.sqrt(gam)*t)*np.exp(-0.5*gam*t**2) for A,t,n in zip(As,trs,ns)]
vl = 0.5*K*tls[-1]**2 
vr = 0.5*K*trs[-1]**2 
plt.style.use('paper')
plt.figure(figsize=(5.95,4.0))
plt.plot(tls[-1],vl,'--',label='$V(x)$')
plt.plot(trs[-1],vr,'--')

for i in range(len(ns)):
    d = dls[i]
    xl = d[:,0]
    rkl = d[:,1]
    wkbl = d[:,3]
    E = Es[i]
    d = drs[i]
    xr = d[:,0]
    rkr = d[:,1]
    wkbr = d[:,3]

    if(i<len(ns)-1):
        plt.plot(tls[i],E + (lanalytics[i]))
        plt.plot(trs[i],E + (ranalytics[i]))
        plt.plot(xl[wkbl==1],E + scale*((rkl[wkbl==1])),'^',ms=3,color='green')
        plt.plot(xl[wkbl==0],E + scale*((rkl[wkbl==0])),'.',color='red')
        plt.plot(xr[wkbr==1],E + scale*((rkr[wkbr==1])),'^',ms=3,color='green')
        plt.plot(xr[wkbr==0],E + scale*((rkr[wkbr==0])),'.',color='red')
    else:
        plt.plot(tls[i],E + (lanalytics[i]))
        plt.plot(trs[i],E + (ranalytics[i]))#,label='true solution')
        plt.plot(xl[wkbl==1],E + scale*((rkl[wkbl==1])),'^',ms=3,color='green')
        plt.plot(xl[wkbl==0],E + scale*((rkl[wkbl==0])),'.',color='red')
        plt.plot(xr[wkbr==0],E + scale*((rkr[wkbr==0])),'.',color='red',label='RK step')
        plt.plot(xr[wkbr==1],E + scale*((rkr[wkbr==1])),'^',color='green',ms=3,label='WKB step')
        plt.text(xr[0]+1.0, E, 'n='+str(ns[i]))
plt.xlabel('$x$')
plt.ylabel('$E_n + 10\Psi_n(x)$')
plt.legend(loc='lower left')
#plt.plot(ts,E*np.ones(len(ts)))
plt.xlim((-15,18))
plt.ylim((-15,160))
plt.tight_layout()
plt.savefig('plots/schrodinger-shooting.pdf')
#plt.show()

# Schrodinger equation with outside-of-potential tails
K=2
n=20
m=1
E=np.sqrt(K/m)*(n-0.5)
gam = np.sqrt(m*K)
nofs=2
boundary=np.sqrt((n-0.5)*2/gam)+5.0
fs = ["test/schrodinger/schrodinger-{}.txt".format(f) for f in ['left','right']]
ds = [np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair}) for f in fs]
ts = np.linspace(-boundary,boundary,1000)
A = (gam/np.pi)**(1/4.0)*1.0/np.sqrt(2.0**(n-1)*math.factorial(n-1))
analytic = [A*sp.eval_hermite(n-1,np.sqrt(gam)*t)*np.exp(-0.5*gam*t**2) for t in ts]
v = 0.5*K*ts**2 
plt.style.use('paper')
#plt.plot(ts,v,'--',label='$V(x)$')
for i in range(nofs):
    d = ds[i]
    x = d[:,0]
    rk = d[:,1]
    wkb = d[:,3]
    plt.plot(x[wkb==1],rk[wkb==1],'.',color='green')
    plt.plot(x[wkb==0],rk[wkb==0],'.',color='red')
#plt.ylim((-2,2))
plt.plot(ts, analytic)
plt.xlim(-boundary,boundary)
plt.savefig("plots/schrodingershooting.pdf")
#plt.show()

# PLOT FOR WILL
# Determining Ni as a function of O_ki to maintain a constant Ntot=60
from scipy.interpolate import interp1d as interp 
f = "test/ms/closed-bg-constn.txt"
data = np.genfromtxt(f,dtype=complex,delimiter=', ')
ok = data[:,0]
ni = data[:,1]
plt.style.use('paper')
plt.ylabel('$N_i$')
plt.xlabel('$\Omega_k^i$')
plt.semilogx(ok,ni,label='$N_{\mathrm{tot}}=60$')
plt.legend()
#plt.savefig('plots/ni-v-oki.pdf')
plt.show() 

# FIGURE 14
# Improved plot showing closed universe spectra with varying curvature
atoday=4.3e4
okis=np.logspace(-3,2.17,300)
labels=np.linspace(0,300,300)
plt.style.use('paper')
fig,axes=plt.subplots(3,2,figsize=(5.90,7.7),sharex=True,sharey='row')
fs = []
fints = []
for num in [0,58,116,174,231,290]:
    print(okis[num])
    fs.append("test/ms/pps-kd-closed-kcts-corr{}.txt".format(num))
    fints.append("test/ms/pps-kd-closed-kint-corr{}.txt".format(num))

ds = [np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair}) for f in fs]
dints = [np.genfromtxt(fint,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair}) for fint in fints]

fig.text(0.5,1.01,"$k \\big[\\big(\\frac{h}{\mathrm{0.7}}\\big) \\big(\\frac{\Omega_{k,\mathrm{0}}}{\mathrm{0.01}}\\big)$ Mpc${}^{-1}$\\big]$ $",ha='center',va='center')
fig.text(0.01,0.5,"$P_{\mathcal{R}}(k)/m^2$",ha='center',va='center',rotation='vertical')
fig.text(0.5,0.00,"comoving $k$",ha='center',va='center')

okis = [1e-3,1e-2,1e-1,1e0,1e1,1e2]
labels = ['(a)','(b)','(c)','(d)','(e)','(f)']
for i,ax,d,dint in zip(range(6),np.ravel(axes),ds,dints):
    ax2=ax.twiny()
    axxs=np.logspace(0,3,4)
    ax2xs=axxs/atoday
    ax2.set_xticks(axxs)
    if(i==1 or i==0):
        ax2.set_xticklabels(ax2xs)
    else:
        ax2.set_xticklabels(['']*len(ax2xs))
    ax2.loglog(ax2xs,np.ones_like(ax2xs),alpha=0)
    if(i!=0 and i!=1):
        ax2.set_xticklabels(['']*len(ax2xs))
    ax.loglog(d[2:,0],d[2:,3],color='red')
    ax.loglog(dint[2:,0],dint[2:,3])
    ax.loglog(dint[2:50,0],dint[2:50,3],'.',ms=2.5)
    ax.set_xlim((5e-1,2e4))
    if(i==0 or i==1):
        ax.set_ylim((1e-2,1e3))
    elif(i==4 or i==5):
        ax.set_ylim((1e0,1e2))
    else:
        ax.set_ylim((1e0,1e3))
    ax.text(0.70,0.1,'$\Omega_k^i=$'+str(okis[i]),transform=ax.transAxes)
    ax.text(0.05,0.9, labels[i], transform=ax.transAxes)

plt.subplots_adjust(wspace=0.1,hspace=0.1)
plt.tight_layout()
plt.savefig("plots/closed-spectra-table-k3-labelled.pdf")
#plt.show()

