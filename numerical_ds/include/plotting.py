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
dsolsrk = datark[:,5]
errs = np.abs(sols - x)/np.abs(sols)
derrs = np.abs(dsols - dx)/np.abs(dsols)
errsrk = np.abs(solsrk - xrk)/np.abs(solsrk)
derrsrk = np.abs(dsolsrk - dxrk)/np.abs(dsolsrk)

# Example solution
# FIGURE 1 

plt.style.use('paper')
plt.semilogx(times,analytic, color='black',label='true solution',lw=0.5)
plt.semilogx(t[wkb==1],x[wkb==1],'x',color='green',label='WKB step')
plt.semilogx(t[wkb==0],x[wkb==0],'x',color='red', label='RK step')
plt.xlim((0,40.0))
plt.ylim((-0.5, 0.7))
plt.xlabel("$t$")
plt.ylabel("$\Re\{x\}$")
plt.legend()
plt.show()
#plt.savefig("plots/airy-example-x-paperstyle.pdf")

# Error progression to late times in the above example
# FIGURE 2

plt.style.use('paper')
plt.loglog(t, errs, '-', color='black')
plt.loglog(t[wkb==1], errs[wkb==1],'x',color='green',label='WKB step')
plt.loglog(t[wkb==0], errs[wkb==0],'x',color='red',label='RK step')
plt.loglog(trk[wkbrk==0], errsrk[wkbrk==0],'x',color='red')
plt.loglog(trk, errsrk, '-', color='black')
plt.ylim((1e-6,1e-2))
plt.xlabel("$t$")
plt.ylabel("relative error, $\\frac{|\Delta x|}{|x|}$")
#plt.axvline(x=3.79)
#plt.text(3.79/(7.0/3.79), 2e-4, 'RK', verticalalignment='center',
#horizontalalignment='right')
#plt.text(7.0, 2e-4, 'WKB', verticalalignment='center',
#horizontalalignment='left')
plt.legend()
#plt.show()
plt.savefig('plots/airy-example-err-paperstyle.pdf')


# Burst equation with w(t), g(t) given analytically
# Single solution at low n 
# FIGURE 3 

f1 = 'plots/burstn40.txt'#'test/burst/solcheck3-optn.txt'#d4w1less_n40.txt' #'plots/burstn40.txt'
n = 40 
data= np.loadtxt(f1,dtype=complex,converters={1:parse_pair, 2:parse_pair})
times = np.linspace(-2*n,2*n,10000)
analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])
t = data[:,0]
x = data[:,1]
dx = data[:,2]
wkb = data[:,3]
sols = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in t])

plt.figure()
plt.style.use('paper')
plt.plot(times,analytic, color='black',label='true solution')
plt.plot(t[wkb==1],x[wkb==1],'x',color='green',label='WKB step')
plt.plot(t[wkb==0],x[wkb==0],'x',color='red',label='RK step')
plt.xlim((-40,40))
plt.ylim((-60,60))
plt.xlabel('$t$')
plt.ylabel('$\Re{\{x(t)\}}$')
plt.legend()
plt.savefig("plots/burstn40_x_paperstyle.pdf")
#plt.show()

# Error progression in this example
# FIGURE 4
errs = np.abs(sols - x)/np.abs(sols)
plt.figure()
plt.style.use('paper')
plt.semilogy(t,errs,color='black')
plt.semilogy(t[wkb==1],errs[wkb==1],'x',color='green',label='WKB step')
plt.semilogy(t[wkb==0],errs[wkb==0],'x',color='red',label='RK step')

plt.ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
plt.xlabel('$t$')
plt.ylim(1e-7,1e-1)
plt.xlim(-60,60)
plt.legend()
plt.savefig("plots/burstn40_err_paperstyle.pdf")
#plt.show()

# Large n example: number of oscillations traversed
# FIGURE 5
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
plt.style.use("paper")
plt.semilogy(ts,oscs,color='black')
plt.semilogy(ts[wkbs==1],oscs[wkbs==1],'.',label='WKB step',color='green')
plt.semilogy(ts[wkbs==0],oscs[wkbs==0],'.',label='RK step',color='red')
plt.xlim((-n,n))
plt.xlabel('$t$')
plt.ylim((1e-2,1e4))
plt.ylabel('oscillations traversed')
plt.legend()
#plt.show()
plt.savefig('plots/burstn1e5_oscs.pdf')

# Changing rtol at n=1e5, effect on relative error progression 
# FIGURE 6
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
plt.style.use('paper')
plt.semilogy(t1,errs1,'-',color='black',label='rtol=$10^{-4}$')
plt.semilogy(t2,errs2,'-.',color='black',label='rtol=$10^{-5}$')
plt.semilogy(t3,errs3,'--',color='black',label='rtol=$10^{-6}$')
plt.ylabel('relative error, $\\frac{|\Delta x|}{|x|}$')
plt.xlabel('$t$')
plt.xlim((-2*n,2*n))
plt.ylim((1e-9, 1e-1))
plt.legend()
#plt.show()
plt.savefig("plots/burstn1e5_rtols.pdf")

# Effect of changing rtol and n on runtime
# FIGURE 7
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
plt.savefig("plots/bursttiming.pdf")

# Effect of changing n, rtol on number of steps
# FIGURE 9
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
plt.loglog(10**ns,wkbs1,'-',color='green',label='WKB steps')
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
# FIGURE 8
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
plt.plot(t[wkb==1],x[wkb==1],'x',color='green',label='WKB step')
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

# FIGURE 10
from cycler import cycler
grays = cycler(color=['0.00','0.60', '0.40', '0.60', '0.80', '0.90'])
#default = cycler(color=['#1f77b4','#ff7f0e', '#d62728', '#9467bd', '#8c564b',
#'#e377c2'])

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

plt.xlabel('$k$/Mpc${}^{-1}$')
plt.ylabel('$P_{\mathcal{R}}(k)$')
plt.xlim((1e-5,1e8))
plt.ylim((6e-10,4e-9))
plt.legend()
plt.savefig("plots/rkwkb-v-bingo-pps.pdf")
#plt.show()

# Timing comparison with BINGO 2 - like-to-like
# FIGURES 11 - 12
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

# FIGURE 11
# tbingo/trkwkb plot
fig,ax=plt.subplots(1,1)
plt.style.use('paper')
ax.set_prop_cycle(grays)
plt.semilogx(k1,tbingo1/t1)
plt.xlabel('$k$')
plt.ylabel('relative runtime, $t_{\mathrm{BINGO}}/t_{\mathrm{RKWKB}}$')
plt.xlim((1e-5,1e8))
#plt.show()
plt.savefig("plots/rkwkb-v-bingo-reltimes-semilog.pdf")


# abs runtimes
plt.semilogx(k1,t1,label='rkwkb, 100')
plt.semilogx(k2,t2,label='rkwkb, 200')
plt.semilogx(kbingo1, tbingo1,label='bingo 100')
plt.semilogx(kbingo2, tbingo2, label='bingo 200')

# relative runtime within trkwkb
# FIGURE 12
tpivot = t1[2502]
plt.style.use('paper')
plt.loglog(k1,t1/tpivot)
plt.loglog([1e-5,k1[2502]],[1,1], 'k:')
plt.loglog([k1[2502],k1[2502]],[0,1],'k:')
plt.xlabel('$k$')
plt.ylabel('relative runtime')
plt.xlim((1e-5,1e8))
#plt.show()
plt.savefig("plots/rkwkb-reltimes.pdf")

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

# FIGURE 13
plt.style.use('paper')
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(n1ref,rk1ref)
ax[0].plot(n1,rk1,'x',color='red',label='RK step')
ax[0].legend()
ax[1].plot(n1ref,rk1ref)
ax[1].plot(n1wkb[rk1steps==0],rk1wkb[rk1steps==0],'x',color='red',label='RK step')
ax[1].plot(n1wkb[rk1steps==1],rk1wkb[rk1steps==1],'x',color='green',label='WKB step')
ax[1].legend(loc='upper right')
fig.text(0.5,0.02,'$N$',ha='center',va='center')
fig.text(0.06,0.5,'$\Re{\{\mathcal{R}_k\}}$',ha='center',va='center',rotation='vertical')
ax[0].set_xlim((n1ref[0],n1ref[-1]))
ax[1].set_xlim((n1ref[0],n1ref[-1]))
plt.savefig('plots/rkwkb-bingo-singlek.pdf')
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

# PPS
# FIGURE 14
#plt.figure(figsize=(2.60,2.375))
f = "test/ms/pps-testcomplexfn.txt" #pps-N-kd-sparsebg2.txt" #sparsebg2.txt
d = np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
p1 = d[:,3] # Bunch-Davies
p2 = d[:,4] # attempted kd
#plt.loglog(k,p1)
plt.style.use('paper')
plt.xlabel("$k$/Mpc${}^{-1}$")
plt.ylabel("$P_{\mathcal{R}}(k)$")
#plt.xlim((1e-3,1e5))
#plt.ylim((6e-10,2e-9))
plt.loglog(k,p1)
#fig.text(0.5,0.02,'$N$',ha='center',va='center')
plt.show()
#plt.savefig("plots/example-pps-kd-N.pdf")

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
#plt.plot(N,phi)
plt.plot(N,np.exp(0.5*(logo_k + np.log(data[:,5]))),color='red' ,label='$1/(aH)^2 $,mimic flat')
#plt.plot(N, data[:,-2],color='red',label='$\gamma$,mimic flat')
#plt.semilogy(N,np.exp(-0.5*logo_k-N))
plt.plot(N,data[:,5],label='$\omega$')
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
f = "test/ms/pps-kd-closed-highk12.txt"
d = np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
p1 = d[:,3] # BD 
p2 = d[:,4] # KD
#plt.loglog(k,p1)
plt.style.use('paper')
plt.xlabel("$k$/Mpc${}^{-1}$")
plt.ylabel("$P_{\mathcal{R}}(k)$")
#plt.xlim((1e-5,1e1))
#plt.ylim((6e-10,2e-9))
plt.loglog(k,p1)
#fig.text(0.5,0.02,'$N$',ha='center',va='center')
plt.show()
#plt.savefig("plots/example-pps-kd-N.pdf")

# single-k solution
f = "test/ms/kd-mimicflat2.txt"
d = np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair})
k = d[:,0]
rk = d[:,1]
wkb = d[:,3]
plt.style.use('default')
plt.plot(k[wkb==1],rk[wkb==1],'x',color='green',label='wkb step')
plt.plot(k[wkb==0],rk[wkb==0],'x',color='red',label='rk step')
#plt.plot(k,rk)

plt.legend()
plt.show()

# Effect of decreasing curvature on a curved universe
# PPS
# FIGURE 15
from cycler import cycler
grays = cycler(color=['0.80','0.50', '0.30', '0.80', '0.50', '0.30'])
nospectra=6
atoday = 4.3e4
fs = ["test/ms/pps-kd-closed-{}intcorr.txt".format(i) for i in range(1,nospectra+1)]
okis = np.logspace(-3,-3+nospectra-1,nospectra) 
ds = [np.genfromtxt(f,delimiter=", ",dtype=complex,converters={1:parse_pair,2:parse_pair}) for f in fs]
ks = [d[:,0]*atoday for d in ds]
ps = [d[:,3] for d in ds]  # BD 
plt.style.use('paper')
fig,ax=plt.subplots(1,1)
ax.set_prop_cycle(grays)
ax2=ax.twiny()
#plt.xlim((1e-5,1e1))
#plt.ylim((6e-10,2e-9))
pivotamp = ps[0][-1]
for i in range(nospectra):
    if(i<nospectra/2):
        ax.loglog(ks[i],ps[i]*pivotamp/ps[i][-1],label='$\Omega_k^i={}$'.format(okis[i]),lw=1.0)
    else:
        ax.loglog(ks[i],ps[i]*pivotamp/ps[i][-1],'--',label='$\Omega_k^i={}$'.format(okis[i]),lw=1.0)
axxs=np.logspace(0,8,9)
#print(axxs)
ax2xs=axxs/atoday
#print(ax2xs)
ax2.set_xticks(axxs)
ax2.set_xticklabels(ax2xs)
ax2.set_xlabel("$k\\big/\\big(\\big(\\frac{h}{\mathrm{0.7}}\\big) \\big(\\frac{\Omega_{k,\mathrm{0}}}{\mathrm{0.01}}\\big)$ Mpc${}^{-1}$\\big)$ $")
ax.set_ylabel("$P_{\mathcal{R}}(k)/m^2$")
ax.set_xlabel("comoving $k$")
ax2.loglog(ax2xs,np.ones_like(ax2xs),alpha=0) 
ax.legend()
#plt.show()
plt.savefig("plots/closed-spectra-log-paper.pdf")

# Schrodinger equation
# FIGURE 16
import scipy.special as sp
K=2
ns=[2,20,40,60,80,100]
m=1
Es=[np.sqrt(K/m)*(n-0.5) for n in ns]
gam = np.sqrt(m*K)
fs = ["test/schrodinger/schrodinger-test{}.txt".format(n) for n in ns]
ds = [np.genfromtxt(f,dtype=complex,converters={1:parse_pair,2:parse_pair}) for f in fs]
ts = [np.linspace(np.real(d[0,0]),np.real(d[-1,0]),1000) for d in ds]
As = [(gam/np.pi)**(1/4.0)*1.0/np.sqrt(2.0**(n-1)*math.factorial(n-1)) for n in ns]
scale = 10
analytics = [scale*A*sp.eval_hermite(n-1,np.sqrt(gam)*t)*np.exp(-0.5*gam*t**2) for A,t,n in zip(As,ts,ns)]
v = 0.5*K*ts[-1]**2 
plt.style.use('paper')
plt.plot(ts[-1],v,'--',label='$V(x)$')
for i in range(len(ns)):
    d = ds[i]
    x = d[:,0]
    rk = d[:,1]
    wkb = d[:,3]
    E = Es[i]
    if(i<len(ns)-1):
        plt.plot(x[wkb==1],E + scale*((rk[wkb==1])),'.',color='green')
        plt.plot(x[wkb==0],E + scale*((rk[wkb==0])),'.',color='red')
        plt.plot(ts[i],E + (analytics[i]))
    else:
        plt.plot(x[wkb==0],E + scale*((rk[wkb==0])),'.',color='red',label='RK step')
        plt.plot(x[wkb==1],E + scale*((rk[wkb==1])),'.',color='green',label='WKB step')
        plt.plot(ts[i],E + (analytics[i]))#,label='true solution')
    plt.text(x[-1]+1.0, E, 'n='+str(ns[i]))
plt.xlabel('$x$')
plt.ylabel('$E_n + 10\Psi_n(x)$')
plt.legend(loc='lower left')
#plt.plot(ts,E*np.ones(len(ts)))
plt.xlim((-25,21))
#plt.ylim((E-1,E+1))
#plt.savefig('plots/schrodinger.pdf')
plt.show()

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


# Determining Ni
from scipy.interpolate import interp1d as interp 
x = np.array([1e-3,1e-2,1e-1,1e0,1e1,1e2])
y = np.array([1.39, 0.25, -0.89, -1.76,-2.10,-2.215])
xnew = np.linspace(-3,2,300)
spl=interp(np.log10(x),y,kind=2)
ysmooth=spl(xnew)
print(ysmooth)
plt.plot(np.log10(x),y,'.')
plt.plot(xnew,ysmooth)
plt.show() 

