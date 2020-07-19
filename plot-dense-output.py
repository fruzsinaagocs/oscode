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

# FIGURE 1, Airy equation
f1 = "test/dense_output_5.txt"
data = np.genfromtxt(f1,dtype=complex,converters={1:parse_pair},delimiter=";",missing_values="",filling_values=-1.0+0j)
times = np.logspace(math.log10(1.0),math.log10(60.0),5000)

wkbs = data[:,-1]
twkb = data[:,0][wkbs==1]
trk = data[:,0][wkbs==0]
tdense = data[:,0][wkbs==-1]
xwkb = data[:,1][wkbs==1]
xrk = data[:,1][wkbs==0]
xdense = data[:,1][wkbs==-1]

analytic = np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in times])

plt.style.use('dense')
plt.figure()
plt.plot(times,analytic,label='analytic solution',color='black',lw=0.7)
plt.plot(trk,xrk,'.',label='oscode',color="C1",ms=7.0)
plt.plot(twkb,xwkb,'.',color="C1",ms=7.0)
plt.plot(tdense,xdense,color="C1",label='oscode dense output')
plt.ylabel('$x$')
plt.xlabel('$t$')

plt.legend()
plt.tight_layout()
plt.savefig("plots/dense-output-airy.pdf")
#plt.show()


# FIGURE 2, burst equation
f1 = "test/dense_output_6.txt"
n=40.0
data = np.genfromtxt(f1,dtype=complex,converters={1:parse_pair},delimiter=";",missing_values="",filling_values=-1.0+0j)
times = np.linspace(-2*n,2*n,8000)

analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])

wkbs = data[:,-1]
twkb = data[:,0][wkbs==1]
trk = data[:,0][wkbs==0]
tdense = data[:,0][wkbs==-1]
xwkb = data[:,1][wkbs==1]
xrk = data[:,1][wkbs==0]
xdense = data[:,1][wkbs==-1]

plt.style.use('dense')
plt.figure()
plt.plot(times,analytic,label='analytic solution',color='black',lw=0.7)
plt.plot(trk,xrk,'.',label='oscode',color="C1",ms=7.0)
plt.plot(twkb,xwkb,'.',color="C1",ms=7.0)
plt.plot(tdense,xdense,color="C1",label='oscode dense output')
plt.xlabel('$t$')
plt.ylabel('$x$')

plt.legend()
plt.xlim((-n/5,n/5))
plt.ylim((-15.0,20.0))
plt.tight_layout()
plt.savefig("plots/dense-output-burst-n40.pdf")
#plt.show()

# FIGURE 3, burst equation, high n
f1 = "test/dense_output_7.txt"
n=1000.0
data = np.genfromtxt(f1,dtype=complex,converters={1:parse_pair},delimiter=";",missing_values="",filling_values=-1.0+0j)
times = np.linspace(-0.6,0.65,8000)

analytic = np.array([100*np.sqrt(1+ti**2)/n * (1j*np.sin(n * np.arctan(ti)) + np.cos(n * np.arctan(ti))) for ti in times ])

wkbs = data[:,-1]
twkb = data[:,0][wkbs==1]
trk = data[:,0][wkbs==0]
tdense = data[:,0][wkbs==-1]
xwkb = data[:,1][wkbs==1]
xrk = data[:,1][wkbs==0]
xdense = data[:,1][wkbs==-1]

plt.style.use('dense')
plt.figure(figsize=(6.0,3.0))
plt.plot(times,analytic,label='analytic solution',color='black',lw=0.7)
plt.plot(trk,xrk,'.',label='oscode',color="C1",ms=7.0)
plt.plot(twkb,xwkb,'.',color="C1",ms=7.0)
plt.plot(tdense,xdense,color="C1",label='oscode dense output')

plt.legend()
#plt.xlim((-n/5,n/5))
#plt.ylim((-15.0,20.0))
plt.tight_layout()
#plt.savefig("plots/dense-output-burst-n40.pdf")
plt.show()


# FIGURE 4, dense output residuals from Airy equation
f1 = "test/dense_output_res_airy_4.txt"
data = np.genfromtxt(f1,dtype=complex,converters={1:parse_pair},delimiter=";",missing_values="",filling_values=-1.0+0j)
times = np.logspace(math.log10(1.0),math.log10(40.0),5000)

wkbs = data[:,-1]
twkb = data[:,0][wkbs==1]
trk = data[:,0][wkbs==0]
tdense = data[:,0][wkbs==-1]
xwkb = data[:,1][wkbs==1]
xrk = data[:,1][wkbs==0]
xdense = data[:,1][wkbs==-1]

anawkb = np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in twkb])
anadense = np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in tdense])
analytic = np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in times])

plt.style.use('dense')
fig,ax = plt.subplots(2,1,sharex=True)

#ax[0].plot(times,analytic,label='analytic solution',color='black',lw=0.7)
#ax[0].plot(trk,xrk,'.',label='oscode',color="C1",ms=7.0)
#ax[0].plot(tdense,xdense,color="C1",label='oscode dense output')
#ax[0].plot(twkb,xwkb,'.',color="C0",ms=7.0)
#ax[0].legend()

ax[0].plot(tdense,abs(xdense-anadense),color='black',lw=0.7,label='dense output')
ax[0].plot(twkb,abs(xwkb-anawkb),'.',ms=7.0,color='C1',label='pyoscode natural step')
ax[0].set_ylabel('$|\Delta x|$')
#ax[0].set_xlabel('$t$')
#ax[0].legend()

ax[1].plot(tdense,abs((xdense-anadense)/anadense),lw=0.7,color='black',label='dense output')
ax[1].plot(twkb,abs((xwkb-anawkb)/anawkb),'.',ms=7.0,color='C1',label='pyoscode natural step')
ax[1].set_ylabel('$|\Delta x/x|$')
ax[1].set_xlabel('$t$')
ax[1].legend()

plt.tight_layout()
#plt.savefig("plots/dense-output-airy-res.pdf")
plt.show()



