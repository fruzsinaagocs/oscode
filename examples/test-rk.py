import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import re
import math

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
def parse_pair(s):
    return complex(*map(float, pair.match(s.decode('utf-8')).groups()))

# FIGURE 1, Airy equation
f1 = "example.txt"
data = np.genfromtxt(f1,dtype=complex,converters={1:parse_pair,2:parse_pair},delimiter=";",missing_values="",filling_values=-1.0+0j,comments='#')
times = np.logspace(math.log10(1.0),math.log10(60.0),5000)

wkbs = data[:,-1]
twkb = data[:,0][wkbs==1]
trk = data[:,0][wkbs==0]
tdense = data[:,0][wkbs==-1]
xwkb = data[:,1][wkbs==1]
xrk = data[:,1][wkbs==0]
xdense = data[:,1][wkbs==-1]

analytic = np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in times])

def ana(t):
    return np.array([sp.airy(-ti)[0] +1j*sp.airy(-ti)[2] for ti in t])

#plt.style.use('dense')
fig,ax=plt.subplots(1,2)
ax[0].plot(times,analytic,label='analytic solution',color='black',lw=0.7)
ax[0].plot(trk,xrk,'.',label='oscode',color="C1",ms=7.0)
#ax[0].plot(twkb,xwkb,'.',color="C1",ms=7.0)
ax[0].plot(tdense,xdense,'.',color="C0",label='oscode dense output')
ax[0].set_ylabel('$x$')
ax[0].set_xlabel('$t$')
ax[0].legend()

ax[1].plot(trk,abs((xrk-ana(trk))/ana(trk)),'.',label='rk step')
ax[1].plot(tdense,abs((xdense-ana(tdense))/ana(tdense)),'.',label='rk dense')
ax[1].legend()

plt.tight_layout()
plt.savefig('rk-dense-output-rtol1e-4.pdf')
#plt.show()


