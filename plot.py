#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
from scipy.special import airy
from ast import literal_eval as le

def plot_airy(nag_input, rkwkb_input):

    nag_data = np.loadtxt(nag_input)
    rkwkb_data = np.genfromtxt(rkwkb_input,dtype=str)
    nag_times = nag_data[:,0][::10]
    rkwkb_times = rkwkb_data[:,0]
    rkwkb_times = np.array([le(datum) for datum in rkwkb_times])
    rkwkb_methods = rkwkb_data[:,-1]
    nag_ai = nag_data[:,-2][::10]
    rkwkb_ai = rkwkb_data[:,1]
    rkwkb_ai = np.array([le(datum)[0]+1j*le(datum)[1] for datum in rkwkb_ai])
    nag_ana = airy(-nag_times)[0]
    rkwkb_ana = np.array([airy(-time)[0] +1j*airy(-time)[2] for time in rkwkb_times])
    rkwkb_times_wkb = np.array([x for i,x in enumerate(rkwkb_times) if rkwkb_methods[i]=='1'])
    rkwkb_times_rk = np.array([x for i,x in enumerate(rkwkb_times) if rkwkb_methods[i]=='0'])
    rkwkb_ai_wkb = np.array([x for i,x in enumerate(rkwkb_ai) if rkwkb_methods[i]=='1'])
    rkwkb_ai_rk = np.array([x for i,x in enumerate(rkwkb_ai) if rkwkb_methods[i]=='0'])
    rkwkb_ana_wkb = np.array([x for i,x in enumerate(rkwkb_ana) if rkwkb_methods[i]=='1'])
    rkwkb_ana_rk = np.array([x for i,x in enumerate(rkwkb_ana) if rkwkb_methods[i]=='0'])
    plt.style.use("fyr")
    plt.xlabel('t')
    plt.ylabel('$|\Delta x/ x|$')
    plt.plot(nag_times, abs((nag_ana - nag_ai)/nag_ana), '.', label='NAG solver')
    plt.loglog(rkwkb_times_wkb, abs((rkwkb_ana_wkb - rkwkb_ai_wkb)/rkwkb_ana_wkb), 'x', color = 'green', label='RKWKB, WKB step')
    plt.loglog(rkwkb_times_rk, abs((rkwkb_ana_rk - rkwkb_ai_rk)/rkwkb_ana_rk), 'x', color = 'red', label="RKWKB, RK step")
    plt.legend()
#    plt.savefig("airy-rerrors-corrected.pdf")
    plt.show()
    plt.plot(rkwkb_times, rkwkb_ai, 'x', color='red')
    ts = np.logspace(1.9,2.5,50000)
    plt.plot(ts, airy(-ts)[0])
    plt.show()

def plot_nag_v_rkwkb(nag_inputfile, inputfile, outputfile=None):

    nag_data = np.loadtxt(nag_inputfile, comments='=')
    nag_times = nag_data[:,0] 
    nag_Rs = nag_data[:,1]
    
    data = np.genfromtxt(inputfile, comments='=',
    dtype=(float,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,int))
    times = [line[0] for line in data]
    methods = [line[-1] for line in data]
    xs = [le(line[1])[0] for line in data]

    plt.style.use("fyr")
    plt.plot([times[i] for i,m in enumerate(methods) if m==0], [xs[i] for
    i,m in enumerate(methods) if m==0], 'x', color='red')
    plt.plot([times[i] for i,m in enumerate(methods) if m==1], [xs[i] for
    i,m in enumerate(methods) if m==1], 'x', color='green')
    plt.plot(nag_times,nag_Rs)
#    plt.savefig(outputfile)
    plt.show()
    
def plot_rkwkb(inputfile):
    
    data = np.genfromtxt(inputfile, comments='=',
    dtype=(float,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple,tuple))
    times = [line[0] for line in data]
    xs = [le(line[1])[0] for line in data]
    grid = np.logspace(np.log10(times[0]),np.log10(times[-1]),num=5000)

    plt.style.use("fyr")
    plt.semilogx(times, xs, 'x', color='red', label='rkwkb')
    plt.legend()
#    plt.savefig("ms-example.pdf")
    plt.show()


def main():

    #pass
    plot_airy("outputs/airy-nag.txt", "outputs/airy-rkwkb.txt")
    #nag_inputfile = sys.argv[1]
    #inputfile = sys.argv[2]
    #if(sys.argv[3]):
    #    outputfile = sys.argv[3]
    #plot_nag_v_rkwkb(nag_inputfile, inputfile)
    #plot_rkwkb(inputfile)

if __name__ == "__main__":
    main()
