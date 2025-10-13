import numpy as np
from scipy.signal import remez, lfilter, freqz
import matplotlib.pyplot as plt
import hampel
import numpy as np


def generateLPF(Fs=2_500_000):
    # MATLAB specs

    Fpass = 1200000;        
    Fstop = 1250000;        
    Dpass = 0.057501127785;  
    Dstop = 0.01;            
    dens  = 20;              

    # Normalized frequencies for remez (0 to Fs/2)
    bands = [0, Fpass, Fstop, Fs/2]
    desired = [1, 0]  # gain in passband and stopband

    # Weights (inverse of ripple)
    weights = [1/Dpass, 1/Dstop]

    # Estimate order (approx, since scipy has no firpmord)
    # In MATLAB: [N, Fo, Ao, W] = firpmord(...)
    # We'll approximate N by trial. Let's use remez directly.
    numtaps = 71 # you may adjust this based on desired ripple

    # Design filter
    b = remez(numtaps, bands, desired, weight=weights, fs=Fs, grid_density=dens)

    # Frequency response
    #w, h = freqz(b, worN=2048, fs=Fs)

    return b





""" b, w, h = generateLPF()

def generate_IQ_Baseband():
    b,w,h= generateLPF()



plt.plot(w/1e6, 20*np.log10(np.unwrap(np.abs(h))))  # unwrap to avoid phase jumps
plt.title("Filter Phase Response")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Phase [radians]")
plt.grid(True)
plt.show() """