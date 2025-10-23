import numpy as np
from scipy.signal import remez
import numpy as np


def generateLPF(Fs=2_500_000):
    Fpass = 1200000;        
    Fstop = 1250000;        
    Dpass = 0.057501127785;  
    Dstop = 0.01;            
    dens  = 20;              
    bands = [0, Fpass, Fstop, Fs/2]
    desired = [1, 0] 
    weights = [1/Dpass, 1/Dstop]
    numtaps = 71 
    b = remez(numtaps, bands, desired, weight=weights, fs=Fs, grid_density=dens)
    return b

