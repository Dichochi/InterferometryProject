import matplotlib.pyplot as plt
import numpy as np
import PreProcessing as prp
from scipy import signal
import os


sample_rate=2.5e6
buffer_duration=0.1
buffer_size=int(buffer_duration*sample_rate)
complete_duration=13*60
num_buffers=int(complete_duration*sample_rate/buffer_size)
Phase_Differences=np.zeros(int(num_buffers))
Correlation_Peaks=np.zeros(int(num_buffers))
t=np.arange(0,780,buffer_duration)
phase_correction= np.exp(-1j*0.005*2*np.pi*t)
idx=0
for  k in range(1, 14):
    if(k<10):
        filename_A = f'SDRChannel1_part00{k}.iq'
        filename_B = f'SDRChannel2_part00{k}.iq'
    else: 
        filename_A = f'SDRChannel1_part0{k}.iq'
        filename_B = f'SDRChannel2_part0{k}.iq'  
    total_samples = os.path.getsize(filename_A)//8  # 8 bytes per complex64
    num_buffers_file = int(np.ceil(total_samples / buffer_size))
    flag=False
    for i in range(1,num_buffers_file):
        offset = i * buffer_size
        if (i == num_buffers_file - 2) and (k != 13):
            num_samps_last=total_samples- (num_buffers_file-1)*buffer_size +buffer_size
            samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=num_samps_last, offset=offset*8)
            samples_Ant_B = np.fromfile(filename_B, dtype=np.complex64, count=num_samps_last, offset=offset*8)
            flag=True
        elif (i == num_buffers_file - 1) and (k == 13):
            num_samps_last=total_samples- (num_buffers_file-1)*buffer_size 
            samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=num_samps_last, offset=offset*8)
            samples_Ant_B = np.fromfile(filename_B, dtype=np.complex64, count=num_samps_last, offset=offset*8)
            flag=True  
        else:
            samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=buffer_size, offset=offset*8)
            samples_Ant_B = np.fromfile(filename_B, dtype=np.complex64, count=buffer_size, offset=offset*8)        
        spectrum_ch1= prp.compute_fft_raw_spectrum(samples_Ant_A,len(samples_Ant_A))
        spectrum_ch2= prp.compute_fft_raw_spectrum(samples_Ant_B,len(samples_Ant_B))
        Freq_Domain_Cross_Spectrum= spectrum_ch1*np.conj(spectrum_ch2)
        #N = len(Freq_Domain_Cross_Spectrum)
        Npad=1
        """Nfft = Npad * N 
        Xspec_shifted = np.fft.fftshift(Freq_Domain_Cross_Spectrum)
        pad_width = (Nfft - N) // 2
        Xspec_padded = np.pad(Xspec_shifted, (pad_width, pad_width), 'constant')
        Xspec_padded = np.fft.ifftshift(Xspec_padded)
        cross_corr = np.fft.ifft(Xspec_padded)"""  
        cross_corr= (np.fft.ifft(Freq_Domain_Cross_Spectrum))
        Ex = np.sum(np.abs(samples_Ant_A)**2)
        Ey = np.sum(np.abs(samples_Ant_B)**2)
        Cxy = cross_corr*Npad / (np.sqrt(Ex) * np.sqrt(Ey))
        Pkidx= np.argmax(np.abs(Cxy))
        Correlation_Peaks[idx]=np.abs(Cxy[Pkidx])
        Phase_Differences[idx]=np.angle(Cxy[Pkidx])

        """time_corr = np.fft.fftshift((Cxy))
        lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate
        plt.figure(figsize=(10, 4))
        plt.plot(lags * 1e6, time_corr)
        plt.xlabel("Lag [μs]")
        plt.ylabel("Normalized Cross-Correlation")
        plt.title(f"Normalized Cross-Correlation vs Lag For a single Chunk of {buffer_duration}s")
        plt.grid(True)
        plt.show()"""
        idx+=1
        if flag:
            break






"""plt.figure(figsize=(10, 4))
np.savez("results_arrays_10ms.npz", 
Correlation_Peaks=Correlation_Peaks, 
Phase_Differences=Phase_Differences)
plt.plot(t, (Phase_Differences1))
plt.title(f"Measured Phase Shift After Correction vs Time for {buffer_duration}s chunk size")
plt.ylabel("Phase Shift [Radians]")
plt.xlabel("Time [s]")
plt.grid(True)


plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks1)
plt.title(f"Measured Correlation Coefficient vs Time for {buffer_duration}s chunk size")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)
Phase_Differences=Phase_Differences[:764]
# Compute amplitude and phase
fft_complex_product= np.fft.fft(Phase_Differences)
freqs = np.fft.fftfreq(len(Phase_Differences), d=1/buffer_duration)  # in Hz
freqs_shifted = np.fft.fftshift(freqs)           # convert to MHz for plotting
fft_abs_product= np.abs(np.fft.fftshift(fft_complex_product))
plt.figure(figsize=(10, 4))
plt.plot(freqs,fft_abs_product)
plt.title("FFT of Phase Difference vs Frequency ")
plt.ylabel("Magnitude")
plt.xlabel("Frequency[Hz]")
plt.grid(True)

plt.show()"""












"""time_corr = np.fft.fftshift(np.abs(Cxy))
lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate
plt.figure(figsize=(10, 4))
plt.plot(lags * 1e6, time_corr)
plt.xlabel("Lag [μs]")
plt.ylabel("Normalized Auto-Correlation")
plt.title(f"Normalized Auto-Correlation vs Lag For a single Chunk of {buffer_duration}s")
plt.grid(True)
plt.show()"""



"""Freq_Domain_Cross_Spectrum=np.fft.fftshift(Freq_Domain_Cross_Spectrum)
freqs= np.fft.fftfreq(len(Freq_Domain_Cross_Spectrum),1/sample_rate)
freqs=np.fft.fftshift(freqs)
plt.figure(figsize=(10, 4))
plt.plot(freqs,10*np.log10(np.abs(Freq_Domain_Cross_Spectrum)/len(Freq_Domain_Cross_Spectrum)))
plt.title(f"Cross Spectrum Magnitude vs Frequency at {(idx+1)*10 +720} s")
plt.ylabel("Magnitude")
plt.xlabel("Frequency[Hz]")
plt.grid(True)
plt.show()""" 