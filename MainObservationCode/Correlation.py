import matplotlib.pyplot as plt
import numpy as np
import PreProcessing as prp
from scipy import signal
import os
sample_rate=2.5e6
chunk_size=int(2*sample_rate)
Freq_Domain_Cross_Spectrum=np.zeros(chunk_size,dtype=np.complex64)
window = np.hanning(len(Freq_Domain_Cross_Spectrum))

samples_Ant_A = np.fromfile('SDRChannel1_part001.iq', dtype=np.complex64, count=chunk_size)

theta_degrees=(np.arange(-90,90,0.001))
theta_radians=np.deg2rad(theta_degrees)
phi=2*np.pi*20*np.sin(theta_radians)/(21e-2)
phi_wrapped=np.angle(np.exp(1j*phi))
plt.figure(figsize=(10, 4))
plt.plot(theta_degrees, np.rad2deg(phi_wrapped))
plt.title("Expected Phase Shiftt vs Angle of arrival")
plt.xlabel("Angle of arrival[deg]")
plt.ylabel("Phase shift[degrees]")
plt.grid(True)
plt.show()







idx=0
for  i in range(1, 2):
    if(i<10):
        filename_A = f'SDRChannel1_part00{i}.iq'
        filename_B = f'SDRChannel2_part00{i}.iq'
    else: 
        filename_A = f'SDRChannel1_part0{i}.iq'
        filename_B = f'SDRChannel2_part0{i}.iq'  
    total_samples = os.path.getsize(filename_A)//8  # 8 bytes per complex64
    num_chunks = int(np.ceil(total_samples / chunk_size))
    for i in range(num_chunks):
        offset = i * chunk_size
        # Read and process 1-second chunks
        samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=chunk_size, offset=offset*8)
    
        samples_Ant_B = np.fromfile(filename_B, dtype=np.complex64, count=chunk_size, offset=offset*8)
        spectrum_ch1= prp.compute_fft_raw_spectrum(samples_Ant_A,chunk_size)
        spectrum_ch2= prp.compute_fft_raw_spectrum(samples_Ant_B,chunk_size)

        #Freq_Domain_Cross_Spectrum+= spectrum_ch2*np.conj(spectrum_ch2)
        Freq_Domain_Cross_Spectrum+= spectrum_ch2*np.conj(spectrum_ch2)
        idx+=1


cross_corr = np.fft.ifft((Freq_Domain_Cross_Spectrum))
lags = np.arange(-chunk_size//2, chunk_size//2) / sample_rate
time_cross_correlation=np.abs(np.fft.fftshift(cross_corr))/np.max(np.abs(cross_corr))


        
plt.figure(figsize=(10, 4))
plt.plot(lags*1e6, time_cross_correlation)
plt.title("Normalized Cross Correlation (Zoomed Around Peak)")
plt.xlabel("Lag [µs]")
plt.ylabel("Normalised Cross Correlation")
plt.grid(True)

freqs = np.fft.fftfreq(chunk_size, d=1/sample_rate)  # in Hz
freqs_shifted = np.fft.fftshift(freqs) / 1e6          # convert to MHz for plotting

# Compute amplitude and phase
cross_power_dB = 10*np.log10(np.abs(np.fft.fftshift(Freq_Domain_Cross_Spectrum))/chunk_size)
phase = np.angle(np.fft.fftshift(Freq_Domain_Cross_Spectrum))
phase_unwrapped = np.unwrap(phase)



plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(freqs_shifted,Freq_Domain_Cross_Spectrum )
plt.title("Cross-Power Spectrum Magnitude Plot")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Magnitude[dB]")
plt.grid(True)

# --- Phase ---
plt.subplot(1, 2, 2)
plt.plot(freqs_shifted, phase_unwrapped)
plt.title("Cross-Power Spectrum Phase (∠V(f))")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Phase [radians]")
plt.grid(True)
plt.tight_layout()
plt.show()

 