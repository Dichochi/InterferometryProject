import matplotlib.pyplot as plt
import numpy as np
import PreProcessing as prp
from scipy import signal
import os


sample_rate=2.5e6
buffer_duration=10
buffer_size=int(buffer_duration*sample_rate)
complete_duration=13*60
num_buffers=int(complete_duration*sample_rate/buffer_size)
Phase_Differences=np.zeros(int(num_buffers))
Cross_Power_products=np.zeros(int(num_buffers),dtype=np.complex64)
t=np.arange(0,780,buffer_duration)
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
    for i in range(num_buffers_file):
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
        """Freq_Domain_Cross_Spectrum=np.fft.fftshift(Freq_Domain_Cross_Spectrum)
        freqs= np.fft.fftfreq(len(Freq_Domain_Cross_Spectrum),1/sample_rate)
        freqs=np.fft.fftshift(freqs)
        plt.figure(figsize=(10, 4))
        plt.plot(freqs,10*np.log10(np.abs(Freq_Domain_Cross_Spectrum)/len(Freq_Domain_Cross_Spectrum)))
        plt.title(f"Cross Spectrum Magnitude vs Frequency at {(idx+1)*10+720} s")
        plt.ylabel("Magnitude")
        plt.xlabel("Frequency[Hz]")
        plt.grid(True)
        plt.show()""" 

        idx_fc=len(Freq_Domain_Cross_Spectrum)//2
        Phase_Differences[idx]=np.angle(Freq_Domain_Cross_Spectrum[idx_fc])
        Cross_Power_products[idx]= Freq_Domain_Cross_Spectrum[idx_fc]

        idx+=1
        if flag:
            break


wavelength=21e-2
baseline=20

theta_degrees=(np.arange(-90,90,0.001))
theta_radians=np.deg2rad(theta_degrees)
phi=2*np.pi*baseline*np.sin(theta_radians)/(wavelength)
phi_wrapped=np.angle(np.exp(1j*phi))
plt.figure(figsize=(10, 4))
plt.plot(theta_degrees, np.rad2deg(phi_wrapped))
plt.title("Angle of arrival vs Expected Phase Shift")
plt.xlabel("Angle of arrival[degrees]")
plt.ylabel("Phase shift[degrees]")
plt.grid(True)



thetafactor= (Phase_Differences*wavelength)/(2*np.pi*baseline)
theta_measured_degrees=np.rad2deg(np.arcsin((thetafactor)))
theta_measured_radians= np.arcsin((thetafactor))
plt.figure(figsize=(10, 4))
plt.plot(theta_measured_degrees, np.rad2deg(Phase_Differences))
plt.title(" Measured Phase Shift vs Angle of arrival")
plt.ylabel("Phase Shift [degrees]")
plt.xlabel("Angle of Arrival [degrees]")
plt.grid(True)


plt.figure(figsize=(10, 4))
plt.plot(t, np.rad2deg(Phase_Differences))
plt.title("Measured Phase Shift vs Time ")
plt.ylabel("Phase Shift [degrees]")
plt.xlabel("Time [s]")
plt.grid(True)

measured_delay= baseline*np.sin(theta_measured_radians)/3E8
plt.figure(figsize=(10, 4))
plt.plot(t, measured_delay)
plt.title("Measured dealay vs Time ")
plt.ylabel("Measured Delay [s]")
plt.xlabel("Time [s]")
plt.grid(True)








# Compute amplitude and phase
fft_complex_product= np.fft.fft(Cross_Power_products)

plt.figure(figsize=(10, 4))
plt.plot(t,np.abs(fft_complex_product))
plt.title("Centre Frequency Magnitude vs Time ")
plt.ylabel("Magnitude")
plt.xlabel("Time[s]")
plt.grid(True)

freqs = np.fft.fftfreq(num_buffers, d=1)  # in Hz
freqs_shifted = np.fft.fftshift(freqs)           # convert to MHz for plotting
fft_abs_product= np.abs(np.fft.fftshift(fft_complex_product))
plt.figure(figsize=(10, 4))
plt.plot(freqs,fft_abs_product)
plt.title("Cross Spectrum Magnitude vs Frequency ")
plt.ylabel("Magnitude")
plt.xlabel("Frequency[Hz]")
plt.grid(True)


Z=np.exp(1j*Phase_Differences) 
Z_fft = np.fft.fft(Z - np.mean(Z))  # remove DC component
Z_fft_shifted = np.fft.fftshift(Z_fft)

plt.figure(figsize=(10, 4))
plt.plot(freqs,np.abs(Z_fft_shifted))
plt.title("Cross Spectrum Magnitude vs Frequency ")
plt.ylabel("Magnitude")
plt.xlabel("Frequency[Hz]")
plt.grid(True)


plt.show()





