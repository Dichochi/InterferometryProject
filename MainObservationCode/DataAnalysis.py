import matplotlib.pyplot as plt
import numpy as np
import PreProcessing as prp
from scipy import signal
import os

buffer_duration=1
t=np.arange(0,780,buffer_duration)
data = np.load("results_arrays_1s.npz")
Correlation_Peaks = data["Correlation_Peaks"]
Phase_Differences= data["Phase_Differences"]


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

Sxy= Correlation_Peaks*np.exp(1j*Phase_Differences)*np.exp(1j*2*np.pi*t*0.495)

#Phase_Differences=Phase_Differences*(-t*0.03)
plt.figure(figsize=(10, 4))
plt.plot(t, np.angle(Sxy))
plt.title(f"Measured Phase Shift After Correction vs Time for {buffer_duration}s chunk size")
plt.ylabel("Phase Shift [Radians]")
plt.xlabel("Time [s]")
plt.grid(True)


plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks)
plt.title(f"Measured Correlation Coefficient vs Time for {buffer_duration}s chunk size")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)
# Compute amplitude and phase
fft_complex_product= np.fft.fft(Sxy)
freqs = np.fft.fftfreq(len(Sxy), d=1/buffer_duration)  # in Hz
freqs_shifted = np.fft.fftshift(freqs)           # convert to MHz for plotting
fft_abs_product= np.abs(np.fft.fftshift(fft_complex_product))
plt.figure(figsize=(10, 4))
plt.plot(freqs,fft_abs_product)
plt.title("FFT of Phase Difference vs Frequency ")
plt.ylabel("Magnitude")
plt.xlabel("Frequency[Hz]")
plt.grid(True)

plt.show()