import matplotlib.pyplot as plt
import numpy as np
import PreProcessing as prp
from scipy import signal
sample_rate=2.5e6
num_samples= int(sample_rate*5)
samples_ch1 = np.fromfile('SDRChannel1_part013.iq', dtype=np.complex64, count=num_samples)
samples_ch2 = np.fromfile('SDRChannel2_part013.iq', dtype=np.complex64, count=num_samples)
#t = np.arange(num_samples)/ sample_rate

FFTfreqs=prp.ComputeFFTFreqs(num_samples,sample_rate)
print(FFTfreqs[len(FFTfreqs)//2])
spectrum_ch1= prp.compute_fft_raw_spectrum(samples_ch1,num_samples)

#prp.PlotFFTMagPhase(FFTmag,FFTphase,FFTfreqs,sample_rate,"Channel 1 Magntitude & Phase Spectrum",False)

spectrum_ch2= prp.compute_fft_raw_spectrum(samples_ch2,num_samples)

#prp.PlotFFTMagPhase(FFTmag,FFTphase,FFTfreqs,sample_rate,"Channel 2 Magntitude & Phase Spectrum",True)


Freq_Domain_Cross_Spectrum= spectrum_ch2*np.conj(spectrum_ch1)
#prp.PlotFFTMagPhase(20*np.log10(np.abs(Freq_Domain_Cross_Spectrum)/num_samples),np.angle(Freq_Domain_Cross_Spectrum),FFTfreqs,sample_rate,"Cross Correlation Freq Domain Result",False)


#Frequency Domain Version
cross_corr = np.fft.ifft((Freq_Domain_Cross_Spectrum))
pkidx=np.argmax(cross_corr)
print(np.angle(Freq_Domain_Cross_Spectrum[len(Freq_Domain_Cross_Spectrum)//2]))
print(np.angle(cross_corr[pkidx]))
lags = np.arange(-num_samples//2, num_samples//2) / sample_rate
time_cross_correlation=np.abs(np.fft.fftshift(cross_corr))/np.max(np.abs(cross_corr))

#Scipy version time domain
#cross_corr=np.abs(signal.correlate(samples_ch1,samples_ch2))
#lags=signal.correlation_lags(num_samples,num_samples)/sample_rate


pkidx = np.argmax(time_cross_correlation)
print("Peak lag:", lags[pkidx])

# Define window around the peak (±500 samples → 1000 total)
window = 400
start = max(pkidx - window, 0)
end   = min(pkidx + window, len(time_cross_correlation))

# Plot only the zoomed section
plt.figure(figsize=(10, 4))
plt.plot(lags[start:end]*1e6, time_cross_correlation[start:end])
plt.title("Normalized Cross Correlation (Zoomed Around Peak)")
plt.xlabel("Lag [µs]")
plt.ylabel("Normalised Cross Correlation")
plt.grid(True)


# Window around peak (±15 samples → 30 total)
half_window =5
start =pkidx - half_window
end   = pkidx + half_window



x = lags[start:end]
y = time_cross_correlation[start:end]

# Fit a 3rd-order polynomial
coeffs = np.polyfit(x, y, 5)
p = np.poly1d(coeffs)

# Fine grid for interpolation
x_fine = np.linspace(x[0], x[-1], 500)
y_fine = p(x_fine)

# Refined peak estimate
fine_idx = np.argmax(y_fine)
refined_lag = x_fine[fine_idx]
print("Refined peak lag:", refined_lag)

# Plot
plt.figure(figsize=(10, 4))

# Convert to microseconds
lags_us = lags * 1e6
refined_lag_us = refined_lag * 1e6
coarse_lag_us = lags[pkidx] * 1e6

# Plot around peak
plt.plot(lags_us[pkidx-20:pkidx+20], 
         time_cross_correlation[pkidx-20:pkidx+20], 
         'o', label="Cross Correlation Samples")

plt.plot(x_fine * 1e6, y_fine, '-', label="4th-order fit")

# Add coarse (original) lag
plt.axvline(coarse_lag_us, color="g", linestyle="--", 
            label=f"Coarse lag = {coarse_lag_us:.3f} µs")

# Add refined (interpolated) lag
plt.axvline(refined_lag_us, color="r", linestyle="--", 
            label=f"Refined lag = {refined_lag_us:.3f} µs")

plt.title("Cross-Correlation with 4th Order Interpolation")
plt.xlabel("Lag [µs]")
plt.ylabel("Normalized Cross Correlation")
plt.legend()
plt.grid(True)
plt.show()
