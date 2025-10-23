import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import seaborn as sns
sample_rate=15e6
num_samples= int(sample_rate*4)
samples_ch1 = np.fromfile('SDRChannel1_part001FirstObservation.iq', dtype=np.complex64, count=num_samples)
samples_ch2 = np.fromfile('SDRChannel2_part001FirstObservation.iq', dtype=np.complex64, count=num_samples)
#t = np.arange(num_samples)/ sample_rate

#FFTfreqs=prp.ComputeFFTFreqs(num_samples,sample_rate)
spectrum_ch1= np.fft.fft(samples_ch1,num_samples)

#prp.PlotFFTMagPhase(FFTmag,FFTphase,FFTfreqs,sample_rate,"Channel 1 Magntitude & Phase Spectrum",False)

spectrum_ch2= np.fft.fft(samples_ch2,num_samples)

#prp.PlotFFTMagPhase(FFTmag,FFTphase,FFTfreqs,sample_rate,"Channel 2 Magntitude & Phase Spectrum",True)


Freq_Domain_Cross_Spectrum= spectrum_ch2*np.conj(spectrum_ch1)
#prp.PlotFFTMagPhase(20*np.log10(np.abs(Freq_Domain_Cross_Spectrum)/num_samples),np.angle(Freq_Domain_Cross_Spectrum),FFTfreqs,sample_rate,"Cross Correlation Freq Domain Result",False)

N = len(Freq_Domain_Cross_Spectrum)
Npad = 1  # zero-padding factor
Nfft = Npad * N

# --- Frequency-domain zero-padding for finer lag resolution ---
Xspec_shifted = np.fft.fftshift(Freq_Domain_Cross_Spectrum)
pad_width = (Nfft - N) // 2
Xspec_padded = np.pad(Xspec_shifted, (pad_width, pad_width), 'constant')
Xspec_padded = np.fft.ifftshift(Xspec_padded)

# --- Cross-correlation computation ---
cross_corr = np.fft.ifft(Xspec_padded)
Ex = np.sum(np.abs(samples_ch1)**2)
Ey = np.sum(np.abs(samples_ch2)**2) 
Cxy = (cross_corr * Npad) / (np.sqrt(Ex) * np.sqrt(Ey))

# --- Time-domain lag axis ---
time_corr = np.fft.fftshift(np.abs(Cxy))
lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'figure.dpi': 100,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.5,
    'lines.linewidth': 1.5,
    'patch.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.6,
    'ytick.minor.width': 0.6,
    'text.usetex': False,
})

sns.set_palette("colorblind")
sns.set_style("whitegrid")

# --- Plot (Thesis-style) ---
fig, ax = plt.subplots(figsize=(10, 4))

ax.plot(lags*1e6, time_corr, color='#1f77b4', linewidth=1.2)

# Reference line at zero
ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)

ax.set_xlabel(r"Time lag $\tau$ [$\mu$s]", fontsize=12)
ax.set_ylabel(r"Normalized Cross-Correlation Magnitude, $|C_{xy}|$", fontsize=12)
ax.set_title(
    r"Cross-Correlation Magnitude between the Two Receive Channels",
    fontsize=13, pad=10
)

# --- Axes Limits & Grid ---
ax.set_xlim(-20, 20)
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.savefig("Cross_Correlation_Calibration.png", dpi=400, bbox_inches='tight')
plt.show()

pkidx = np.argmax(time_corr)
half_window =5
start =pkidx - half_window
end   = pkidx + half_window



x = lags[start:end]
y = time_corr[start:end]

# Fit a 3rd-order polynomial
coeffs = np.polyfit(x, y, 4)
p = np.poly1d(coeffs)

# Fine grid for interpolation
x_fine = np.linspace(x[0], x[-1], 500)
y_fine = p(x_fine)

# Refined peak estimate
fine_idx = np.argmax(y_fine)
refined_lag = x_fine[fine_idx]
print("Refined peak lag:", refined_lag)

# Plot
fig, ax = plt.subplots(figsize=(10, 4))

# --- Convert units ---
lags_us = lags * 1e6
refined_lag_us = refined_lag * 1e6
coarse_lag_us = lags[pkidx] * 1e6

# --- Plot Cross-Correlation Samples ---
ax.plot(
    lags_us[pkidx-20:pkidx+20],
    time_corr[pkidx-20:pkidx+20],
    'o',
    color='#1f77b4',
    markersize=4,
    label=r"Cross-Correlation Samples"
)

# --- Plot 4th-Order Polynomial Fit ---
ax.plot(
    x_fine * 1e6,
    y_fine,
    '-',
    color='#ff7f0e',
    linewidth=1.4,
    label=r"4$^{th}$-Order Polynomial Fit"
)

# --- Add Coarse Lag (Original Estimate) ---
ax.axvline(
    coarse_lag_us,
    color='g',
    linestyle='--',
    linewidth=1,
    alpha=0.8,
    label=fr"Coarse Lag = {coarse_lag_us:.3f} $\mu$s"
)

# --- Add Refined Lag (Interpolated Estimate) ---
ax.axvline(
    refined_lag_us,
    color='r',
    linestyle='--',
    linewidth=1,
    alpha=0.8,
    label=fr"Refined Lag = {refined_lag_us:.3f} $\mu$s"
)

# --- Axis Labels & Title ---
ax.set_xlabel(r"Time Lag $\tau$ [$\mu$s]", fontsize=12)
ax.set_ylabel(r"Normalized Cross-Correlation Magnitude, $|C_{xy}|$", fontsize=12)
ax.set_title(
    r"Cross-Correlation Peak Refinement using 4$^{th}$-Order Polynomial Fit",
    fontsize=13, pad=10
)

# --- Axes Limits & Grid ---
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)

# --- Legend & Layout ---
ax.legend(fontsize=10, loc='best', frameon=True)

plt.tight_layout()
plt.savefig("Cross_Correlation_Refined_Peak.png", dpi=400, bbox_inches='tight')
plt.show()









""" plt.figure(figsize=(10, 4))
plt.plot(lags*1e6, np.abs(np.fft.ifftshift(Cxy)))
plt.title("Normalized Cross Correlation (Zoomed Around Peak)")
plt.xlabel("Lag [µs]")
plt.xlim(-20,20)
plt.ylabel("Normalised Cross Correlation")
plt.grid(True)
plt.show()




N=len(Freq_Domain_Cross_Spectrum)
cross_corr = np.fft.ifft(Freq_Domain_Cross_Spectrum)
Ex = (1/N)*np.sum(np.abs(spectrum_ch1)**2)
Ey = (1/N)*np.sum(np.abs(spectrum_ch2)**2)
Cxy = cross_corr/ (np.sqrt(Ex) * np.sqrt(Ey))


lags = np.arange(-num_samples//2, num_samples//2) / sample_rate
 """



















""" #Frequency Domain Version
cross_corr = np.fft.ifft((Freq_Domain_Cross_Spectrum))
lags = np.arange(-num_samples//2, num_samples//2) / sample_rate

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


peak_idx = np.argmax(np.abs(Freq_Domain_Cross_Spectrum))
    
# Extract window around peak
start = max(0, peak_idx - window_size//2)
end   = min(len(Freq_Domain_Cross_Spectrum), peak_idx + window_size//2)
spectrum_window = Freq_Domain_Cross_Spectrum[start:end]

# Zero-pad window
N = len(spectrum_window)
Nfft = N * Npad
pad_width = (Nfft - N) // 2
spectrum_padded = np.pad(spectrum_window, (pad_width, pad_width), 'constant')

# Compute cross-correlation
cross_corr = np.fft.ifft(spectrum_padded)

# Normalization
Ex = np.sum(np.abs(samples_Ant_A)**2)
Ey = np.sum(np.abs(samples_Ant_B)**2)
Cxy = cross_corr * Npad / (np.sqrt(Ex) * np.sqrt(Ey))

# Shift for plotting
time_corr = np.fft.fftshift(np.abs(Cxy))
lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate

# Plot
plt.figure(figsize=(10,4))
plt.plot(lags*1e6, time_corr)  # convert to microseconds for readability
plt.title("Cross-Correlation (Peak-Centered Zero-Padded)")
plt.xlabel("Lag [µs]")
plt.ylabel("Normalized Magnitude")
plt.grid(True)
plt.show()



 """



