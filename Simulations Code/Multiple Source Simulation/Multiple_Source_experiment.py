import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import LowPassFilter as LPF
from scipy.signal import remez, freqz, lfilter


b= LPF.generateLPF()


def Phase_Differneces(tau_instr:float,tau_geo_s1:np.ndarray,tau_geo_s2:np.ndarray,i_s1:np.ndarray, q_s1:np.ndarray, i_s2:np.ndarray, q_s2:np.ndarray, n_samples:int,fs=2.5e6):
    # Create complex baseband signal
    np.random.seed(None)
    isntr_noise_stdv= 20*np.sqrt(5)
    noise_ant1=np.random.normal(0,isntr_noise_stdv,n_samples) + 1j*np.random.normal(0,isntr_noise_stdv,n_samples)
    noise_ant2= np.random.normal(0,isntr_noise_stdv,n_samples) + 1j*np.random.normal(0,isntr_noise_stdv,n_samples)
    s1 = (i_s1 + 1j * q_s1)  
    #s2=(i_s2+ 1j*q_s2)
    #phase_shift = np.exp(1j * 2 * np.pi * np.random.rand())  # random complex phase
    s2 = (i_s2 + 1j*q_s2)

    band_limited_noise_ant1= lfilter(b,1,noise_ant1)
    band_limited_noise_ant2= lfilter(b,1,noise_ant2)
    baseband__s1= lfilter(b,1,s1)
    baseband__s2= lfilter(b,1,s2)

    freqs=np.fft.fftfreq(n_samples,d=1/fs)
    band_limited_noise_ant1_spectrum= np.fft.fft(band_limited_noise_ant1)
    band_limited_noise_ant2_spectrum= np.fft.fft(band_limited_noise_ant2)
    baseband__s1_spectrum= np.fft.fft(baseband__s1)
    baseband__s2_spectrum= np.fft.fft(baseband__s2)


    ant1_spectrum = baseband__s1_spectrum +baseband__s2_spectrum +band_limited_noise_ant1_spectrum
    baseband__s1_spectrum_delayed=baseband__s1_spectrum*np.exp(-1j*2*np.pi*1420e6*tau_geo_s1)*np.exp(-1j*2*np.pi*freqs*tau_geo_s1)
    baseband__s2_spectrum_delayed=baseband__s2_spectrum*np.exp(-1j*2*np.pi*1420e6*tau_geo_s2)*np.exp(-1j*2*np.pi*freqs*tau_geo_s2)
    ant2_spectrum= baseband__s1_spectrum_delayed+baseband__s2_spectrum_delayed +band_limited_noise_ant2_spectrum
    #ant2_spectrum*=np.exp(-1j*2*np.pi*1420e6*tau_instr)*np.exp(-1j*2*np.pi*freqs*tau_instr)

    Freq_Domain_Cross_Spectrum= ant1_spectrum*np.conj(ant2_spectrum)
    N=len(Freq_Domain_Cross_Spectrum)
    cross_corr = np.fft.ifft(Freq_Domain_Cross_Spectrum)
    Ex = (1/N)*np.sum(np.abs(ant1_spectrum)**2)
    Ey = (1/N)*np.sum(np.abs(ant2_spectrum)**2)
    Cxy = cross_corr/ (np.sqrt(Ex) * np.sqrt(Ey))
    Pkidx= np.argmax(np.abs(Cxy))
    """time_corr = np.fft.fftshift(np.abs(Cxy))
    lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / fs
    plt.figure(figsize=(10, 4))
    plt.plot(lags * 1e6, time_corr)
    plt.xlabel("Lag [μs]")
    plt.ylabel("Normalized Cross-Correlation")
    plt.title(f"Normalized Cross-Correlation vs Lag For a single Chunk of 1s")
    plt.grid(True)
    plt.show()"""""

    return np.angle(Cxy[Pkidx]) , np.abs(Cxy[Pkidx])


"""
    N=len(Freq_Domain_Cross_Spectrum)
    Npad=10
    Nfft = Npad * N 
    Xspec_shifted = np.fft.fftshift(Freq_Domain_Cross_Spectrum)
    pad_width = (Nfft - N) // 2
    Xspec_padded = np.pad(Xspec_shifted, (pad_width, pad_width), 'constant')
    Xspec_padded = np.fft.ifftshift(Xspec_padded)
    cross_corr = np.fft.ifft(Xspec_padded)
    Ex = (1/N)*np.sum(np.abs(ant1_spectrum)**2)
    Ey = (1/N)*np.sum(np.abs(ant2_spectrum)**2)
    Cxy = cross_corr*Npad / (np.sqrt(Ex) * np.sqrt(Ey))

time_corr = np.fft.fftshift(np.abs(Cxy))
lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / fs
plt.figure(figsize=(10, 4))
plt.plot(lags * 1e6, time_corr)
plt.xlabel("Lag [μs]")
plt.ylabel("Normalized Cross-Correlation")
plt.title(f"Normalized Cross-Correlation vs Lag For a single Chunk of 1s")
plt.grid(True)
plt.show()"""





""" time_corr = np.fft.fftshift((Cxy))
lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / fs
plt.figure(figsize=(10, 4))
plt.plot(lags * 1e6, time_corr)
plt.xlabel("Lag [μs]")
plt.ylabel("Normalized Cross-Correlation")
plt.title(f"Normalized Cross-Correlation vs Lag For a single Chunk of {duration}s")
plt.grid(True)
plt.show() """









""" # --- Plot magnitude spectrum ---
plt.figure(figsize=(10, 4))
plt.plot(freqs_shifted/1e6, 20*np.log10(np.abs(spectrum_shifted)/len(Baseband_signal)))
plt.title("Frequency Response of Complex Baseband Signal")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Magnitude [dB]")
plt.grid(True)
plt.tight_layout()
plt.show()


phi=2*np.pi*baseline*np.sin(theta_radians)/(wavelength)
phi_wrapped=np.angle(np.exp(1j*phi))
plt.figure(figsize=(10, 4))
plt.plot(theta_degrees, np.rad2deg(phi_wrapped))
plt.title("Angle of arrival vs Expected Phase Shift")
plt.xlabel("Angle of arrival[degrees]")
plt.ylabel("Phase shift[degrees]")
plt.grid(True) """