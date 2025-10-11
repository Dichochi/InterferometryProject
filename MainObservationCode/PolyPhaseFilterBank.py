import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.fft import fft, fftfreq, fftshift

M=8 #Number of Taps
P=1024 #Number of Polyphase Branches= Number of fft channels

def zero_pad(x):
    block_size = M * P
    L = len(x)
    remainder = L % block_size
    if remainder != 0:
        pad_len = block_size - remainder
        x = np.concatenate([x, np.zeros(pad_len)])
    return x

def generate_win_coeffs(M, P, window_fn="hamming"):
    win_coeffs = scipy.signal.get_window(window_fn, M*P)
    sinc       = scipy.signal.firwin(M * P, cutoff=1.0/P, window="rectangular")
    win_coeffs *= sinc
    return win_coeffs



def pfb_fir_frontend(x, win_coeffs, M, P):
    x=zero_pad(x)
    W = int(x.shape[0] / M / P)
    x_p = x.reshape((W*M, P)).T
    h_p = win_coeffs.reshape((M, P)).T
    x_summed = np.zeros((P, M * W - M))
    for t in range(0, M*W-M):
        x_weighted = x_p[:, t:t+M] * h_p
        x_summed[:, t] = x_weighted.sum(axis=1)
    return x_summed.T

def pfb_fft(x_p, P, axis=1):
    return (np.fft.fft(x_p, P, axis=axis))

def pfb_filterbank(x, win_coeffs, M, P):
    x_fir = pfb_fir_frontend(x, win_coeffs, M, P)
    x_pfb = pfb_fft(x_fir, P)
    return x_pfb


def compute_Channelised(x:np.ndarray,):
    win_coeffs= generate_win_coeffs(M,P)
    return pfb_filterbank(x,win_coeffs,M,P)
    

def correctdelay(X_Channelised:np.ndarray,tau_inst:float,fs:int):
    freqs = np.fft.fftfreq(P, d=1/fs)  
    phase_ramp = np.exp(1j * 2 * np.pi * freqs * tau_inst)
    return X_Channelised * phase_ramp[None, :]


def PFBSpectrometer(spectrum:np.ndarray,fs:int,n_int:int= 10):  
    x_psd = np.abs(spectrum)**2
    x_psd = x_psd[:int(x_psd.shape[0]//n_int)*n_int] # Trim array so we can do time integration
    x_psd = x_psd.reshape(x_psd.shape[0]//n_int, n_int, x_psd.shape[1])# Integrate over time by reshaping and averaging (efficient)
    pfb_integrated = x_psd.mean(axis=1)  # Average over n_int blocks
    avg_spectrum = np.mean(pfb_integrated, axis=0)  # Average across all integrated blocks
    # Convert to dB (using integrated power directly)
    spectrum_db = 10 * np.log10(avg_spectrum / len(avg_spectrum))
    # Create frequency axis
    freqs = np.fft.fftshift(np.fft.fftfreq(P, 1/fs)) / 1e6
    plt.figure(figsize=(12, 6))
    plt.plot(freqs, np.fft.fftshift(spectrum_db))  # spectrum_db is now 1D
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (dB)')
    plt.title('PFB Power Spectrum (Time Integrated)')
    plt.grid(True)
    plt.show()



"""




 

def pfb_delayed(x_pfb:np.ndarray,tau:float):
    freq_bins = np.fft.fftfreq(P, 1/fs)
    delay_arr= np.exp(-1j * 2 * np.pi * freq_bins * tau)
    X_delayed = x_pfb * delay_arr[None,:]
    return X_delayed 


def two_antenna_correlation(X1, X2, normalize=True, axis=0):
    # Cross-correlation
    prod = X1 * np.conj(X2)   # (num_blocks, P)
    V12 = np.mean(prod, axis=axis)
    return V12

def plot_interferometer_output(V12):

    Plot amplitude and phase of cross-correlation across all channels.
    
    # Frequency axis
    title="Two-Antenna Correlator Output"
    freqs = np.fft.fftfreq(P, 1/fs)

    # Amplitude & Phase
    amp = np.abs(V12)
    phase = np.unwrap(np.angle(V12))

    plt.figure(figsize=(12,6))

    # Amplitude
    plt.subplot(2,1,1)
    plt.plot(amp)
    plt.title("Cross-correlation Amplitude per Channel")

    plt.subplot(2,1,2)
    plt.plot(phase)
    plt.title("Cross-correlation Phase per Channel")
    plt.tight_layout()
    plt.show() """

""" x_signal = sdata.generate_Hydrogen_signal()  # shape: (num_samples,)

# 2. Prepare PFB
win_coeffs = generate_win_coeffs(M, P, window_fn="hamming")
x = zero_pad(x_signal)
X = pfb_filterbank(x, win_coeffs, M, P)  # shape: (num_blocks, P)

X_delayed=pfb_delayed(X,1e-6) 


V12=two_antenna_correlation(X,X_delayed)
plot_interferometer_output(V12) """


""" 
def TestDelay(tau_frac=1e-6):
    
    Test fractional delay using the PFB.
    
    tau_frac: fractional delay in seconds
    Fs: sampling rate in Hz
    
    # 1. Generate test signal
    x_signal = sdata.generate_Hydrogen_signal()  # shape: (num_samples,)

    # 2. Prepare PFB
    win_coeffs = generate_win_coeffs(M, P, window_fn="hamming")
    x = zero_pad(x_signal)
    X = pfb_filterbank(x, win_coeffs, M, P)  # shape: (num_blocks, P)

    X_delayed=pfb_delayed(X,1e-6) 

    # 3. Choose a channel to test
    channel_idx = 1800
    tstch = X[:, channel_idx]  # complex samples in frequency domain

    # 4. Compute channel frequency relative to baseband
    delta_f = fs / P
    f_k = channel_idx * delta_f  # frequency of this channel

    # 5. Apply fractional delay (phase ramp)
    tstch_delayed = X_delayed[:, channel_idx] 

    # 6. Compute expected phase shift
    expected_phase = (-2 * np.pi * f_k * tau_frac ) # radians

    # 7. Measured phase shift
    # Take mean phase difference between original and delayed channel
    measured_phase = np.angle(np.mean(tstch_delayed * np.conj(tstch)))

    print(f"Channel {channel_idx} frequency: {f_k:.2f} Hz")
    print(f"Expected phase shift: {expected_phase:.4f} rad")
    print(f"Measured phase shift: {measured_phase:.4f} rad")
 """
# Example usage:
#TestDelay(tau_frac=1e-6)


