
import numpy as np
from scipy.signal import remez, lfilter, freqz
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import remez, freqz




def ComputeFFTFreqs(n_samples,sample_rate: int)-> np.ndarray:
    return np.fft.fftshift(np.fft.fftfreq(n_samples, 1 / sample_rate)) / 1e6



def compute_fft_raw_spectrum(x_array: np.ndarray,N:int):
 
    spectrum2 = (np.fft.fft(x_array,N))

    return spectrum2




   
def compute_fft(x_array: np.ndarray):
    n_samples = len(x_array)
    window = np.hanning(n_samples)
    x_window = x_array * window
    # Compute FFT and shift zero frequency to center
    spectrum = np.fft.fftshift(np.fft.fft(x_window))
    spectrum2 = (np.fft.fft(x_array))
    # Magnitude in dB normalized to reference
    spectrum_db = 20 * np.log10(np.abs(spectrum) / n_samples)

    # Phase in radians
    spectrum_phase = np.unwrap(np.angle(spectrum))

    return  spectrum_db, spectrum_phase



def PlotFFTMag(spectrum_db: np.ndarray, freqs: np.ndarray,
               sample_rate: int, name: str, show: bool):
    # ---- Magnitude Plot ----
    plt.figure(figsize=(10, 5))
    plt.plot(freqs, spectrum_db, lw=1.2)
    print(np.max(spectrum_db))
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Magnitude [dB]")
    plt.title(name)
    plt.ylim(-250, 50)
    #plt.xlim(0, sample_rate * 1e-6 / 2)
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    if show:
        plt.show()

def PlotFFTPhase(spectrum_phase: np.ndarray, freqs: np.ndarray,
               sample_rate: int, name: str, show: bool):
    # ---- Magnitude Plot ----
    plt.figure(figsize=(10, 5))
    plt.plot(freqs, spectrum_phase, lw=1.2)
    plt.ylabel("Unwarpped Phase [Radians]")
    plt.xlabel("Frequency[MHz]")
    plt.title(name)
    #plt.xlim(0, sample_rate * 1e-6 / 2)
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.tight_layout()
    if show:
        plt.show()


def PlotFFTMagPhase(spectrum_db: np.ndarray, spectrum_phase: np.ndarray, freqs: np.ndarray,
                    sample_rate: int, name: str, show: bool):
    """
    Plot magnitude and phase as stacked subplots
    Annotate phase at ~100 kHz and show local slope around it
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # ---- Magnitude Plot ----
    ax1.plot(freqs, spectrum_db, lw=1.2, color='tab:blue')
    ax1.set_ylabel("Magnitude [dB]")
    ax1.set_title(name)
    #ax1.set_ylim(-250, 50)
    ax1.grid(True, which="both", ls="--", alpha=0.6)
    
    # ---- Phase Plot ----
    #phase_unwrapped = np.unwrap(spectrum_phase)
    ax2.plot(freqs, spectrum_phase, lw=1.2, color='tab:red')
    ax2.set_xlabel("Frequency [MHz]")
    ax2.set_ylabel("Phase [Radians]")
    ax2.grid(True, which="both", ls="--", alpha=0.6)

    plt.tight_layout()
    
    if show:
        plt.show()

    
#b, w, h = generateLPF()
# plt.plot(w/1e6, np.unwrap(np.angle(h)))  # unwrap to avoid phase jumps
# plt.title("Filter Phase Response")
# plt.xlabel("Frequency [MHz]")
# plt.ylabel("Phase [radians]")
# plt.grid(True)
# plt.show()
#ant1_signal=(sdata.generate_Hydrogen_signal())
#PolyPhaseFilterBank.PlotFrequencySpectrum(ant1_signal)
#FFTMagnitude, FFTPhase = compute_fftMag(ant1_signal)
#PlotFFTMag(FFTMagnitude,FFTPhase ,ComputeFFTFreqs(sdata.n_samples,sdata.sample_rate))
# # Using defaults
# plt.figure(figsize=(10,5))
# plt.plot(w/1e6, 20*np.log10(np.abs(h)))
# plt.title('Equiripple FIR Bandpass Frequency Response')
# plt.xlabel('Frequency [MHz]')
# plt.ylabel('Amplitude [dB]')
# plt.grid(True)
# plt.show()

