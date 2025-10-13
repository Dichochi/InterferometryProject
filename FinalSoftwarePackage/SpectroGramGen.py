import numpy as np
import matplotlib.pyplot as plt

def compute_spectrogram(y, fs, N, overlap_factor=0.5):
    # Ensure y is a column vector
    y = np.asarray(y, dtype=np.complex64).flatten()
    ylen = len(y)

    # Window function (Hamming)
    window = np.hamming(N).astype(np.float32)

    # Compute overlap and hop size
    overlap_samples = int(np.floor(overlap_factor * N))
    hop = N - overlap_samples

    # Number of frames
    num_frames = 1 + int(np.floor((ylen - N) / hop))

    # Initialize STFT matrix
    S = np.zeros((N, num_frames), dtype=np.complex64)

    # Loop through frames
    for i in range(num_frames):
        start = i * hop
        end = start + N
        frame = y[start:end]

        # Apply window and remove DC
        frame = frame * window
        frame = frame - np.mean(frame)

        # FFT + fftshift
        Y = np.fft.fftshift(np.fft.fft(frame, n=N))
        S[:, i] = Y

    # Frequency and time axes
    f = np.linspace(-fs/2, fs/2, N, endpoint=False)
    t = (np.arange(num_frames) * hop + N/2) / fs

    return S, f, t

samples_Ant_A = np.fromfile("SDRChannel1_part001.iq", dtype=np.complex64)

fs = 2.5e6  # Sampling frequency (Hz)
duration = 12.0  # seconds
t = np.arange(0, duration, 1/fs)

# Spectrogram parameters
N = 256
overlap_factor = 0.5

# Compute spectrogram
S, f, t_axis = compute_spectrogram(samples_Ant_A, fs, N, overlap_factor)

# Convert to magnitude in dB
S_db = 20 * np.log10(np.abs(S) + 1e-12)

# Plot spectrogram
plt.figure(figsize=(10, 5))
plt.pcolormesh(t_axis, f*1e6, S_db, shading='auto', cmap='viridis')
plt.title("Spectrogram (Manual STFT)")
plt.xlabel("Time [s]")
plt.ylabel("Frequency [MHz]")
plt.colorbar(label="Magnitude [dB]")
plt.tight_layout()
plt.show()