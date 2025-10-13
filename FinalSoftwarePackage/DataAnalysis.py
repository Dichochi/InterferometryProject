import numpy as np
import os
from scipy.signal import lfilter
import LowPassFilter as LPF  # Your custom filter module

def simulate_interferometer(duration=60, sample_rate=2.5e6, baseline=12, isntr_noise_stdv=20):
    c = 3e8
    wearth = 2 * np.pi / (24 * 60 * 60)
    buffer_duration = 1
    buffer_size = int(buffer_duration * sample_rate)

    num_buffers = int(duration * sample_rate / buffer_size)
    Phase_Differences = np.zeros(num_buffers)
    Correlation_Peaks = np.zeros(num_buffers)

    b = LPF.generateLPF()  # Low-pass filter coefficients
    theta_current = 0
    idx = 0


    num_files = int(np.ceil(duration / 60))

    for k in range(1, num_files + 1):
        # Construct filename
        if k < 10:
            filename_A = f"SDRChannel1_part00{k}.iq"
        else:
            filename_A = f"SDRChannel1_part0{k}.iq"

        total_samples = os.path.getsize(filename_A) // 8  # 8 bytes per complex64
        num_buffers_file = int(np.ceil(total_samples / buffer_size))
        flag = False

        for i in range(1, num_buffers_file):
            offset = i * buffer_size

            if (i == num_buffers_file - 2 and k != num_files):
                num_samps_last = total_samples - (num_buffers_file - 1) * buffer_size + buffer_size
                samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=num_samps_last, offset=offset * 8)
                flag = True
            elif (i == num_buffers_file - 1 and k == num_files):
                num_samps_last = total_samples - (num_buffers_file - 1) * buffer_size
                samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=num_samps_last, offset=offset * 8)
                flag = True
            else:
                samples_Ant_A = np.fromfile(filename_A, dtype=np.complex64, count=buffer_size, offset=offset * 8)

            # --- Spectrum and synthetic second channel ---
            spectrum_ch1 = np.fft.fft(samples_Ant_A)
            N = len(spectrum_ch1)
            freqs = np.fft.fftfreq(N, d=1 / sample_rate)

            noise_ant2 = np.random.normal(0, isntr_noise_stdv, N) + 1j * np.random.normal(0, isntr_noise_stdv, N)
            band_limited_noise_ch2 = lfilter(b, 1, noise_ant2)
            band_limited_noise_ch2_spectrum = np.fft.fft(band_limited_noise_ch2)

            # --- Apply geometric delay ---
            theta_radians = np.linspace(theta_current, theta_current + buffer_duration * wearth, N, endpoint=False)
            tau_geo = (baseline / c) * np.sin(theta_radians)

            spectrum_ch2 = (
                spectrum_ch1
                * np.exp(-1j * 2 * np.pi * 1420e6 * tau_geo)
                * np.exp(-1j * 2 * np.pi * freqs * tau_geo)
                + band_limited_noise_ch2_spectrum
            )

            # --- Cross correlation ---
            Freq_Domain_Cross_Spectrum = spectrum_ch1 * np.conj(spectrum_ch2)
            cross_corr = np.fft.ifft(Freq_Domain_Cross_Spectrum)
            Ex = (1 / N) * np.sum(np.abs(spectrum_ch1) ** 2)
            Ey = (1 / N) * np.sum(np.abs(spectrum_ch2) ** 2)
            Cxy = cross_corr / (np.sqrt(Ex) * np.sqrt(Ey))

            Pkidx = np.argmax(np.abs(Cxy))
            Correlation_Peaks[idx] = np.abs(Cxy[Pkidx])
            Phase_Differences[idx] = np.angle(Cxy[Pkidx])

            theta_current += buffer_duration * wearth
            idx += 1
            if flag:
                break

    return Phase_Differences, Correlation_Peaks


simulate_interferometer(60)