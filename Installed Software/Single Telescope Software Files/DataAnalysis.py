import numpy as np
import os
from scipy.signal import lfilter
import LowPassFilter as LPF  # Your custom filter module
import matplotlib.pyplot as plt
def simulate_interferometer(duration=60, sample_rate=2.5e6,buffer_duration=1, baseline=14, isntr_noise_stdv=40):
    c = 3e8
    wearth = 2 * np.pi / (24 * 60 * 60)
    buffer_duration = buffer_duration
    buffer_size = int(buffer_duration * sample_rate)

    num_buffers = int(duration * sample_rate / buffer_size)
    Phase_Differences = np.zeros(num_buffers+1)
    Correlation_Peaks = np.zeros(num_buffers+1)
    Cxy_window = np.zeros(201)
    lags_window = np.zeros(201)

    b = LPF.generateLPF()  # Low-pass filter coefficients
    theta_current = 0
    idx = 0
    t=np.arange(0,duration+buffer_duration,buffer_duration)

    num_files = int(np.ceil(duration / 60))
    save_dir = r"C:\Users\Utente\Documents\Interferometer Data"

    for k in range(1, num_files + 1):
        # Construct filename
        if k < 10:
            filename_A = os.path.join(save_dir, f"SDRChannel1_part00{k}.iq")
        else:
            filename_A = os.path.join(save_dir, f"SDRChannel1_part0{k}.iq")

        total_samples = os.path.getsize(filename_A) // 8  # 8 bytes per complex64
        num_buffers_file = int(np.ceil(total_samples / buffer_size))
        flag = False
    

        for i in range(num_buffers_file):
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
            Pkidx = np.argmax(np.abs(cross_corr))
            Correlation_Peaks[idx] = np.abs(cross_corr[Pkidx])
            Phase_Differences[idx] = np.angle(cross_corr[Pkidx])



            if(i==0) and (k==1):
                Npad=8
                Nfft = Npad * N 
                Ex = (1 / N) * np.sum(np.abs(spectrum_ch1) ** 2)
                Ey = (1 / N) * np.sum(np.abs(spectrum_ch2) ** 2)
                Xspec_shifted = np.fft.fftshift(Freq_Domain_Cross_Spectrum)
                pad_width = (Nfft - N) // 2
                Xspec_padded = np.pad(Xspec_shifted, (pad_width, pad_width), 'constant')
                Xspec_padded = np.fft.ifftshift(Xspec_padded)
                cross_corr = np.fft.ifft(Xspec_padded) 
                Cxy = cross_corr*Npad / (np.sqrt(Ex) * np.sqrt(Ey))
                
                # Shift once and use consistently
                half_window=74
                time_corr = np.fft.fftshift(Cxy)
                lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate

                main_lobe_Pk = np.argmax(np.abs(time_corr))

                start_idx = main_lobe_Pk - half_window
                end_idx = main_lobe_Pk + half_window + 1

                Cxy_window = time_corr[start_idx:end_idx]
                lags_window = lags[start_idx:end_idx]


            theta_current += buffer_duration * wearth
            idx += 1
            if flag:
                break

    return Cxy_window, lags_window, Phase_Differences[:-2], Correlation_Peaks[:-2],t[:-2]


