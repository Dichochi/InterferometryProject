import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import seaborn as sns
import os


sample_rate=2.5e6
buffer_duration=1
buffer_size=int(buffer_duration*sample_rate)
complete_duration=13*60
num_buffers=int(complete_duration*sample_rate/buffer_size)
Phase_Differences=np.zeros(int(num_buffers))
Correlation_Peaks=np.zeros(int(num_buffers))
Unnormlised_Peaks=np.zeros(int(num_buffers))
t=np.arange(0,780,buffer_duration)
tau_instr=-1.8*(0.102e-6)
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

        spectrum_ch1= np.fft.fft(samples_Ant_A)  
        spectrum_ch2= np.fft.fft(samples_Ant_B)
        N = len(spectrum_ch1)
        freqs=np.fft.fftfreq(N,1/sample_rate)
        phase_correction=np.exp(2*np.pi*1j*tau_instr*freqs)*np.exp(2*np.pi*1j*tau_instr*1420e6)
        spectrum_ch1*=phase_correction
        Freq_Domain_Cross_Spectrum= spectrum_ch1*np.conj(spectrum_ch2)
        Npad=1

        cross_corr= (np.fft.ifft(Freq_Domain_Cross_Spectrum))
        Ex = np.sum(np.abs(samples_Ant_A)**2)
        Ey = np.sum(np.abs(samples_Ant_B)**2)
        Cxy = cross_corr*Npad / (np.sqrt(Ex) * np.sqrt(Ey))
        Pkidx= np.argmax(np.abs(Cxy))
        Correlation_Peaks[idx]=np.abs(Cxy[Pkidx])
        Phase_Differences[idx]=np.angle(Cxy[Pkidx])
        Unnormlised_Peaks[idx]=np.abs(cross_corr[Pkidx])

        #time_corr = np.correlate(samples_Ant_A, samples_Ant_B,"full") / np.sqrt(Ex * Ey)

        #time_cross_corr=np.abs(signal.correlate(samples_Ant_A,samples_Ant_B))/ (np.sqrt(Ex) * np.sqrt(Ey))
        #print(np.max(time_cross_corr))
        #lags=signal.correlation_lags(2.5e6,2.5e6)/sample_rate

    
        idx+=1
        if flag:
            break







np.savez("results_arrays_corrected_1_7.npz", 
Correlation_Peaks=Correlation_Peaks, 
Phase_Differences=Phase_Differences,
Unnormlised_Peaks=Unnormlised_Peaks)

"""
        N=len(spectrum_ch1)
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
        half_window=74
        time_corr = np.fft.fftshift(Cxy)
        lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate

        main_lobe_Pk = np.argmax(np.abs(time_corr))

        start_idx = main_lobe_Pk - half_window
        end_idx = main_lobe_Pk + half_window + 1

        Cxy_window = time_corr[start_idx:end_idx]
        lags_window = lags[start_idx:end_idx]
        plt.figure(figsize=(10, 4))
        plt.plot(lags_window,Cxy_window)
        plt.grid(True)
        plt.show
       

 Nfft = Npad * N 
        Xspec_shifted = np.fft.fftshift(Freq_Domain_Cross_Spectrum)
        pad_width = (Nfft - N) // 2
        Xspec_padded = np.pad(Xspec_shifted, (pad_width, pad_width), 'constant')
        Xspec_padded = np.fft.ifftshift(Xspec_padded)
        cross_corr = np.fft.ifft(Xspec_padded)


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

        ax.plot(lags, np.abs(time_corr), color='#1f77b4', linewidth=1.2)

        # Reference line at zero
        ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)

        # --- Axis Labels & Title ---
        ax.set_xlabel(r"Time Lag $\tau$ [$\mu$s]", fontsize=12)
        ax.set_ylabel(r"Normalized Cross-Correlation Magnitude, $|C_{xy}|$", fontsize=12)
        ax.set_title(
            f"Normalized Cross-Correlation vs Lag For a single Chunk of {buffer_duration}s",
            fontsize=13, pad=10
        )

        # --- Axes Limits & Grid ---
        ax.set_xlim(-30, 30)
        ax.minorticks_on()

        ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
        ax.tick_params(axis='both', which='major', labelsize=10)

        plt.tight_layout()
        plt.savefig("Cross_Correlation_Chunk.png", dpi=400, bbox_inches='tight')







        time_corr = np.fft.fftshift((Cxy))
        lags = np.arange(-len(time_corr)//2, len(time_corr)//2) / sample_rate
        plt.figure(figsize=(10, 4))
        plt.plot(lags * 1e6,   np.abs(time_corr))
        plt.xlabel("Lag [μs]")
        plt.xlim(-30,30)
        plt.ylabel("Normalized Cross-Correlation")
        plt.title(f"Normalized Cross-Correlation vs Lag For a single Chunk of {buffer_duration}s")
        plt.grid(True)
        plt.show()"""