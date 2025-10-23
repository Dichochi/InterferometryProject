import adi
import numpy as np
import matplotlib.pyplot as plt


duration = 60*60  
num_samps_per_buffer = 2**18
sample_rate = 2.5e6
rx_lo = 1420e6
rx_gain0 = 30
rx_gain1 = 30
rx_mode = "manual"
uri = "usb:1.1.5"


num_samps = int(sample_rate * duration)
num_buffers = int(np.ceil(num_samps / num_samps_per_buffer))


sdr = adi.ad9361(uri=uri)
sdr.rx_enabled_channels = [0, 1]
sdr.sample_rate = int(sample_rate)
sdr.rx_lo = int(rx_lo)
sdr.gain_control_mode = rx_mode
sdr.rx_hardwaregain_chan0 = int(rx_gain0)
sdr.rx_hardwaregain_chan1 = int(rx_gain1)
sdr.rx_buffer_size = int(num_samps_per_buffer)


phase_array = np.zeros(num_buffers, dtype=np.float64)
PeaksArray = np.zeros(num_buffers, dtype=np.float64)


window = np.hanning(num_samps_per_buffer)

for r in range(num_buffers):
    data = sdr.rx() 
    ch1 = data[0] * window
    ch2 = data[1] * window


    X1 = np.fft.fft(ch1)
    X2 = np.fft.fft(ch2)


    S12 = X1 * np.conj(X2)

    k = np.argmax(np.abs(S12))
    phase_array[r] = np.angle(S12[k])

  
    Cxy = np.fft.ifft(S12)
    Ex = np.sum(np.abs(ch1)**2)
    Ey = np.sum(np.abs(ch2)**2)
    Cxy_norm = np.abs(Cxy) / np.sqrt(Ex * Ey)

    PeaksArray[r] = np.max(Cxy_norm)


plt.figure(figsize=(10,4))
plt.plot(np.unwrap(phase_array))
plt.xlabel('Chunk index')
plt.ylabel('Phase (rad)')
plt.title('Phase at Peak of Cross-Spectrum')
plt.grid(True)
plt.show()


np.savez("Phase_Peaks_60minGain30dB.npz", phase_array=phase_array, PeaksArray=PeaksArray)



