import adi
import time
import matplotlib.pyplot as plt
import numpy as np
class SDRSystem:
    def __init__(self, 
                 duration=3, 
                 num_samps_per_buffer=2**18, 
                 sample_rate=2.5e6, 
                 rx_lo=1420e6,
                 rx_gain0=60,
                 rx_gain1=60,
                 rx_mode="manual",
                 uri="usb:1.1.5"):
        # Store parameters
        self.duration = duration
        self.num_samps_per_buffer = num_samps_per_buffer
        self.sample_rate = sample_rate
        self.rx_lo = rx_lo
        self.rx_gain0 = rx_gain0
        self.rx_gain1 = rx_gain1
        self.rx_mode = rx_mode
        self.uri = uri

        # Compute dependent parameters
        self.num_samps = int(self.sample_rate * self.duration)
        self.num_buffers = int(np.ceil(self.num_samps / self.num_samps_per_buffer))

        # Initialize SDR object
        self.sdr = None


# Create radio
    def connect(self):
        """Instantiate and configure SDR. Returns True if successful."""
        try:
            self.sdr = adi.ad9361(uri=self.uri)
            self.sdr.rx_enabled_channels = [0, 1]
            self.sdr.sample_rate = int(self.sample_rate)
            self.sdr.rx_lo = int(self.rx_lo)
            self.sdr.gain_control_mode = self.rx_mode
            self.sdr.rx_hardwaregain_chan0 = int(self.rx_gain0)
            self.sdr.rx_hardwaregain_chan1 = int(self.rx_gain1)
            self.sdr.rx_buffer_size = int(self.num_samps_per_buffer)

            # Optional: test read
            _ = self.sdr.rx()
            return True
        except Exception as e:
            print(f"SDR connection failed: {e}")
            return False


    def record_data(self, sample_rate, duration,ch1_gain,ch2_gain):
            """Record IQ samples to file."""
            num_samps = int(sample_rate * duration)
            self.sdr.sample_rate= sample_rate
            self.sdr.rx_hardwaregain_chan0=ch1_gain
            self.sdr.rx_hardwaregain_chan0=ch2_gain
            num_samps_per_buffer= self.num_samps_per_buffer
            num_buffers = int(np.ceil(num_samps / num_samps_per_buffer))
            numbuffers_in_file=0
            file_idx=1
            max_buffers_per_file=int(np.ceil(sample_rate*60/num_samps_per_buffer))
            num_samps_in_last_buffer = num_samps % num_samps_per_buffer
            f1 = open(f'SDRChannel1_part{file_idx:03d}.iq', 'wb')
            f2 = open(f'SDRChannel2_part{file_idx:03d}.iq', 'wb')
            for r in range(num_buffers):
                data = self.sdr.rx()
                # Create buffer for this chunk
                if r==(num_buffers-1):
                    num_samps_in_buffer=num_samps_in_last_buffer   
                else:
                    num_samps_in_buffer=num_samps_per_buffer
                buffer_ch1 = np.array(data[0][:num_samps_in_buffer], dtype=np.complex64)
                buffer_ch2 = np.array(data[1][:num_samps_in_buffer], dtype=np.complex64)
                buffer_ch1.tofile(f1)
                buffer_ch2.tofile(f2) 
                numbuffers_in_file+=1
                if numbuffers_in_file>=max_buffers_per_file and r != num_buffers - 1:
                    f1.close()
                    f2.close()
                    file_idx+=1
                    numbuffers_in_file=0
                    f1 = open(f'SDRChannel1_part{file_idx:03d}.iq', 'wb')
                    f2 = open(f'SDRChannel2_part{file_idx:03d}.iq', 'wb')

            f1.close()
            f2.close()










 



# Rx_0=data[0]
# Rx_1=data[1]
# Rx_total = Rx_0 
# NumSamples = len(Rx_total)
# win = np.hamming(NumSamples)
# y = Rx_total * win
# sp = np.absolute(np.fft.fft(y))
# sp = sp[1:-1]
# sp = np.fft.fftshift(sp)
# s_mag = np.abs(sp) / (np.sum(win)/2)    # Scale FFT by window and /2 since we are using half the FFT spectrum
# s_dbfs = 20*np.log10(s_mag/(2**12))     # Pluto is a 12 bit ADC, so use that to convert to dBFS
# xf = np.fft.fftfreq(NumSamples, 1/samp_rate)
# xf = np.fft.fftshift(xf[1:-1])/1e6
# plt.plot(xf, s_dbfs)
# plt.xlabel("frequency [MHz]")
# plt.ylabel("dBfs")
# plt.draw()
# plt.show()



