import adi
import numpy as np
duration = 13*60
num_samps_per_buffer = 2**18
sample_rate = 2.5e6
rx_lo = 1420e6
rx_gain0 = 60
rx_gain1 = 60
rx_mode = "manual"
uri = "usb:1.1.5"
num_samps = int(sample_rate * duration)
num_buffers = int(np.ceil(num_samps / num_samps_per_buffer))

sdr = adi.ad9361(uri=uri)
'''Configure Rx properties'''
sdr.rx_enabled_channels = [0, 1]
sdr.sample_rate = int(sample_rate)
sdr.rx_lo = int(rx_lo)
sdr.gain_control_mode = rx_mode
sdr.rx_hardwaregain_chan0 = int(rx_gain0)
sdr.rx_hardwaregain_chan1 = int(rx_gain1)
sdr.rx_buffer_size = int(num_samps_per_buffer)

flag_end=False
numbuffers_in_file=0
file_idx=1
max_buffers_per_file=int(np.ceil(sample_rate*30/num_samps_per_buffer))
num_samps_in_last_buffer = num_samps % num_samps_per_buffer
f1 = open(f'SDRChannel1_part{file_idx:03d}.iq', 'wb')
f2 = open(f'SDRChannel2_part{file_idx:03d}.iq', 'wb')
for r in range(num_buffers):
    data = sdr.rx()
    if r==(num_buffers-1):
        num_samps_in_buffer=num_samps_in_last_buffer   
    else:
        num_samps_in_buffer=num_samps_per_buffer
    buffer_ch1 = np.array(data[0][:num_samps_in_buffer], dtype=np.complex64)
    buffer_ch2 = np.array(data[1][:num_samps_in_buffer], dtype=np.complex64)
    buffer_ch1.tofile(f1)
    buffer_ch2.tofile(f2) 
    numbuffers_in_file+=1
    if numbuffers_in_file>=max_buffers_per_file:
        f1.close()
        f2.close()
        file_idx+=1
        numbuffers_in_file=0
        f1 = open(f'SDRChannel1_part{file_idx:03d}.iq', 'wb')
        f2 = open(f'SDRChannel2_part{file_idx:03d}.iq', 'wb')

f1.close()
f2.close()
