import adi
import matplotlib.pyplot as plt
import numpy as np

'''Setup'''
samp_rate = 1e6   
NumSamples = 2**10
rx_lo = 1420e6
rx_mode = "manual"  
rx_gain0 = 60
rx_gain1 = 60
fc0 = int(100e3)
phase_cal = 0
num_scans = 4000
d = 0.095         
sdr = adi.ad9361(uri='usb:1.1.5')

'''Configure properties for the Radio'''
sdr.rx_enabled_channels = [0, 1]
sdr.sample_rate = int(samp_rate)
sdr.rx_lo = int(rx_lo)
sdr.gain_control_mode = rx_mode
sdr.rx_hardwaregain_chan0 = int(rx_gain0)
sdr.rx_hardwaregain_chan1 = int(rx_gain1)
sdr.rx_buffer_size = int(NumSamples)
sdr._rxadc.set_kernel_buffers_count(1)  
ts = 1 / float(samp_rate)
xf = np.fft.fftfreq(NumSamples, ts)
xf = np.fft.fftshift(xf)/1e6
signal_start = int(NumSamples*(samp_rate/2+fc0/2)/samp_rate)
signal_end = int(NumSamples*(samp_rate/2+fc0*2)/samp_rate)

def calcTheta(phase):
    arcsin_arg = np.deg2rad(phase)*3E8/(2*np.pi*rx_lo*d)
    arcsin_arg = max(min(1, arcsin_arg), -1)    
    calc_theta = np.rad2deg(np.arcsin(arcsin_arg))
    return calc_theta

def dbfs(raw_data):
    NumSamples = len(raw_data)
    win = np.hamming(NumSamples)
    y = raw_data * win
    s_fft = np.fft.fft(y) / np.sum(win)
    s_shift = np.fft.fftshift(s_fft)
    s_dbfs = 20*np.log10(np.abs(s_shift)/(2**11))     
    return s_dbfs

'''Collect Data'''
for i in range(20):  
    data = sdr.rx()


fig = plt.figure(figsize=(3, 3))
ax = plt.subplot(111, polar=True)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_thetamin(-90)
ax.set_thetamax(90)
ax.set_rlim(bottom=-20, top=0)
ax.set_yticklabels([])


line = ax.vlines(0, 0, -20)
text = ax.text(-2, -14, "0 deg")

plt.ion()  
plt.show()

for i in range(num_scans):
    data = sdr.rx()
    Rx_0 = data[0]
    Rx_1 = data[1]
    
    peak_sum = []
    delay_phases = np.arange(-180, 180, 2)
    
    for phase_delay in delay_phases:
        delayed_Rx_1 = Rx_1 * np.exp(1j * np.deg2rad(phase_delay + phase_cal))
        delayed_sum = dbfs(Rx_0 + delayed_Rx_1)
        peak_sum.append(np.max(delayed_sum[signal_start:signal_end]))
    
    peak_dbfs = np.max(peak_sum)
    peak_delay_index = np.where(peak_sum == peak_dbfs)
    peak_delay = delay_phases[peak_delay_index[0][0]]
    
    if i == 1:
        phase_cal = peak_delay
        with open("PhasedelaysPower.txt", "a") as f:
            f.write(f"{phase_cal}\n")
    steer_angle = int(calcTheta(peak_delay))
    

    line.remove()
    text.remove()
    line = ax.vlines(np.deg2rad(steer_angle), 0, -20)
    text = ax.text(-2, -14, f"{steer_angle} deg")
    
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    plt.pause(0.05)

plt.ioff()  
if i>40: print('\a')  