import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.signal import remez, freqz, lfilter
import TwoAntennaExperiment as Fringe
fs=2.5e6
baseline=12
wavelength= 21e-2
c=3e8
duration=13*60
chunk_size=1
numchunks= int(np.ceil(duration/chunk_size))
tau_instr=0.102e-6
#3dB Beamwidth 5 deg
wearth=2*np.pi/(24*60*60)
t=np.arange(0,duration,chunk_size)
theta_current=0
n_samples = int(fs * chunk_size)
# Generate white Gaussian noise for I and Q channels
np.random.seed(42)  # For reproducibility
i_signal = np.random.normal(0,1,n_samples)
q_signal = np.random.normal(0,1,n_samples) #Explain how std is related to power
var_sun_signal=10000
Phase_Differences=np.zeros(numchunks)
Correlation_Peaks=np.zeros(numchunks)
for i in range(numchunks):
    # generate theta values for this chunk
    theta_radians = np.linspace(theta_current, theta_current + chunk_size*wearth, n_samples,endpoint=False)
    tau_geo = (baseline / c) * np.sin(theta_radians)
    tau_total = tau_geo + tau_instr
    i_signal_reduced_power= np.sqrt(var_sun_signal*(1-0.1*theta_current*180/np.pi))*i_signal
    q_signal_reduced_power= np.sqrt(var_sun_signal*(1-0.1*theta_current*180/np.pi))*q_signal
    Phase_Differences[i],Correlation_Peaks[i] = Fringe.Phase_Differneces(tau_total, i_signal_reduced_power, q_signal_reduced_power, n_samples)
    # advance the current angle for next chunk
    theta_current += chunk_size * wearth

  
    
np.savez(f"simulated_results_arrays_1s_chucnk_size_baseline_12m_Isntr_Delay_0_102us.npz", 
Correlation_Peaks=Correlation_Peaks,
Phase_Differences=Phase_Differences)



plt.figure(figsize=(10, 4))
plt.plot(t, Phase_Differences)
plt.title(f"Simulated Phase Shift vs Time for {chunk_size}s buffer size With Instrumental Delay of {tau_instr*1e6} μs")
plt.ylabel("Phase Shift [Radians]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.show()



plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks)
plt.title(f"Measured Correlation Coefficient vs Time for {chunk_size}s buffer size With Instrumental Delay of {tau_instr*1e6}s")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)
plt.show()
