import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.signal import remez, freqz, lfilter
import Multiple_Source_experiment as Fringe
fs=2.5e6
baseline=14
wavelength= 21e-2
c=3e8
duration=13*60
chunk_size=1
numchunks= int(np.ceil(duration/chunk_size))
tau_instr=0
#3dB Beamwidth 5 deg
wearth=2*np.pi/(24*60*60)
t=np.arange(0,duration,chunk_size)
theta_current=0
theta_current_s2=0.5*np.pi/180
n_samples = int(fs * chunk_size)
# Generate white Gaussian noise for I and Q channels

#First Source
np.random.seed(42)  # For reproducibility
i_s1 = np.random.normal(0,1,n_samples)
q_s1 = np.random.normal(0,1,n_samples) #Explain how std is related to power



np.random.seed(43)  # For reproducibility
i_s2 = np.random.normal(0,1,n_samples)
q_s2= np.random.normal(0,1,n_samples) #Explain how std is related to power

scale=0.1772143622620896
var_s1=10000
var_s2=12000
Phase_Differences=np.zeros(numchunks)
Correlation_Peaks=np.zeros(numchunks)
for i in range(numchunks):
    # generate theta values for this chunk
    theta_radians_s1 = np.linspace(theta_current, theta_current+ chunk_size*wearth, n_samples,endpoint=False)
    tau_geo_s1 = (baseline / c) * np.sin(theta_radians_s1)
    power_factor_s1= np.sinc(scale * np.rad2deg(theta_radians_s1))**2
    theta_radians_s2 = np.linspace(theta_current_s2, theta_current_s2 + chunk_size*wearth, n_samples,endpoint=False)
    
    tau_geo_s2 = (baseline / c) * np.sin(theta_radians_s2)
    power_factor_s2= np.sinc(scale * np.rad2deg(theta_radians_s2))**2

    i_s1_dropped_power= np.sqrt(var_s1*power_factor_s1)*i_s1
    q_s1_dropped_power= np.sqrt(var_s1*power_factor_s1)*q_s1

    i_s2_dropped_power= np.sqrt(var_s2*power_factor_s2)*i_s2
    q_s2_dropped_power= np.sqrt(var_s2*power_factor_s2)*q_s2
    Phase_Differences[i],Correlation_Peaks[i] = Fringe.Phase_Differneces(tau_instr,tau_geo_s1,tau_geo_s2, i_s1_dropped_power,q_s1_dropped_power,i_s2_dropped_power, q_s2_dropped_power, n_samples)
    # advance the current angle for next chunk
    theta_current += chunk_size * wearth
    theta_current_s2+=chunk_size * wearth

  
    
np.savez(f"Two_source_simulated_results_14m_baseline.npz", 
Correlation_Peaks=Correlation_Peaks, 
Phase_Differences=Phase_Differences)
plt.figure(figsize=(10, 4))
plt.plot(t, Phase_Differences)
plt.title(f"Simulated Phase Shift vs Time for 2 Point Sources Offset by 0.5 degrees")
plt.ylabel("Phase Shift [Degrees]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.show()




plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks)
plt.title(f"Simulated Correlation Coefficient vs Time for 2 Point Sources Offset by 0.5 degreess")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)
plt.show()

