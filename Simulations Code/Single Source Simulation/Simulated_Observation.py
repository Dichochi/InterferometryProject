import matplotlib.pyplot as plt
import numpy as np
import TwoAntennaExperiment as Fringe


fs=2.5e6
baseline=14
wavelength= 21e-2
c=3e8
duration=13*60
chunk_size=1
numchunks= int(np.ceil(duration/chunk_size))
tau_instr=0.102e-6
t=np.arange(0,duration,chunk_size)
n_samples = int(fs * chunk_size)
np.random.seed(42)  # For reproducibility
i_signal = np.random.normal(0,1,n_samples)
q_signal = np.random.normal(0,1,n_samples) 
var_sun_signal=10000
Phase_Differences=np.zeros(numchunks)
Correlation_Peaks=np.zeros(numchunks)
wearth=2*np.pi/(24*60*60)
theta_current=0
scale=0.1772143622620896 


"""""
def equation(scale):
    return np.sin(np.pi*scale*2.5)/(np.pi*scale*2.5) - 0.707
scale_solution = fsolve(equation, 0.1)[0]
"""

for i in range(numchunks):
    theta_radians = np.linspace(theta_current, theta_current + chunk_size*wearth, n_samples,endpoint=False)
    tau_geo = (baseline / c) * np.sin(theta_radians)
    tau_total = tau_geo + tau_instr
    power_factor= np.sinc(scale * np.rad2deg(theta_radians))**2
    i_signal_reduced_power= np.sqrt(var_sun_signal*power_factor)*i_signal
    q_signal_reduced_power= np.sqrt(var_sun_signal*power_factor)*q_signal
    Phase_Differences[i],Correlation_Peaks[i] = Fringe.Phase_Differneces(tau_total, i_signal_reduced_power, q_signal_reduced_power, n_samples)
    theta_current += chunk_size * wearth

  
    
np.savez(f"simulated_results_arrays_1s_chucnk_size_Isntr_Delay_0_102us_baseline_14m.npz", 
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
