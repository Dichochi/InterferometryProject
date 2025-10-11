import numpy as np
import matplotlib.pyplot as plt
x = 10**2*np.random.normal(0, 1, int(2.5e6))
print(np.var(x, ddof=1))

n_samples=int(2.5e6)
i_signal = np.random.normal(0,1,n_samples)
q_signal = np.random.normal(0,1,n_samples) #Explain how std is related to power
var_sun_signal=10000

theta_current=5*np.pi/180

i_signal= np.sqrt(var_sun_signal*(1-0.2*theta_current*180/np.pi))*i_signal
q_signal= np.sqrt(var_sun_signal*(1-0.2*theta_current*180/np.pi))*q_signal


print(np.var(i_signal, ddof=1))  



chunk_size=1
t=np.arange(0,780,chunk_size)
data = np.load("simulated_results_arrays_1s_chucnk_size_Isntr_Delay_0us.npz")
Correlation_Peaks = data["Correlation_Peaks"]
Phase_Differences= data["Phase_Differences"]



plt.figure(figsize=(10, 4))
plt.plot(t, Phase_Differences)
plt.title(f"Simulated Phase Shift vs Time for {chunk_size}s buffer size With Instrumental Delay of {0*1e6} μs")
plt.ylabel("Phase Shift [Radians]")
plt.xlabel("Time [s]")
plt.grid(True)




plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks)
plt.title(f"Measured Correlation Coefficient vs Time for {chunk_size}s buffer size With Instrumental Delay of {0*1e6}s")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)
plt.show()