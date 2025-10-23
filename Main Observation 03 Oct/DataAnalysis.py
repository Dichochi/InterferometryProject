import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
data = np.load("results_arrays_1s.npz")
Correlation_Peaks = data["Correlation_Peaks"]
Phase_Differences= data["Phase_Differences"]
Unnormlised_Peaks=data["Unnormlised_Peaks"]
t=np.arange(0,13*60,1)


unwrapped_phase = np.unwrap(Phase_Differences)
print(Phase_Differences[-2])
num_wraps = (unwrapped_phase[-1] - unwrapped_phase[0]) / (2 * np.pi)

print(f"Total number of 2π wraps: {num_wraps:.4f}")

print(np.std(Correlation_Peaks[:779]))
print(np.mean(Correlation_Peaks[:779]))
print(np.max(Correlation_Peaks[:779]))
print(np.min(Correlation_Peaks[:779]))


# --- Plot Styling ---
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

wearth=2*np.pi/(24*60*60)
Sxy=Correlation_Peaks*np.exp(1j*(Phase_Differences))
N=20*len(Sxy)
theta=np.rad2deg(np.linspace(-wearth*20*60,wearth*20*60+wearth*15.5,N))

fig, ax = plt.subplots(figsize=(10, 4))

# --- Plot Phase Differences in Degrees ---
ax.plot(theta, np.abs(np.fft.fftshift(np.fft.ifft(Sxy,N))), color='#1f77b4', linewidth=1.2)


ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)


ax.set_xlabel(r"Angular Span of Observation $\theta$°", fontsize=12)
ax.set_ylabel(r"Magnitude", fontsize=12)
ax.set_title(
    r"1D Reconstructed Image For 13 minute solar observation",
    fontsize=13, pad=10
)


ax.set_xlim(-0.2, 0.2)
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.savefig("1D_Image_main_Observ.png", dpi=400, bbox_inches='tight')
plt.show()






fig, ax = plt.subplots(figsize=(10, 4))

ax.plot(t, Correlation_Peaks, color='#1f77b4', linewidth=1.2)


ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)


ax.set_xlabel(r"Time [s]", fontsize=12)
ax.set_ylabel(r"Correlation Coefficient ($r$)", fontsize=12)
ax.set_title(
    r"Observed Correlation Coefficient during 13-Minute Solar Observation",
    fontsize=13, pad=10
)

# --- Axes Limits & Grid ---
ax.set_xlim(t[0], t[-2])
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.savefig("Cross_Correlation_Coefficients.png", dpi=400, bbox_inches='tight')





fig, ax = plt.subplots(figsize=(10, 4))

ax.plot(t, Unnormlised_Peaks, color='#1f77b4', linewidth=1.2)


ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)


max_val=Unnormlised_Peaks[0]
half_power = max_val / 2

half_idx = np.argmin(np.abs(Unnormlised_Peaks[:779] - half_power))
half_time = t[half_idx]


ax.axvline(half_time, color='r', linestyle='--', linewidth=1, label=f"Half-Power Point = {half_time:.2f} s")
ax.text(half_time + 2, half_power, 'Half-Power Point', color='r', fontsize=10, rotation=90, va='bottom')


ax.set_xlabel(r"Time [s]", fontsize=12)
ax.set_ylabel(r"Cross-Correlation Magnitude [arb. units]", fontsize=12)
ax.set_title(
    r"Observed Unnormalized Cross-Correlation Magnitude during 13-Minute Solar Observation",
    fontsize=13, pad=10
)


ax.set_xlim(t[0], t[778])
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.legend()

plt.tight_layout()
plt.savefig("UnnormalisedCross_Correlation_Coefficients.png", dpi=400, bbox_inches='tight')
plt.show()





fig, ax = plt.subplots(figsize=(10, 4))

# --- Plot Phase Differences in Degrees ---
ax.plot(t, np.rad2deg(Phase_Differences), color='#1f77b4', linewidth=1.2)

# --- Reference Line at Zero ---
ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)

# --- Axis Labels & Title ---
ax.set_xlabel(r"Time [s]", fontsize=12)
ax.set_ylabel(r"Fringe Phase Difference $\phi$°", fontsize=12)
ax.set_title(
    r"Observed Fringe Phase Variation during 13-Minute Solar Observation",
    fontsize=13, pad=10
)

# --- Axes Limits & Grid ---
ax.set_xlim(t[0], t[-1])
ax.set_ylim(-182,182)
ax.minorticks_on()

ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.6)
ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='gray', alpha=0.4)
ax.tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.savefig("Fringe_Pattern.png", dpi=400, bbox_inches='tight')
















plt.figure(figsize=(10, 4))
plt.plot(t, Correlation_Peaks)
plt.title(f"Measured Correlation Coefficient vs Time ")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)

plt.figure(figsize=(10, 4))
plt.plot(t, Phase_Differences)
plt.title(f"Measured Phase Differences vs Time ")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)

plt.figure(figsize=(10, 4))
plt.plot(t, Unnormlised_Peaks)
plt.title(f"Measured unnormalised Correlation Coefficient vs Time ")
plt.ylabel("Correlation Coefficient")
plt.xlabel("Time [s]")
plt.grid(True)

Sxy= np.abs(Unnormlised_Peaks)*np.exp(1j*Phase_Differences)

N=len(Sxy)*10


lambda_=21e-2
b=12
image_1D = np.fft.ifftshift(np.fft.ifft(Sxy,N))
theta = np.linspace(-0.5, 0.5, N) * (lambda_/b)
# Plot magnitude of reconstructed 1D image
plt.figure(figsize=(10, 4))
plt.plot(np.rad2deg(theta), np.abs(image_1D))
plt.title("1D Reconstructed Image (Magnitude)")
plt.xlabel("Spatial Bin")
plt.ylabel("Amplitude")
plt.grid(True)

plt.show()
