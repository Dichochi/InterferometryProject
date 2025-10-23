import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
import astropy.units as u
import seaborn as sns
# --- Observer location: Cape Town ---
latitude = -33.933056
longitude = 18.476389
elevation = 25
location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=elevation*u.m)

# --- Time setup: 15:00 local SAST → 13:00 UTC ---
start_time_utc = '2025-10-03 13:00:00'
n_seconds = 780  # 780s observation
time = Time(start_time_utc) + np.arange(n_seconds) * u.s

# --- AltAz frame ---
altaz_frame = AltAz(obstime=time, location=location)

# --- Sun coordinates ---
sun = get_sun(time)
sun_altaz = sun.transform_to(altaz_frame)

# --- Store Sun alt/az in degrees ---
sun_data = np.zeros((n_seconds, 2))  # columns: [alt, az]
sun_data[:,0] = sun_altaz.alt.deg
sun_data[:,1] = sun_altaz.az.deg

# --- Convert Sun az/el to radians ---
az = np.radians(sun_data[:, 1])
el = np.radians(sun_data[:, 0])

# --- ENU unit vectors ---
E = np.cos(el) * np.sin(az)
N = np.cos(el) * np.cos(az)
U = np.sin(el)
enu_vectors = np.column_stack((E, N, U))

# --- Baseline ---
baseline = np.array([12, 15.97, 0])

# --- Projected baseline along Sun direction ---
path_length = np.dot(enu_vectors, baseline)
# --- Convert to phase (radians) ---
wavelength = 0.21  # meters
phase_shift = (2 * np.pi / wavelength) * path_length

# --- Wrap phase to [-pi, pi] ---
phase_shift_wrapped = (phase_shift + np.pi) % (2 * np.pi) - np.pi

# --- Count number of 2pi wraps using unwrapped phase ---
unwrapped_phase = np.unwrap(phase_shift_wrapped)
num_wraps = (unwrapped_phase[-1] - unwrapped_phase[0]) / (2*np.pi)
print("Total number of 2π wraps:", num_wraps)

# --- Time array in seconds for plotting ---
t = np.arange(n_seconds)

# --- Plot ---

fringe_rate = num_wraps / (n_seconds * u.s)  # wraps per second
print(f"Fringe rate: {fringe_rate:.4f} Hz")

plt.figure(figsize=(12, 4))
plt.plot(t, -phase_shift_wrapped)
plt.title("Measured Phase Shift After Correction vs Time")
plt.ylabel("Phase Shift [Radians]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.tight_layout()
plt.show()










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



fig, ax = plt.subplots(figsize=(10, 4))

# --- Plot Phase Differences in Degrees ---
ax.plot(t, np.rad2deg(-phase_shift_wrapped), color='#1f77b4', linewidth=1.2)

# --- Reference Line at Zero ---
ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.4)

# --- Axis Labels & Title ---
ax.set_xlabel(r"Time [s]", fontsize=12)
ax.set_ylabel(r"Fringe Phase, $\phi(t)$ [°]", fontsize=12)
ax.set_title(
    r"Simulated Solar Fringe Pattern over a 13-Minute Observation as Seen from the Telescope Array Location",
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
plt.savefig("Solar_Fringe_Pattern.png", dpi=400, bbox_inches='tight')
plt.show()



