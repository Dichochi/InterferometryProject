import tkinter as tk
from tkinter import filedialog, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import numpy as np
import SDRConnection
import tkinter.messagebox as messagebox
import threading
import DataAnalysis

class InterferometerApp:
    def __init__(self, master):
        self.master = master
        master.title("UCT Radio Astronomy Interferometer System")
        master.geometry("1200x800")

        # --- Main Layout: Left controls panel + Right plots panel ---
        self.left_frame = tk.Frame(master, width=250, bg="#f0f0f0", relief=tk.RAISED, borderwidth=2)
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        self.right_frame = tk.Frame(master)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

         # --- Status Block ---
        self.status_label = tk.Label(self.left_frame, text="System Status", font=("Arial", 12, "bold"), bg="#f0f0f0")
        self.status_label.pack(pady=10)

        self.pluto_status_var = tk.StringVar(value="Pluto: Disconnected")
        self.acquisition_status_var = tk.StringVar(value="Acquisition: Idle")  # <-- added

        self.pluto_status_label = tk.Label(self.left_frame,textvariable=self.pluto_status_var,bg="#f0f0f0",fg="red")
        self.pluto_status_label.pack(anchor="w", padx=10)
        tk.Label(self.left_frame, textvariable=self.acquisition_status_var, bg="#f0f0f0").pack(anchor="w", padx=10)

        # --- System Data Block ---
        sys_label = tk.Label(self.left_frame, text="System Data", font=("Arial", 12, "bold"), bg="#f0f0f0")
        sys_label.pack(pady=5)

        self.sample_rate_var = tk.StringVar(value="Sample Rate: --")
        self.duration_var = tk.StringVar(value="Duration: --")
        self.Freq_Resolution_var = tk.StringVar(value="Frequency Resolution: --")

        tk.Label(self.left_frame, textvariable=self.sample_rate_var, bg="#f0f0f0").pack(anchor="w", padx=10)
        tk.Label(self.left_frame, textvariable=self.duration_var, bg="#f0f0f0").pack(anchor="w", padx=10)
        tk.Label(self.left_frame, textvariable=self.Freq_Resolution_var, bg="#f0f0f0").pack(anchor="w", padx=10)
          # Example FFT channel selection

        # --- Controls Block ---



        ctrl_label = tk.Label(self.left_frame, text="Controls", font=("Arial", 12, "bold"), bg="#f0f0f0")
        ctrl_label.pack(pady=10)



        self.Pluto_Connect_button = tk.Button(self.left_frame, text="Connect to Pluto SDR", command=self.Pluto_Connect)
        self.Pluto_Connect_button.pack(fill=tk.X, padx=20, pady=5)


        observ_label = tk.Label(self.left_frame, text="Observation", font=("Arial", 12, "bold"), bg="#f0f0f0")
        observ_label.pack(pady=10)
         # Sample Rate Entry
        tk.Label(self.left_frame, text="Sample Rate (MHz):", bg="#f0f0f0").pack(anchor="w", padx=10)
        self.sample_rate_entry = tk.Entry(self.left_frame)
        self.sample_rate_entry.pack(fill=tk.X, padx=20, pady=2)
        self.sample_rate_entry.insert(0, "2.5")  # default value

        # Duration Entry
        tk.Label(self.left_frame, text="Duration (minutes):", bg="#f0f0f0").pack(anchor="w", padx=10)
        self.duration_entry = tk.Entry(self.left_frame)
        self.duration_entry.pack(fill=tk.X, padx=20, pady=2)
        self.duration_entry.insert(0, "1")  # default value

        self.start_button = tk.Button(self.left_frame, text="Start Observation", command=self.start_observation)
        self.start_button.pack(fill=tk.X, padx=20, pady=5)

        observ_label = tk.Label(self.left_frame, text="Observation Data Analysis", font=("Arial", 12, "bold"), bg="#f0f0f0")
        observ_label.pack(pady=10)

        tk.Label(self.left_frame, text="FFT Window Size:", bg="#f0f0f0").pack(anchor="w", padx=10, pady=5)
        self.fft_size = ttk.Combobox(self.left_frame, values=["256", "512", "1024", "2048"])
        self.fft_size.current(2)  # default 1024
        self.fft_size.pack(fill=tk.X, padx=20)

        self.observe_button = tk.Button(self.left_frame, text="Analyse Data", command=self.analyse_data)
        self.observe_button.pack(fill=tk.X, padx=20, pady=5)



        # --- Create Figures for Plots (Right Panel) ---
        self.fig = Figure(figsize=(8, 8), dpi=100)

        # Top: Fringe Pattern
        self.ax_fringe = self.fig.add_subplot(3, 1, 1)
        self.ax_fringe.set_title("Fringe Pattern")
        self.ax_fringe.set_xlabel("Time [s]")
        self.ax_fringe.set_ylabel("Amplitude")

        # Middle: Correlation Coefficient
        self.ax_corr = self.fig.add_subplot(3, 1, 2)
        self.ax_corr.set_title("Correlation Coefficient")
        self.ax_corr.set_xlabel("Time [s]")
        self.ax_corr.set_ylabel("r")

        # Bottom: 2 spectrograms side-by-side
        # Create a 2x2 grid and only use the bottom row for both spectrograms
        gs = self.fig.add_gridspec(3, 2, height_ratios=[1, 1, 2])

        self.ax_spec_ch1 = self.fig.add_subplot(gs[2, 0])
        self.ax_spec_ch1.set_title("Spectrogram – Channel 1")
        self.ax_spec_ch1.set_xlabel("Time [s]")
        self.ax_spec_ch1.set_ylabel("Frequency [Hz]")

        self.ax_spec_ch2 = self.fig.add_subplot(gs[2, 1])
        self.ax_spec_ch2.set_title("Spectrogram – Channel 2")
        self.ax_spec_ch2.set_xlabel("Time [s]")
        self.ax_spec_ch2.set_ylabel("Frequency [Hz]")

        # Adjust spacing
        self.fig.tight_layout(pad=4.0)
        

        # Attach to Tkinter canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # --- Placeholder Functions ---
    def Pluto_Connect(self):
        # Instantiate with defaults
        global sdr_system
        sdr_system = SDRConnection.SDRSystem()

        # Connect SDR
        if sdr_system.connect():
            # Record data
            self.pluto_status_var.set("Pluto: Connected!")
            # Update the label color (assumes you have a Label bound to this var)
             
            self.pluto_status_label.config(fg="green")
            self.acquisition_status_var.set("Mode:Idle")
        else:
            self.pluto_status_var.set("Pluto: Connection Failed")
            self.pluto_status_label.config(fg="red")


    def update_plots(self, data):
        # Correlator
        corr = np.correlate(data, data, mode="full")
        self.ax_corr.clear()
        self.ax_corr.plot(corr)
        self.ax_corr.set_title("Correlator Output")

        # Visibility Spectrum
        vis = np.abs(np.fft.fft(data))
        self.ax_vis.clear()
        self.ax_vis.plot(vis)
        self.ax_vis.set_title("Visibility Spectrum")

        # FFT Spectrogram
        self.ax_spec.clear()
        self.ax_spec.specgram(data, NFFT=int(self.fft_size.get()), Fs=1.0, noverlap=128)
        self.ax_spec.set_title("FFT Spectrogram")
        # Adjust layout spacing
        self.fig.tight_layout(pad=3.0)
        self.canvas.draw()

    def start_observation(self):
        if self.pluto_status_var.get() != "Pluto: Connected!":
            messagebox.showwarning("Connection Error", "Connection Error! Connect to Pluto first")
        else:
            # Disable the start button to prevent multiple recordings
            self.start_button.config(state='disabled')
            
            self.acquisition_status_var.set("Data Acquisition In Progress")
            
            # Update system data display
            sample_rate_mhz = float(self.sample_rate_entry.get())
            duration_min = float(self.duration_entry.get())
            
            # Convert to standard units
            sample_rate_var = sample_rate_mhz * 1e6  # MHz → Hz
            duration_var = duration_min * 60  # minutes → seconds
            
            self.sample_rate_var.set(f"Sample Rate: {sample_rate_var}")
            self.duration_var.set(f"Duration: {duration_var}s")
            self.master.update_idletasks()
            
            def record_thread():
                try:
                    sdr_system.record_data(sample_rate=sample_rate_var, duration=duration_var)
                    # Update GUI from the thread safely
                    self.master.after(0, lambda: self.acquisition_status_var.set("Data Acquisition Complete"))
                except Exception as e:
                    # Handle errors safely
                    self.master.after(0, lambda: self.acquisition_status_var.set("Recording Error"))
                    self.master.after(0, lambda: messagebox.showerror("Recording Error", f"Error: {str(e)}"))
                finally:
                    # Re-enable the button when done
                    self.master.after(0, lambda: self.start_button.config(state='normal'))
            
            thread = threading.Thread(target=record_thread, daemon=True)
            thread.start()

    def analyse_data(self):
        self.acquisition_status_var.set("Analysing Data")

        self.mode_status_var.set("Mode: Idle")


if __name__ == "__main__":
    root = tk.Tk()
    app = InterferometerApp(root)
    root.mainloop()

