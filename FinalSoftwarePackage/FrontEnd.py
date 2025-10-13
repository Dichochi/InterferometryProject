import tkinter as tk
from tkinter import filedialog, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import numpy as np
import SDRConnection
import tkinter.messagebox as messagebox
import threading
import DataAnalysis as danalysis
from PIL import Image, ImageTk, ImageEnhance # Pillow for better PNG support

class InterferometerApp:
    def __init__(self, master):
        self.master = master
        master.title("UCT Radio Astronomy Interferometer System")
        master.geometry("1500x900")

        # --- Main Layout: Left controls panel + Right plots panel ---
        self.left_frame = tk.Frame(master, width=250, bg="#f0f0f0", relief=tk.RAISED, borderwidth=2)
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        self.right_frame = tk.Frame(master)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        logo_img = Image.open("UCT_Logo.png")
        logo_img = logo_img.resize((70, 70), Image.Resampling.LANCZOS)
        self.logo_photo = ImageTk.PhotoImage(logo_img)
        logo_label = tk.Label(self.left_frame, image=self.logo_photo, bg="#f0f0f0")
        logo_label.pack(pady=2)

        self.status_label = tk.Label(self.left_frame, text="System Status", font=("Arial", 12, "bold"), bg="#f0f0f0")
        self.status_label.pack(pady=5)

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
        self.win_Resolution_var = tk.StringVar(value="Buffer Size: --")

        tk.Label(self.left_frame, textvariable=self.sample_rate_var, bg="#f0f0f0").pack(anchor="w", padx=10)
        tk.Label(self.left_frame, textvariable=self.duration_var, bg="#f0f0f0").pack(anchor="w", padx=10)
        tk.Label(self.left_frame, textvariable=self.win_Resolution_var, bg="#f0f0f0").pack(anchor="w", padx=10)
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

        # --- Channel 1 Gain ---
        tk.Label(self.left_frame, text="Channel 1 Gain [dB]:", bg="#f0f0f0").pack(anchor="w", padx=10, pady=5)
        self.gain_ch1 = ttk.Combobox(self.left_frame, values=[str(g) for g in range(0, 61, 10)])  # 0–60 dB in 1 dB steps
        self.gain_ch1.current(6)  # default 30 dB
        self.gain_ch1.pack(fill=tk.X, padx=20)

        # --- Channel 2 Gain ---
        tk.Label(self.left_frame, text="Channel 2 Gain [dB]:", bg="#f0f0f0").pack(anchor="w", padx=10, pady=5)
        self.gain_ch2 = ttk.Combobox(self.left_frame, values=[str(g) for g in range(0, 61, 10)])
        self.gain_ch2.current(6)  # default 30 dB
        self.gain_ch2.pack(fill=tk.X, padx=20)

        self.start_button = tk.Button(self.left_frame, text="Start Observation", command=self.start_observation)
        self.start_button.pack(fill=tk.X, padx=20, pady=5)

        observ_label = tk.Label(self.left_frame, text="Observation Data Analysis", font=("Arial", 12, "bold"), bg="#f0f0f0")
        observ_label.pack(pady=10)

        tk.Label(self.left_frame, text="Window Size (Resolution) in s:", bg="#f0f0f0").pack(anchor="w", padx=10, pady=5)
        self.win_size = ttk.Combobox(self.left_frame, values=["0.1", "1"])
        self.win_size.current(1)  # default 1024
        self.win_size.pack(fill=tk.X, padx=20)

        self.analyse_button = tk.Button(self.left_frame, text="Analyse Data", command=self.analyse_data)
        self.analyse_button.pack(fill=tk.X, padx=20, pady=5)



        # --- Create Figures for Plots (Right Panel) ---
        self.fig = Figure(figsize=(8, 8), dpi=100)

        # Top: Fringe Pattern
        self.ax_fringe = self.fig.add_subplot(3, 1, 1)
        self.ax_fringe.set_title("Fringe Pattern: Relative Phase Shift Between Antenna Channels")
        self.ax_fringe.set_xlabel("Time [s]")
        self.ax_fringe.set_ylabel("Radians")  # optional clarification

        # Middle: Correlation Coefficient
        self.ax_corr = self.fig.add_subplot(3, 1, 2)
        self.ax_corr.set_title("Correlation Coefficient (Magnitude of Normalized Cross-Correlation)")
        self.ax_corr.set_xlabel("Time [s]")
        self.ax_corr.set_ylabel("Correlation Coefficient")

        # Bottom: Cross-Correlation of 1-second chunk
        self.ax_crosscorr = self.fig.add_subplot(3, 1, 3)
        self.ax_crosscorr.set_title("Cross-Correlation of a 1-Second Chunk")
        self.ax_crosscorr.set_xlabel("Lag [s]")
        self.ax_crosscorr.set_ylabel("|Cxy| (Normalized Magnitude)")


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


    def update_plots(self):
        # Use the outputs from simulate_interferometer
        Cxy_window, lags_window, Phase_Differences, Correlation_Peaks,t = self.simulated_data  # stored after analyse_data

        # --- 1. Correlation Coefficient over full observation ---
        self.ax_corr.clear()
        self.ax_corr.plot(t,Correlation_Peaks)
        self.ax_corr.set_title("Correlation Coefficient (Magnitude of Normalized Cross-Correlation)")
        self.ax_corr.set_xlabel("Time [s]")
        self.ax_corr.set_ylabel("Correlation Coefficient")
        self.ax_corr.grid(True)

        # --- 2. Cross-Correlation of a single 1-second chunk ---
        self.ax_crosscorr.clear()
        self.ax_crosscorr.plot(lags_window*1e6, np.abs(Cxy_window))
        self.ax_crosscorr.set_title("Cross-Correlation of a 1-Second Chunk")
        self.ax_crosscorr.set_xlabel("Lag [μs]")
        self.ax_crosscorr.set_ylabel("|Cxy| (Normalized Magnitude)")
        self.ax_crosscorr.grid(True)

        # --- 3. Phase Differences over time (optional, for fringe) ---
        self.ax_fringe.clear()
        self.ax_fringe.plot(Phase_Differences)
        self.ax_fringe.set_title("Fringe Pattern: Relative Phase Shift Between Antenna Channels")
        self.ax_fringe.set_xlabel("Time [s]")
        self.ax_fringe.set_ylabel("Radians")  
        self.ax_fringe.grid(True)

        # Adjust layout and redraw
        self.fig.tight_layout(pad=2.0)
        self.canvas.draw()


    def start_observation(self):
        global GLOBAL_SAMPLE_RATE, GLOBAL_DURATION
        if self.pluto_status_var.get() != "Pluto: Connected!":
            messagebox.showwarning("Connection Error", "Connection Error! Connect to Pluto first")
        else:
            # Disable the start button to prevent multiple recordings
            self.start_button.config(state='disabled')
            
            self.acquisition_status_var.set("Data Acquisition In Progress")
            
            # Update system data display
            sample_rate_mhz = float(self.sample_rate_entry.get())
            duration_min = float(self.duration_entry.get())
            ch1_gain = int(self.gain_ch1.get())  # returns 0, 10, 20, ..., 60
            ch2_gain = int(self.gain_ch2.get())
            
            # Convert to standard units
            GLOBAL_SAMPLE_RATE = sample_rate_mhz * 1e6   # MHz → Hz
            GLOBAL_DURATION = duration_min * 60          # minutes → seconds
            
            # --- Update display ---
            self.sample_rate_var.set(f"Sample Rate: {sample_rate_mhz:.2f} MHz")
            self.duration_var.set(f"Duration: {duration_min:.1f} min")
            self.master.update_idletasks()
            
            def record_thread():
                try:
                    sdr_system.record_data(sample_rate=GLOBAL_SAMPLE_RATE, duration=GLOBAL_DURATION,ch1_gain=ch1_gain, ch2_gain=ch2_gain)
                    # Update GUI from the thread safely
                except Exception as e:
                    # Handle errors safely
                    self.master.after(0, lambda: self.acquisition_status_var.set("Recording Error"))
                    self.master.after(0, lambda: messagebox.showerror("Recording Error", f"Error: {str(e)}"))
                finally:
                    # Re-enable the button when done
                    self.master.after(0, lambda: self.acquisition_status_var.set("Data Acquisition Complete"))
                    self.master.after(0, lambda: self.start_button.config(state='normal'))
            
            thread = threading.Thread(target=record_thread, daemon=True)
            thread.start()

    def analyse_data(self):
        global GLOBAL_DURATION, GLOBAL_SAMPLE_RATE  # declare globals

        # Set default values if needed
        GLOBAL_DURATION = 13 * 60
        GLOBAL_SAMPLE_RATE = 2.5e6
        self.sample_rate_var.set(f"Sample Rate: {GLOBAL_SAMPLE_RATE*1e-6:.2f} MHz")
        self.duration_var.set(f"Duration: {GLOBAL_DURATION/60:.1f} min")
        self.win_Resolution_var.set(f"Buffer size: 1s")

        self.acquisition_status_var.set("Analysing Data")
        self.analyse_button.config(state='disabled')
        self.start_button.config(state='disabled')
        self.Pluto_Connect_button.config(state='disabled')

        buffer_duration = float(self.win_size.get())

        # --- Run analysis in a separate thread ---
        def analysis_thread():
            try:
                simulated_data = danalysis.simulate_interferometer(
                    duration=GLOBAL_DURATION,
                    sample_rate=GLOBAL_SAMPLE_RATE,
                    buffer_duration=buffer_duration
                )

                # Update GUI safely from main thread
                self.master.after(0, lambda: self.update_after_analysis(simulated_data))

            except Exception as e:
                self.master.after(0, lambda: messagebox.showerror(
                    "Analysis Error", f"Error during analysis: {str(e)}"))
            finally:
                # Re-enable button in main thread
                self.master.after(0, lambda: (
                self.analyse_button.config(state='normal'),
                self.start_button.config(state='normal'),
                self.Pluto_Connect_button.config(state='normal')
                ))


        thread = threading.Thread(target=analysis_thread, daemon=True)
        thread.start()


    def update_after_analysis(self, simulated_data):
        """Update plots and GUI after analysis is done."""
        self.simulated_data = simulated_data
        self.update_plots()
        self.acquisition_status_var.set("Analysis Complete")



if __name__ == "__main__":
    root = tk.Tk()
    app = InterferometerApp(root)
    root.mainloop()

