"""
Global Seismic Phase Scanning Program
Python implementation of seismic phase detection using sweep signal correlation
This code performs seismic signal processing to detect phase arrivals through
dual filter bank processing and cross-correlation analysis.
autor:zs zhang
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

# Start timing for performance measurement
start_time = time.time()

# =============================================================================
# Parameter Configuration
# =============================================================================
spweepint = 100        # Scanning interval in samples
sweeplen = 4000        # Total scanning length in samples  
fs = 100               # Sampling frequency (Hz)
d = 0.02               # Damping coefficient
dfo = 0.02             # Frequency offset for filter banks

# Calculate filter parameters
A0 = 2 * np.pi * d / fs  # Normalized damping parameter
r = 1 / (1 + A0)         # Filter stability coefficient
a2 = r * r               # Squared coefficient for filter denominator

# =============================================================================
# Generate Sweep Signal (Source Signal)
# =============================================================================
f1 = 5                  # Start frequency (Hz)
f2 = 10                 # End frequency (Hz)
T = 300                 # Duration of each sweep segment (seconds)
N = fs * T              # Number of samples per segment
n = np.arange(N)        # Sample indices
t = n / fs              # Time vector

# Generate linear frequency sweep profiles
f01 = f1 + (f2 - f1) / (T) * t  # Upsweep: frequency increasing from f1 to f2
f02 = f2 + (f1 - f2) / (T) * t  # Downsweep: frequency decreasing from f2 to f1

# Generate sweep signals using frequency modulation
# The (2*pi*f)^2 term provides amplitude scaling based on frequency
x1 = (2 * np.pi * f01) ** 2 * np.cos(2 * np.pi * (f1 * t + (f2 - f1) / (2 * T) * t ** 2))
x2 = (2 * np.pi * f02) ** 2 * np.cos(2 * np.pi * (f2 * t + (f1 - f2) / (2 * T) * t ** 2))

# Combine multiple sweep segments to form complete source signal
xd = np.concatenate([x1, x2, x1, x2, x1, x2, x1, x2, x1, x2])

# =============================================================================
# Noise Generation Function
# =============================================================================
def noisegen(signal_data, snr_db):
    """
    Generate noisy signal with specified Signal-to-Noise Ratio (SNR)
    
    Parameters:
    signal_data: numpy array, original clean signal
    snr_db: float, signal-to-noise ratio in decibels (dB)
    
    Returns:
    noisy_signal: numpy array, signal with added noise
    noise: numpy array, generated noise component
    """
    # Calculate signal power (mean squared amplitude)
    signal_power = np.mean(signal_data ** 2)
    
    # Convert SNR from dB to linear scale and calculate noise power
    # SNR_linear = 10^(SNR_dB/10) = signal_power / noise_power
    noise_power = signal_power / (10 ** (snr_db / 10))
    
    # Generate Gaussian white noise with calculated power
    noise = np.random.normal(0, np.sqrt(noise_power), len(signal_data))
    
    # Add noise to original signal
    return signal_data + noise, noise

# =============================================================================
# Generate Synthetic Receiver Data
# =============================================================================
# Extended time vector for receiver data
t = np.arange(N * 5 * 2) / fs
buzero = 38 * fs  # Buffer length for zero padding

# Create receiver data with multiple seismic arrivals at different times
# Each arrival is a time-shifted and scaled version of the source signal
recdata = (np.concatenate([np.zeros(5 * fs), xd, np.zeros(buzero - 5 * fs)])) + \
          (np.concatenate([np.zeros(11 * fs), xd, np.zeros(buzero - 11 * fs)])) * 1 + \
          (np.concatenate([np.zeros(int(11.2 * fs)), xd, np.zeros(buzero - int(11.2 * fs))])) * 1.0 + \
          (np.concatenate([np.zeros(35 * fs), xd, np.zeros(buzero - 35 * fs)])) * 1.05 + \
          (np.concatenate([np.zeros(36 * fs), xd, np.zeros(buzero - 36 * fs)])) * 1.00

# Add noise to simulate real seismic data (-10 dB SNR)
Y, NOISE = noisegen(recdata, -10)
recdata = Y

# Extract source signal and frequency profile for processing
src = xd[:60000]                    # Source signal for correlation
yi_linear = np.concatenate([f01, f02])  # Complete frequency sweep profile

# =============================================================================
# Initialize Processing Variables
# =============================================================================
dd = 0               # Additional time shift parameter
xx = np.zeros(2)     # Input buffer for filter (current and previous samples)
yy6 = np.zeros(60003) # Output buffer for stage 1 filtering
j = 0                # Counter for processed shifts

ym = []              # List to store maximum correlation values
yy_storage = []      # List to store processed signals

# =============================================================================
# Main Scanning Loop - Process Different Time Shifts
# =============================================================================
print("Starting seismic phase scanning process...")

for k in range(0, sweeplen + 1, spweepint):
    """
    For each time shift k, apply the complete filter bank processing chain:
    1. Forward filtering through lower and higher frequency bands
    2. Reverse filtering for phase compensation  
    3. Cross-correlation with source signal
    """
    
    # Reset filter states for each iteration
    xx = np.zeros(2)
    yy6 = np.zeros(60003)
    
    # =========================================================================
    # Stage 1: Lower Frequency Band Filtering (Forward Direction)
    # =========================================================================
    # Purpose: Extract signal components in the lower frequency band
    # This enhances signals near (sweep_frequency - dfo)
    for ii in range(60000):
        # Get current frequency from sweep profile
        ff = yi_linear[ii]
        ffL = ff - dfo  # Lower band frequency
        
        # Calculate filter coefficients and parameters
        w1 = 2 * np.pi * ffL / fs  # Normalized angular frequency (lower band)
        w0 = 2 * np.pi * ff / fs   # Normalized angular frequency (center)
        a11 = -2 * r * np.cos(2 * np.pi * ffL / fs)  # Filter coefficient
        
        # Calculate gain normalization factors
        k1 = (1 - 2 * r * np.cos(w1) * np.cos(w0) + r * r * np.cos(2 * w0)) ** 2 + \
             (2 * r * np.sin(w0) * np.cos(w1) - r * r * np.sin(2 * w0)) ** 2
        k3 = (1 - np.cos(2 * w0)) ** 2 + (np.sin(2 * w0)) ** 2
        
        # Construct filter coefficients
        aa = [1, a11, a2]           # Denominator coefficients (AR part)
        bb1 = np.sqrt(k1 / k3)      # Gain normalization factor
        bb = [bb1, 0, -bb1]         # Numerator coefficients (MA part)
        
        # Update input buffer and apply IIR filter
        xx[1] = xx[0]  # Shift current to previous
        xx[0] = recdata[ii + dd + k]  # New input sample with time shift
        
        # IIR filter difference equation:
        # y[n] = b0*x[n] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
        yy6[ii + 3] = xx[0] * bb[0] + xx[1] * bb[2] - aa[1] * yy6[ii + 2] - aa[2] * yy6[ii + 1]
    
    # =========================================================================
    # Stage 2: Higher Frequency Band Filtering (Forward Direction)  
    # =========================================================================
    # Purpose: Extract signal components in the higher frequency band
    # This enhances signals near (sweep_frequency + dfo)
    xx = np.zeros(2)
    yy7 = np.zeros(60003)
    
    for ii in range(60000):
        ff = yi_linear[ii]
        ffH = ff + dfo  # Higher band frequency
        
        # Calculate filter parameters for higher band
        w2 = 2 * np.pi * ffH / fs
        w0 = 2 * np.pi * ff / fs
        
        # Calculate gain normalization
        k2 = (1 - 2 * r * np.cos(w2) * np.cos(w0) + r * r * np.cos(2 * w0)) ** 2 + \
             (2 * r * np.sin(w0) * np.cos(w2) - r * r * np.sin(2 * w0)) ** 2
        
        # Construct filter coefficients
        a21 = -2 * r * np.cos(2 * np.pi * ffH / fs)
        aa = [1, a21, a2]
        bb1 = np.sqrt(k2)
        bb = [bb1, 0, 0]  # Different coefficients than lower band
        
        # Apply IIR filter (input from previous stage)
        xx[0] = yy6[ii + 3]
        yy7[ii + 3] = xx[0] * bb[0] - aa[1] * yy7[ii + 2] - aa[2] * yy7[ii + 1]
    
    # =========================================================================
    # Phase Removal: Reverse Processing
    # =========================================================================
    # Purpose: Remove phase distortions introduced by forward filtering
    # This is done by processing the signal in reverse direction
    
    # Reverse the signal and frequency profile for backward processing
    YY7 = yy7[::-1]             # Reverse time order
    ryi_linear = yi_linear[::-1]  # Reverse frequency sweep profile
    
    # Stage 3: Lower Frequency Band Filtering (Reverse Direction)
    xx = np.zeros(2)
    yy8 = np.zeros(60003)
    
    for ii in range(60000):
        ff = ryi_linear[ii]  # Frequency from reversed profile
        ffL = ff - dfo
        
        # Calculate filter coefficients (same as Stage 1 but with reversed frequencies)
        w1 = 2 * np.pi * ffL / fs
        w0 = 2 * np.pi * ff / fs
        a11 = -2 * r * np.cos(2 * np.pi * ffL / fs)
        
        k1 = (1 - 2 * r * np.cos(w1) * np.cos(w0) + r * r * np.cos(2 * w0)) ** 2 + \
             (2 * r * np.sin(w0) * np.cos(w1) - r * r * np.sin(2 * w0)) ** 2
        k3 = (1 - np.cos(2 * w0)) ** 2 + (np.sin(2 * w0)) ** 2
        
        aa = [1, a11, a2]
        bb1 = np.sqrt(k1 / k3)
        bb = [bb1, 0, -bb1]
        
        # Apply IIR filter in reverse direction
        xx[1] = xx[0]
        xx[0] = YY7[ii]
        yy8[ii + 3] = xx[0] * bb[0] + xx[1] * bb[2] - aa[1] * yy8[ii + 2] - aa[2] * yy8[ii + 1]
    
    # Stage 4: Higher Frequency Band Filtering (Reverse Direction)
    xx = np.zeros(2)
    yy9 = np.zeros(60003)
    
    for ii in range(60000):
        ff = ryi_linear[ii]
        ffH = ff + dfo
        
        # Calculate filter coefficients (same as Stage 2 but with reversed frequencies)
        w2 = 2 * np.pi * ffH / fs
        w0 = 2 * np.pi * ff / fs
        
        k2 = (1 - 2 * r * np.cos(w2) * np.cos(w0) + r * r * np.cos(2 * w0)) ** 2 + \
             (2 * r * np.sin(w0) * np.cos(w2) - r * r * np.sin(2 * w0)) ** 2
        
        a21 = -2 * r * np.cos(2 * np.pi * ffH / fs)
        aa = [1, a21, a2]
        bb1 = np.sqrt(k2)
        bb = [bb1, 0, 0]
        
        # Apply final IIR filter
        xx[0] = yy8[ii + 3]
        yy9[ii + 3] = xx[0] * bb[0] - aa[1] * yy9[ii + 2] - aa[2] * yy9[ii + 1]
    
    # Store processed signal and prepare for correlation
    yy_storage.append(yy9)
    YY9 = yy9[::-1]  # Reverse back to original time order
    
    # =========================================================================
    # Cross-correlation Analysis
    # =========================================================================
    # Purpose: Detect seismic arrivals by measuring similarity between 
    # processed signal and source signal at different time lags
    
    # Compute cross-correlation using FFT for efficiency
    Ym = signal.correlate(YY9, src, mode='full')
    
    # Store maximum absolute correlation value
    # High values indicate strong similarity (potential seismic arrival)
    ym.append(np.max(np.abs(Ym)))
    j += 1
    
    # Progress reporting
    if k % 1000 == 0:
        print(f"Processed time shift: {k}/{sweeplen} samples")

print("Scanning process completed")

# =============================================================================
# Save Results to File
# =============================================================================
# Convert sample shifts to time in seconds
n = np.arange(0, sweeplen + 1, spweepint)
tm = n / 100  # Time in seconds

# Save correlation maxima to file for further analysis
filename = f'sweep_{spweepint}_{sweeplen}{d:.2f}.swp'
with open(filename, 'w') as fid:
    for i in range(len(tm)):
        # Save time shift and corresponding maximum correlation value
        fid.write(f'{tm[i]:.4f} {ym[i]:.4f}\n')

print(f"Results saved to: {filename}")

# =============================================================================
# Visualization of Results
# =============================================================================
print("Generating visualization plots...")

# Create figure with two subplots
fig = plt.figure(figsize=(15, 8))

# Subplot 1: Receiver data with arrival markers
ax1 = plt.subplot(2, 1, 1)
t_plot = np.arange(300000) / fs  # Time vector for plotting

# Plot receiver data
ax1.plot(t_plot, recdata[:300000], color=[156/255, 169/255, 194/255], 
         linewidth=0.8, label='Receiver Data')

# Mark known arrival times with vertical dashed lines
yylim = ax1.get_ylim()
arrival_times = [5, 11, 11.2, 35, 36]
for arr_time in arrival_times:
    ax1.plot([arr_time, arr_time], yylim, '--', 
             color=[194/255, 69/255, 83/255], linewidth=1.0, alpha=0.7)

# Configure plot appearance
ax1.set_xlim([0, 60])
ax1.text(0.8, 0.9, 'SNR = -10 dB', transform=ax1.transAxes, fontsize=14)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude')
ax1.grid(True, alpha=0.3)
ax1.set_title('Receiver Data with Known Arrival Times')

# Subplot 2: Correlation maxima vs time shift
ax2 = plt.subplot(2, 1, 2)
# Plot correlation results as red plus markers
ax2.plot(tm, ym, '+', color=[194/255, 69/255, 83/255], 
         markersize=6, label='Correlation Maxima')

ax2.set_xlim([0, 60])
ax2.set_xlabel('Time Shift (s)')
ax2.set_ylabel('$CC_{max}$')  # LaTeX formatting for subscript
ax2.grid(True, alpha=0.3)
ax2.set_title('Seismic Phase Detection - Correlation Maxima vs Time Shift')

# Improve layout and display plot
plt.tight_layout()
plt.show()

# =============================================================================
# Performance Summary
# =============================================================================
end_time = time.time()
execution_time = end_time - start_time

print("\n" + "="*50)
print("EXECUTION SUMMARY")
print("="*50)
print(f"Total processing time: {execution_time:.2f} seconds")
print(f"Number of time shifts processed: {len(ym)}")
print(f"Time shift range: 0 to {sweeplen/fs:.1f} seconds")
print(f"Maximum correlation value: {max(ym):.4f}")

# Simple peak detection to identify possible arrivals
# Note: This uses a fixed threshold - may need adjustment for different data
peaks, properties = signal.find_peaks(ym, height=0.3*max(ym), distance=10)
if len(peaks) > 0:
    print(f"Detected potential arrivals at: {tm[peaks]} seconds")
    print(f"Peak correlation values: {[ym[p] for p in peaks]}")
else:
    print("No clear arrivals detected above threshold")

print("Program completed successfully!")