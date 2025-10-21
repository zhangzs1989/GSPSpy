A Python implementation of seismic phase detection using sweep signal correlation analysis for automatic seismic wave arrival time identification.

Project Overview
This project implements an advanced seismic signal processing algorithm that detects seismic wave phase arrivals through dual filter bank processing and cross-correlation analysis. The program can:

Generate linear sweep signals as source signals

Simulate multi-path seismic wave propagation and noise environments

Use adaptive filter banks for signal enhancement

Precisely detect phase arrival times through cross-correlation analysis

Visualize processing results and detection performance

Features
Dual Filter Bank Processing: Uses high and low frequency band filters for signal separation and enhancement

Phase Compensation: Eliminates phase distortion introduced by filters through reverse processing

Efficient Correlation Analysis: Utilizes FFT to accelerate cross-correlation calculations

Noise Robustness: Effective detection even under low signal-to-noise ratio conditions

Visualization Output: Provides intuitive display of data waveforms and detection results

Installation Requirements
System Requirements
Python 3.7+

Supported OS: Windows, Linux, macOS

Dependencies
bash
pip install numpy matplotlib scipy
Quick Installation
bash
git clone https://github.com/your-username/seismic-phase-scanning.git
cd seismic-phase-scanning
pip install -r requirements.txt
Usage
Basic Execution
python
python seismic_phase_detection.py
The program will automatically:

Generate sweep source signals

Create synthetic receiver data containing multiple arrivals

Add noise to simulate real environments

Execute phase scanning detection

Save results and generate visualizations

Parameter Configuration
Modify these key parameters in the code:

python
# Scanning parameters
spweepint = 10        # Scanning interval (samples)
sweeplen = 8000       # Total scanning length (samples)
fs = 100              # Sampling frequency (Hz)

# Filter parameters
d = 0.02              # Damping coefficient
dfo = 0.02            # Frequency offset

# Sweep signal parameters
f1 = 5                # Start frequency (Hz)
f2 = 10               # End frequency (Hz)
T = 300               # Sweep segment duration (seconds)
Output Files
The program generates:

.swp file: Contains time shifts and maximum correlation values

Visualization charts:

Receiver data waveform (marked with known arrival times)

Maximum correlation values vs. time shifts

Algorithm Principle
Signal Generation
Generate linear up-sweep and down-sweep chirp signals

Combine multiple sweep segments to form complete source signal

Noise Simulation
Add Gaussian white noise to simulate real observation environments

Supports configurable Signal-to-Noise Ratio (SNR)

Four-Stage Filter Processing
Forward Low-Frequency Filtering: Extract signals near (sweep_frequency - dfo)

Forward High-Frequency Filtering: Extract signals near (sweep_frequency + dfo)

Reverse Low-Frequency Filtering: Phase compensation processing

Reverse High-Frequency Filtering: Complete phase recovery

Correlation Detection
Use FFT-accelerated cross-correlation calculations

Detect maximum correlation values as phase arrival indicators

Result Interpretation
Output Chart Description
Figure 1: Receiver Data

Gray curve: Noisy receiver signal

Red dashed lines: Known phase arrival time positions

Figure 2: Correlation Detection Results

Red "+" markers: Maximum correlation values at each time shift

Peak positions: Detected potential phase arrival times

Performance Metrics
Program outputs at completion:

Total processing time

Number of processed time shifts

Detected peak positions and correlation values

Application Scenarios
Seismic monitoring and early warning systems

Geophysical exploration data processing

Academic research and teaching demonstrations

Seismic wave propagation characteristics research

Customization and Extension
Adding New Signal Types
Modify the xd generation section to use different source signals

Adjusting Detection Sensitivity
Modify threshold parameters for correlation peak detection:

python
peaks, properties = signal.find_peaks(ym, height=0.3*max(ym), distance=10)
Processing Real Data
Replace recdata with actual seismic observation data

Performance Optimization
For large-scale data processing, consider:

Increase spweepint to improve speed by reducing scanning density

Decrease sweeplen to narrow the scanning range

Use more efficient numerical computation libraries (e.g., optimized NumPy versions)

License
This project is licensed under the MIT License. See LICENSE file for details.

Contributing
Issues and Pull Requests are welcome to improve this project.

Citation
If you use this code in your research, please cite:

text
Global Seismic Phase Scanning Program. GitHub Repository. https://github.com/your-username/seismic-phase-scanning
Support
Please submit GitHub Issues for problems or contact the maintainer.
