# 🌋 Seismic Phase Detection

Python implementation of seismic phase detection using sweep signal correlation analysis.

## Features

- **Dual Filter Bank Processing** - High/low frequency signal separation
- **Phase Compensation** - Reverse processing eliminates distortion  
- **Noise Robustness** - Effective detection under low SNR conditions
- **FFT-Accelerated Correlation** - Efficient computation for large datasets
- **Comprehensive Visualization** - Results interpretation and analysis

## Installation

```bash
git clone https://github.com/zhangzs1989/GSPSpy.
cd GSPSpy
pip install numpy matplotlib scipy
```
## Configuration
```bash
Key parameters in code:
spweepint = 10    # Scanning interval (samples)
sweeplen = 8000   # Total scanning length (samples)
fs = 100          # Sampling frequency (Hz)
d = 0.02          # Damping coefficient
f1, f2 = 5, 10    # Sweep frequency range (Hz)
```
## Quick Start
```bash
python seismic_phase_detection.py
```
The program will:
- Generate sweep source signals
- Create synthetic receiver data with multiple arrivals
- Add noise simulation (-10 dB SNR)
- Perform phase scanning detection
- Save results and generate visualizations
