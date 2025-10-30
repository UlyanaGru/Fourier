# Fourier Signal Processing Toolkit

A comprehensive Python-based signal processing toolkit for analyzing wave propagation in materials using Fourier analysis, correlation methods, and spectral decomposition. This tool is specifically designed for processing thickness measurement data and extracting wave characteristics.

## Features

- **Temporal Signal Analysis**: Visualize time-domain signals from multiple measurement points
- **Frequency Domain Analysis**: FFT-based spectral decomposition to identify dominant frequencies
- **Spatial Wave Analysis**: Extract dominant wavelengths from thickness profiles
- **Cross-Correlation**: Compute time delays between signals using correlation methods
- **Phase Velocity Calculation**: Determine wave propagation speed between measurement points
- **Customizable Visualization**: Professional-quality plots with customizable fonts and formatting
- **Multi-format Output**: Generate publication-ready figures in various formats

## Installation

### Requirements
```bash
pip install numpy pandas matplotlib scipy
```
<style>
.image-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 20px;
    justify-items: center;
    max-width: 1200px;
    margin: 0 auto;
    padding: 20px;
}

.image-grid img {
    width: 100%;
    max-width: 350px;
    height: 250px;
    object-fit: contain;
    border: 1px solid #ddd;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}
</style>

<div class="image-grid">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/mesh2_with_phase.png" alt="Mesh 2 with phase">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/mesh3_fklog10.png" alt="Mesh 3 FK log10">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/figout/freq_.jpg" alt="Frequency">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/mesh3_ampl.png" alt="Mesh 3 amplitude">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/figout/3Dsurf.png" alt="3D surface">
    <img src="https://github.com/UlyanaGru/Fourier/blob/master/mesh3_time_length_delta_filtred.png" alt="Time length delta filtered">
</div>