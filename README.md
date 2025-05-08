# Amiet Tools - Airfoil Turbulence Analysis

This repository contains code and results for analyzing turbulence-airfoil interaction using the amiet_tools Python package.

## Repository Structure

- **code/**
  - **amiet_tools_scripts/**: Original test scripts from the amiet_tools package
  - **modified_scripts/**: Modified versions of the test scripts with enhanced features
  
- **results/**
  - **pressure_plots/**: Plots showing surface pressure distributions, directivity patterns, and chordwise magnitude plots
  - **beamforming_plots/**: Beamforming analysis plots showing the acoustic source distribution
  - **endsem_results/**: Additional results from end-semester analysis
  - **spectral_analysis/**: Results from spectral analysis of turbulence
  - **velocity_analysis/**: Velocity field analysis results
  - **whyresults/**: Additional output visualizations

## Key Features

1. **Turbulence-Airfoil Interaction Analysis**:
   - Single gust plots showing pressure distributions
   - Multiple gust scenarios with varying wavenumbers
   - Directivity patterns for different frequencies
   
2. **Beamforming Analysis**:
   - 3D visualization of the measurement setup
   - Beamforming maps at multiple frequencies (kc = 5, 10, 20)
   - Visualizations with different rendering settings

3. **Visualization Improvements**:
   - Smooth pressure contour plots using interpolation
   - Multiple export formats (PNG, SVG, PDF) for presentations
   - Consistent color mapping and improved legends

## Dependencies

- Python 3.x
- NumPy
- Matplotlib
- amiet_tools package
- SciPy (for interpolation)

## References

This work is based on the amiet_tools package:
https://github.com/fchirono/amiet_tools 