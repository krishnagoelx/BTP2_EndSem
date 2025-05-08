#!/usr/bin/env python3
"""
Analysis of Turbulent Velocity Wavenumber Spectrum using amiet_tools package

This script generates plots to provide insights into the physics of turbulence spectra:
1. 2D contour (heat-map) of Φ(kₓ,kᵧ)
2. Effect of integral length scale Λ
3. Comparison of von Kármán vs Liepmann spectrum models
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import scipy.special as ss  # For gamma function

import amiet_tools as AmT

# Close any existing plots to ensure a clean environment
plt.close('all')

# Create results directory
results_dir = "velocity_analysis_revised"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# -----------------------------------------------------------------------
# Load the DARP2016 setup values - using the same values as TestScript2_MultGustsPlots.py
# -----------------------------------------------------------------------
# Using .json files which are definitely in the current directory
setup_file = 'DARP2016_TestSetup.json'
airfoil_file = 'DARP2016_AirfoilGeom.json'

DARP2016Setup = AmT.loadTestSetup(setup_file)
(c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
 flow_param, dipole_axis) = DARP2016Setup.export_values()

# Load airfoil geometry
DARP2016Airfoil = AmT.loadAirfoilGeom(airfoil_file)
(b, d, Nx, Ny, XYZ_airfoil, dx, dy) = DARP2016Airfoil.export_values()

# Frequency of analysis - same as in TestScript2_MultGustsPlots.py
kc = 20     # approx 7.2 kHz
f0 = kc*c0/(2*np.pi*(2*b))
FreqVars = AmT.FrequencyVars(f0, DARP2016Setup)
(k0, Kx, Ky_crit) = FreqVars.export_values()

# Print key parameters
print("DARP2016 Setup Values:")
print(f"Flow velocity (Ux): {Ux} m/s")
print(f"Turbulence intensity: {turb_intensity}")
print(f"Integral length scale: {length_scale} m")
print(f"Mach number: {Mach}")
print(f"Frequency: {f0:.1f} Hz")

# -----------------------------------------------------------------------
# Setup for plots
# -----------------------------------------------------------------------

# Energy-containing wavenumber
ke = (np.sqrt(np.pi)/length_scale)*(ss.gamma(5./6)/ss.gamma(1./3))
print(f"Energy-containing wavenumber (ke): {ke:.2f} 1/m")

# Create wavenumber ranges for 2D plots
num_points = 70  # Reduced for faster computation
kx_range = np.linspace(-5*ke, 5*ke, num_points)
ky_range = np.linspace(-5*ke, 5*ke, num_points)

# For comparison plots
kx_extended = np.linspace(-5*ke, 5*ke, 500)
ky_zero = np.array([0])  # For 1D slices with ky = 0

# -----------------------------------------------------------------------
# PLOT 1: 2D contour (heat-map) of Φ(kₓ,kᵧ)
# -----------------------------------------------------------------------

# Create 2D mesh grid for plotting
KX, KY = np.meshgrid(kx_range, ky_range)

# Calculate spectrum values - using same model as in TestScript2_MultGustsPlots.py
print("Calculating 2D spectrum (this may take a moment)...")

# Simplified calculation approach
Phi_2D_K = np.zeros((len(kx_range), len(ky_range)))

for i, kx in enumerate(kx_range):
    for j, ky in enumerate(ky_range):
        # Calculate spectrum value for this point
        Phi_2D_K[i, j] = AmT.Phi_2D(np.array([kx]), np.array([ky]), Ux, turb_intensity, length_scale, model='K')[0, 0]
    
    # Print progress every 10%
    if i % (num_points // 10) == 0:
        print(f"Progress: {i / num_points * 100:.0f}%")

print("Calculation complete! Creating plots...")

# Create plot
plt.figure(figsize=(10, 8))
# Use logarithmic color scale for better visualization
# Handle potential zero or negative values in the spectrum
Phi_2D_K_plot = np.maximum(Phi_2D_K.T, 1e-20)  # Ensure all values are positive for log scale
contour = plt.pcolormesh(KX, KY, Phi_2D_K_plot, 
                         cmap='viridis', shading='auto', norm=LogNorm())
plt.colorbar(contour, label='Φ(kₓ,kᵧ) [m²/s²]')
plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('Wavenumber kᵧ [1/m]')
plt.title('2D von Kármán Turbulent Velocity Spectrum Φ(kₓ,kᵧ)')
plt.axis('equal')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "turbulence_spectrum_2D_contour.png"), dpi=200)
plt.close()

# -----------------------------------------------------------------------
# PLOT 2: Effect of integral length scale Λ
# -----------------------------------------------------------------------

# Calculate spectra for different integral length scales
length_scale_1 = length_scale         # Original length scale 
length_scale_2 = length_scale * 2     # Double the length scale
length_scale_3 = length_scale * 0.5   # Half the length scale

ke_1 = (np.sqrt(np.pi)/length_scale_1)*(ss.gamma(5./6)/ss.gamma(1./3))
ke_2 = (np.sqrt(np.pi)/length_scale_2)*(ss.gamma(5./6)/ss.gamma(1./3))
ke_3 = (np.sqrt(np.pi)/length_scale_3)*(ss.gamma(5./6)/ss.gamma(1./3))

# Create a good wavenumber range to show the effect
k_range = np.logspace(np.log10(ke_3/10), np.log10(ke_3*10), 500)

Phi_L1 = AmT.Phi_2D(k_range, ky_zero, Ux, turb_intensity, length_scale_1, model='K')[:, 0]
Phi_L2 = AmT.Phi_2D(k_range, ky_zero, Ux, turb_intensity, length_scale_2, model='K')[:, 0]
Phi_L3 = AmT.Phi_2D(k_range, ky_zero, Ux, turb_intensity, length_scale_3, model='K')[:, 0]

plt.figure(figsize=(10, 6))
plt.loglog(k_range, Phi_L1, 'b-', linewidth=2, 
          label=f'Λ = {length_scale_1*1000:.1f} mm (Original)')
plt.loglog(k_range, Phi_L2, 'r-', linewidth=2, 
          label=f'Λ = {length_scale_2*1000:.1f} mm (2× Original)')
plt.loglog(k_range, Phi_L3, 'g-', linewidth=2, 
          label=f'Λ = {length_scale_3*1000:.1f} mm (0.5× Original)')

# Mark the energy-containing wavenumbers
plt.axvline(x=ke_1, color='b', linestyle='--', alpha=0.7, 
           label=f'kₑ = {ke_1:.1f} 1/m (Original)')
plt.axvline(x=ke_2, color='r', linestyle='--', alpha=0.7,
           label=f'kₑ = {ke_2:.1f} 1/m (2× Original)')
plt.axvline(x=ke_3, color='g', linestyle='--', alpha=0.7,
           label=f'kₑ = {ke_3:.1f} 1/m (0.5× Original)')

plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('Φ(kₓ, 0) [m²/s²]')
plt.title('Effect of Integral Length Scale on von Kármán Spectrum')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "turbulence_spectrum_length_scale_effect.png"), dpi=200)
plt.close()

# -----------------------------------------------------------------------
# PLOT 3: Comparison of von Kármán vs Liepmann models
# -----------------------------------------------------------------------

# Calculate spectra for both models along kx with ky=0
k_range = np.logspace(np.log10(ke/10), np.log10(ke*10), 500)
Phi_K = AmT.Phi_2D(k_range, ky_zero, Ux, turb_intensity, length_scale, model='K')[:, 0]
Phi_L = AmT.Phi_2D(k_range, ky_zero, Ux, turb_intensity, length_scale, model='L')[:, 0]

plt.figure(figsize=(10, 6))
plt.loglog(k_range, Phi_K, 'b-', linewidth=2, label='von Kármán Model')
plt.loglog(k_range, Phi_L, 'r-', linewidth=2, label='Liepmann Model')

# Create power law lines for comparison
# Reference point in high-k region for slope comparison
k_ref_idx = int(0.7 * len(k_range))
k_ref = k_range[k_ref_idx]
ref_K = Phi_K[k_ref_idx]
ref_L = Phi_L[k_ref_idx]

# Power law lines for comparison
k_slope = np.logspace(np.log10(k_ref), np.log10(k_range[-1]), 50)
k73_line = ref_K * (k_slope/k_ref)**(-7/3)  # -7/3 power law (von Kármán)
k52_line = ref_L * (k_slope/k_ref)**(-5/2)  # -5/2 power law (Modified slope)

plt.loglog(k_slope, k73_line, 'b--', linewidth=1.5, label='k^(-7/3) (von Kármán)')
plt.loglog(k_slope, k52_line, 'r--', linewidth=1.5, label='k^(-5/2) (Modified slope)')

plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('Φ(kₓ, 0) [m²/s²]')
plt.title('Comparison of von Kármán and Liepmann Spectrum Models')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "turbulence_spectrum_model_comparison.png"), dpi=200)
plt.close()

# Add an additional useful plot: Energy distribution across wavenumbers
plt.figure(figsize=(10, 6))
plt.semilogx(k_range, k_range*Phi_K, 'b-', linewidth=2, label='von Kármán Model')
plt.semilogx(k_range, k_range*Phi_L, 'r-', linewidth=2, label='Liepmann Model')
plt.axvline(x=ke, color='k', linestyle='--', alpha=0.7, 
           label=f'Energy-containing wavenumber kₑ = {ke:.1f} 1/m')
plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('kₓ·Φ(kₓ, 0) [m/s²]')
plt.title('Energy Distribution Across Wavenumbers')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "turbulence_energy_distribution.png"), dpi=200)
plt.close()

print(f"All plots have been saved to the {results_dir} directory.") 