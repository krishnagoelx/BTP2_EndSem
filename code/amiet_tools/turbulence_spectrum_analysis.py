#!/usr/bin/env python3
"""
Analysis of Turbulent Velocity Wavenumber Spectrum using amiet_tools package

This script generates four plots to provide insights into the physics of turbulence spectra:
1. 2D contour (heat-map) of Φ(kₓ,kᵧ)
2. Radial cut (Φ vs |k|) to verify the -5/3 or -7/3 slope
3. Effect of integral length scale Λ
4. Comparison of von Kármán vs Liepmann spectrum models
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import scipy.special as ss  # Import scipy.special for gamma function

import amiet_tools as AmT

# Create results directory
results_dir = "whyresults"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# -----------------------------------------------------------------------
# Load the DARP2016 setup values
DARP2016Setup = AmT.loadTestSetup('DARP2016_TestSetup.json')
(c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
 flow_param, dipole_axis) = DARP2016Setup.export_values()

# Load airfoil geometry
DARP2016Airfoil = AmT.loadAirfoilGeom('DARP2016_AirfoilGeom.json')
(b, d, Nx, Ny, XYZ_airfoil, dx, dy) = DARP2016Airfoil.export_values()

print("DARP2016 Setup Values:")
print(f"Flow velocity (Ux): {Ux} m/s")
print(f"Turbulence intensity: {turb_intensity}")
print(f"Integral length scale: {length_scale} m")
print(f"Mach number: {Mach}")
print(f"Beta parameter: {beta}")

# -----------------------------------------------------------------------
# Setup for plots
# -----------------------------------------------------------------------

# Create a grid of wavenumbers for 2D spectrum plots
# Use normalized wavenumbers for better visualization
ke = (np.sqrt(np.pi)/length_scale)*(ss.gamma(5./6)/ss.gamma(1./3))
print(f"Energy-containing wavenumber (ke): {ke:.2f} 1/m")

# Create wavenumber ranges - roughly from 0.1*ke to 10*ke
kx_range = np.linspace(-5*ke, 5*ke, 200)
ky_range = np.linspace(-5*ke, 5*ke, 200)

# For 1D slice, create higher resolution in radial direction
k_radial = np.logspace(-1, 1, 100) * ke

# -----------------------------------------------------------------------
# PLOT 1: 2D contour (heat-map) of Φ(kₓ,kᵧ)
# -----------------------------------------------------------------------

# Create 2D mesh grid
KX, KY = np.meshgrid(kx_range, ky_range)
# Calculate spectrum values
Phi_2D_K = AmT.Phi_2D(kx_range, ky_range, Ux, turb_intensity, length_scale, model='K')

# Create plot
plt.figure(figsize=(10, 8))
# Use logarithmic color scale for better visualization of dynamic range
contour = plt.pcolormesh(KX/ke, KY/ke, np.log10(Phi_2D_K.T), 
                         cmap='viridis', shading='auto')
plt.colorbar(contour, label='log₁₀(Φ(kₓ,kᵧ)) [m²/s²]')
plt.xlabel('Normalized kₓ/kₑ')
plt.ylabel('Normalized kᵧ/kₑ')
plt.title('2D von Kármán Turbulent Velocity Spectrum Φ(kₓ,kᵧ)')
plt.axis('equal')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{results_dir}/turbulence_spectrum_2D_contour.png", dpi=200)

# -----------------------------------------------------------------------
# PLOT 2: Radial cut to verify -7/3 slope
# -----------------------------------------------------------------------

# Calculate spectrum along radial direction (kx=k, ky=0)
Phi_radial_K = AmT.Phi_2D(k_radial, np.array([0]), Ux, turb_intensity, length_scale, model='K')[:, 0]

# Create power law line for comparison
power_law = lambda k, k_ref: Phi_radial_K[np.abs(k_radial - k_ref).argmin()] * (k/k_ref)**(-7/3)
k_ref = ke  # Reference wavenumber for power law
k_line = np.logspace(np.log10(k_ref), np.log10(5*ke), 50)
power_line = power_law(k_line, k_ref)

plt.figure(figsize=(10, 6))
plt.loglog(k_radial/ke, Phi_radial_K, 'b-', linewidth=2, label='von Kármán Spectrum')
plt.loglog(k_line/ke, power_line, 'r--', linewidth=2, label=r'k$^{-7/3}$ Power Law')

plt.xlabel('Normalized Wavenumber k/kₑ')
plt.ylabel('Φ(k, 0) [m²/s²]')
plt.title('Radial Cut of von Kármán Spectrum with -7/3 Power Law')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{results_dir}/turbulence_spectrum_radial_cut.png", dpi=200)

# -----------------------------------------------------------------------
# PLOT 3: Effect of integral length scale Λ
# -----------------------------------------------------------------------

# Calculate spectra for different integral length scales
length_scale_1 = length_scale  # Original length scale from DARP2016
length_scale_2 = length_scale * 2  # Double the length scale

ke_1 = (np.sqrt(np.pi)/length_scale_1)*(ss.gamma(5./6)/ss.gamma(1./3))
ke_2 = (np.sqrt(np.pi)/length_scale_2)*(ss.gamma(5./6)/ss.gamma(1./3))

# Create wavenumber range to capture both spectra
kx_extended = np.linspace(-5*min(ke_1, ke_2), 5*max(ke_1, ke_2), 500)
ky_zero = np.array([0])  # Just look at ky = 0 slice

Phi_L1 = AmT.Phi_2D(kx_extended, ky_zero, Ux, turb_intensity, length_scale_1, model='K')[:, 0]
Phi_L2 = AmT.Phi_2D(kx_extended, ky_zero, Ux, turb_intensity, length_scale_2, model='K')[:, 0]

plt.figure(figsize=(10, 6))
plt.semilogx(kx_extended, Phi_L1, 'b-', linewidth=2, 
             label=f'Λ = {length_scale_1*1000:.1f} mm (Original)')
plt.semilogx(kx_extended, Phi_L2, 'r-', linewidth=2, 
             label=f'Λ = {length_scale_2*1000:.1f} mm (2× Original)')

# Mark the energy-containing wavenumbers
plt.axvline(x=ke_1, color='b', linestyle='--', alpha=0.7, 
            label=f'kₑ = {ke_1:.1f} 1/m (Original)')
plt.axvline(x=ke_2, color='r', linestyle='--', alpha=0.7,
            label=f'kₑ = {ke_2:.1f} 1/m (2× Original)')

plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('Φ(kₓ, 0) [m²/s²]')
plt.title('Effect of Integral Length Scale on von Kármán Spectrum')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{results_dir}/turbulence_spectrum_length_scale_effect.png", dpi=200)

# -----------------------------------------------------------------------
# PLOT 4: Comparison of von Kármán vs Liepmann models
# -----------------------------------------------------------------------

# Calculate spectra for both models along kx with ky=0
Phi_K = AmT.Phi_2D(kx_extended, ky_zero, Ux, turb_intensity, length_scale, model='K')[:, 0]
Phi_L = AmT.Phi_2D(kx_extended, ky_zero, Ux, turb_intensity, length_scale, model='L')[:, 0]

plt.figure(figsize=(10, 6))
plt.loglog(kx_extended, Phi_K, 'b-', linewidth=2, label='von Kármán Model')
plt.loglog(kx_extended, Phi_L, 'r-', linewidth=2, label='Liepmann Model')

# Power law lines for comparison
k_range = np.logspace(np.log10(ke), np.log10(5*ke), 50)
k73_line = power_law(k_range, ke)  # -7/3 power law
k5_line = k73_line[0] * (k_range/ke)**(-5)  # -5 power law (approximate Liepmann roll-off)

plt.loglog(k_range, k73_line, 'b--', linewidth=1.5, label=r'k$^{-7/3}$ (von Kármán)')
plt.loglog(k_range, k5_line, 'r--', linewidth=1.5, label=r'k$^{-5}$ (Liepmann approximation)')

plt.xlabel('Wavenumber kₓ [1/m]')
plt.ylabel('Φ(kₓ, 0) [m²/s²]')
plt.title('Comparison of von Kármán and Liepmann Spectrum Models')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{results_dir}/turbulence_spectrum_model_comparison.png", dpi=200)

print(f"All plots have been saved to the {results_dir} directory.") 