"""
Script to plot a 2D contour of turbulent velocity spectrum Phi_2D

This script creates a 2D contour plot of the turbulent velocity spectrum
Phi_2D from the amiet_tools package and saves the resulting image.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import amiet_tools as AmT

# Disable LaTeX to avoid potential issues
plt.rc('text', usetex=False)
plt.close('all')

# Get absolute paths to the setup files
current_dir = os.path.dirname(os.path.abspath(__file__))
setup_file = os.path.join(current_dir, 'DARP2016_TestSetup.json')
airfoil_file = os.path.join(current_dir, 'DARP2016_AirfoilGeom.json')

# Load test setup from file
DARP2016Setup = AmT.loadTestSetup(setup_file)

# Export variables to current namespace
(c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
 flow_param, dipole_axis) = DARP2016Setup.export_values()

# Load airfoil geometry from file
DARP2016Airfoil = AmT.loadAirfoilGeom(airfoil_file)
(b, d, Nx, Ny, XYZ_airfoil, dx, dy) = DARP2016Airfoil.export_values()
XYZ_airfoil_calc = XYZ_airfoil.reshape(3, Nx*Ny)

# Frequency of analysis
kc = 20     # approx 7.2 kHz (chordwise normalised frequency = k0*(2*b))
f0 = kc*c0/(2*np.pi*(2*b))

FreqVars = AmT.FrequencyVars(f0, DARP2016Setup)
(k0, Kx, Ky_crit) = FreqVars.export_values()

# Get vector of spanwise hydrodynamic gusts 'Ky'
Ky = AmT.ky_vector(b, d, k0, Mach, beta, method='AcRad')

# Print key parameters
print(f"Flow velocity (Ux): {Ux} m/s")
print(f"Turbulence intensity: {turb_intensity}")
print(f"Integral length scale: {length_scale} m")
print(f"Frequency: {f0:.1f} Hz")
print(f"Streamwise wavenumber (Kx): {Kx:.2f} 1/m")
print(f"Critical spanwise wavenumber (Ky_crit): {Ky_crit:.2f} 1/m")
# Create wavenumber ranges centered on actual airfoil and flow parameters
num_points = 1000  # Number of points in each direction

# Extract extents directly from airfoil geometry
x_coords = XYZ_airfoil[0].flatten()
y_coords = XYZ_airfoil[1].flatten()
x_min, x_max = np.min(x_coords), np.max(x_coords)
y_min, y_max = np.min(y_coords), np.max(y_coords)
# Use existing Ky vector from AmT.ky_vector for more accuracy
# We already have Ky = AmT.ky_vector(b, d, k0, Mach, beta, method='AcRad')
# Let's center our kx range on Kx (aerodynamic wavenumber)
kx_half_range = max(3 * abs(Kx), 2 * np.max(np.abs(Ky)))
kx_range = np.linspace(Kx - kx_half_range, Kx + kx_half_range, num_points)
# Create a symmetric range for ky based on the Ky vector from AmT.ky_vector
ky_max_val = max(np.max(np.abs(Ky)), 2 * Ky_crit)
ky_range = np.linspace(-ky_max_val, ky_max_val, num_points)

# Create meshgrid for plotting
KX, KY = np.meshgrid(kx_range, ky_range)

# Calculate spectrum values for all points in the grid
print("Calculating 2D spectrum (this may take a moment)...")
print(f"Airfoil dimensions: x=[{x_min:.3f}, {x_max:.3f}], y=[{y_min:.3f}, {y_max:.3f}]")
print(f"Wavenumber ranges: kx=[{kx_range[0]:.1f}, {kx_range[-1]:.1f}], ky=[{ky_range[0]:.1f}, {ky_range[-1]:.1f}]")
print(f"Key wavenumbers: Kx={Kx:.2f}, Ky_crit={Ky_crit:.2f}")

# Calculate spectrum directly using AmT.Phi_2D function
# This is more efficient than calculating point by point
Phi_2D_values = AmT.Phi_2D(kx_range, ky_range, Ux, turb_intensity, length_scale, model='K')

print("Calculation complete! Creating plot...")

# Create plot
plt.figure(figsize=(10, 8))

# Use logarithmic color scale for better visualization
# Handle potential zero values in the spectrum
min_value = 1e-20  # Small positive value
Phi_2D_plot = np.maximum(Phi_2D_values, min_value)  # Replace zeros with small value

# Plot with logarithmic scale
contour = plt.pcolormesh(KX, KY, Phi_2D_plot.T, 
                         norm=LogNorm(),
                         cmap='viridis', shading='auto')

# Add colorbar
cbar = plt.colorbar(contour)
cbar.set_label('Φ(kₓ, kᵧ) [m²/s²]', fontsize=12)

# Add labels and title
plt.xlabel('Streamwise Wavenumber kₓ [1/m]', fontsize=12)
plt.ylabel('Spanwise Wavenumber kᵧ [1/m]', fontsize=12)
plt.title('2D Turbulent Velocity Spectrum Φ(kₓ, kᵧ)', fontsize=14)


# Set grid and make it pretty
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the figure with a more descriptive filename
output_dir = "spectrum_analysis"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
plot_filename = os.path.join(output_dir, f"phi_2D_contour_kc{kc}_f{f0:.0f}Hz.png")
plt.savefig(plot_filename, dpi=300)

print(f"Plot saved as '{plot_filename}'")

# Print detailed information about wavenumbers
print("\nWavenumber Information:")
print(f"Streamwise wavenumber (Kx): {Kx:.2f} 1/m")
print(f"Critical spanwise wavenumber (Ky_crit): {Ky_crit:.2f} 1/m")
print(f"kx range: [{kx_range[0]:.2f}, {kx_range[-1]:.2f}] 1/m with {len(kx_range)} points")
print(f"ky range: [{ky_range[0]:.2f}, {ky_range[-1]:.2f}] 1/m with {len(ky_range)} points")

# Print first few values of Ky vector from AmT.ky_vector
print("\nFirst 5 values of Ky vector from AmT.ky_vector:")
for i, k in enumerate(Ky[:5]):
    print(f"Ky[{i}] = {k:.2f} 1/m")

# Print statistics on the spectrum values
print("\nSpectrum Statistics:")
print(f"Maximum spectrum value: {np.max(Phi_2D_values):.6e}")
print(f"Minimum non-zero spectrum value: {np.min(Phi_2D_values[Phi_2D_values > 0]):.6e}")
print(f"Mean spectrum value: {np.mean(Phi_2D_values):.6e}")

plt.show() 