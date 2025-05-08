"""Airfoil turbulence analysis script."""

import numpy as np
import os
import datetime

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.lines as mlines

import amiet_tools as AmT

plt.rc('text', usetex=False)
plt.close('all')

# Create output directory for plots with timestamp
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"../beamforming_results_{timestamp}"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

# Always save figures
save_fig = True

# Multiple frequencies to analyze
frequencies = [5, 10, 20]  # kc values (normalized frequencies)

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Get absolute paths to the setup files
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(current_dir)
setup_file = os.path.join(root_dir, 'DARP2016_TestSetup.json')
airfoil_file = os.path.join(root_dir, 'DARP2016_AirfoilGeom.json')

# load test setup from file
DARP2016Setup = AmT.loadTestSetup(setup_file)
# load airfoil geometry from file
DARP2016Airfoil = AmT.loadAirfoilGeom(airfoil_file)

# export variables to current namespace
(c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
 flow_param, dipole_axis) = DARP2016Setup.export_values()

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# define airfoil points over the whole chord

# load airfoil geometry from file
(b, d, Nx, Ny, XYZ_airfoil, dx, dy) = DARP2016Airfoil.export_values()
XYZ_airfoil_calc = XYZ_airfoil.reshape(3, Nx*Ny)


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# create DARP array

# DARP2016 spiral microphone array coordinates
XYZ_array, array_cal = AmT.DARP2016_MicArray()

# Number of mics
M = XYZ_array.shape[1]


# obtain propag time and shear layer crossing point for every source-mic pair
# (forward problem - frequency independent!)
T_sl_fwd, XYZ_sl_fwd = AmT.ShearLayer_matrix(XYZ_airfoil_calc, XYZ_array, z_sl, Ux, c0)

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Create grid of scan points

scan_sides = np.array([.65, .65])       # scan plane side length
scan_spacings = np.array([0.01, 0.01])  # scan points spacing

scan_xy = AmT.rect_grid(scan_sides, scan_spacings)

# Reshape grid points for 2D plotting
plotting_shape = (scan_sides/scan_spacings+1)[::-1].astype(int)
scan_x = scan_xy[0, :].reshape(plotting_shape)
scan_y = scan_xy[1, :].reshape(plotting_shape)

# Number of grid points
N = scan_xy.shape[1]

# create array with (x, y, z) coordinates of the scan points
scan_xyz = np.concatenate((scan_xy, np.zeros((1, N))))


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Plot the mics and grid points as 3D scatter plot

fig_grid = plt.figure(figsize=(10, 8))
ascan_x = fig_grid.add_subplot(111, projection='3d')
plot_grid1 = ascan_x.scatter(XYZ_array[0], XYZ_array[1], XYZ_array[2], c='r',
                             marker='o', s=50, label='Mic Array')
plot_grid2 = ascan_x.scatter(scan_xyz[0], scan_xyz[1], scan_xyz[2], c='b',
                             marker='^', s=20, alpha=0.5, label='Grid Points')

# Plot airfoil outline
airfoil_x = [-b, b, b, -b, -b]
airfoil_y = [-d, -d, d, d, -d]
airfoil_z = [0, 0, 0, 0, 0]
ascan_x.plot(airfoil_x, airfoil_y, airfoil_z, 'k-', linewidth=2)

ascan_x.set_xlabel('x [m]', fontsize=14)
ascan_x.set_ylabel('y [m]', fontsize=14)
ascan_x.set_zlabel('z [m]', fontsize=14)
ascan_x.set_title('Measurement Setup: Microphone Array and Scan Grid', fontsize=16)

# Create proxy artist to add legend
# --> numpoints = 1 to get only one dot in the legend
# --> linestyle= "none" So there is no line drawn in the legend
scatter1_proxy = mlines.Line2D([0], [0], linestyle="none", c='r', marker='o', markersize=10)
scatter2_proxy = mlines.Line2D([0], [0], linestyle="none", c='b', marker='^', markersize=10)
airfoil_proxy = mlines.Line2D([0], [0], linestyle="-", c='k', linewidth=2)
ascan_x.legend([scatter1_proxy, scatter2_proxy, airfoil_proxy], 
               ['Mic Array', 'Grid Points', 'Airfoil Outline'],
               numpoints=1, fontsize=12)

# Set better viewing angle
ascan_x.view_init(30, 225)
plt.tight_layout()

# Save in multiple formats
if save_fig:
    # Save as PNG (good for slides)
    plt.savefig(f"{output_dir}/3D_setup.png", dpi=300, bbox_inches='tight')
    # Save as SVG (vector format, good for scaling)
    plt.savefig(f"{output_dir}/3D_setup.svg", format='svg', bbox_inches='tight')
    # Save as PDF (vector format, good for documents)
    plt.savefig(f"{output_dir}/3D_setup.pdf", format='pdf', bbox_inches='tight')
    print(f"Saved 3D setup plot in multiple formats")


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Preparations for beamforming calculations

# Dynamic range for plotting
dynamic_range = 15      # [dB]

# obtain propag time and shear layer crossing point for every scan-mic pair
T_sl, XYZ_sl = AmT.ShearLayer_matrix(scan_xyz, XYZ_array, z_sl, Ux, c0)

# Loop through different frequencies
for kc in frequencies:
    # frequency [Hz]
    f0 = kc*c0/(2*np.pi*(2*b))
    freq_hz = int(np.round(f0))
    
    print(f"\nProcessing frequency: kc = {kc}, f = {freq_hz} Hz")
    
    FreqVars = AmT.FrequencyVars(f0, DARP2016Setup)
    (k0, Kx, Ky_crit) = FreqVars.export_values()
    
    # %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Calculate airfoil acoustic source strength CSM
    
    # vector of spanwise gust wavenumbers
    Ky = AmT.ky_vector(b, d, k0, Mach, beta)
    
    # Turbulence spectrum (von Karman)
    Phi2 = AmT.Phi_2D(Kx, Ky, Ux, turb_intensity, length_scale, model='K')[0]
    
    # calculate source CSM
    Sqq, Sqq_dxy = AmT.calc_airfoil_Sqq(DARP2016Setup, DARP2016Airfoil, FreqVars, Ky, Phi2)
    
    # apply weighting for airfoil grid areas
    Sqq *= Sqq_dxy
    
    # %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Create mic array CSM
    
    # create fwd transfer function
    G_fwd = AmT.dipole_shear(XYZ_airfoil_calc, XYZ_array, XYZ_sl_fwd, T_sl_fwd, k0, c0, Mach)
    
    # calculate mic array CSM
    CSM = (G_fwd @ Sqq @ G_fwd.conj().T)*4*np.pi
    
    # %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Apply classical beamforming algorithm
    
    # Creates steering vector and beamforming filters
    G_grid = np.zeros((M, N), 'complex')
    W = np.zeros((M, N), 'complex')
    
    # dipole grid with shear layer correction
    G_grid = AmT.dipole_shear(scan_xyz, XYZ_array, XYZ_sl, T_sl, k0, c0, Mach)
    
    # calculate beamforming filters
    for n in range(N):
        W[:, n] = G_grid[:, n]/(np.linalg.norm(G_grid[:, n], ord=2)**2)
    
    # vector of source powers
    A = np.zeros(N)
    
    # apply the beamforming algorithm
    for n in range(N):
        A[n] = (W[:, n].conj().T @ CSM @ W[:, n]).real
    
    # Reshape grid points for 2D plotting
    A_grid = np.zeros(plotting_shape, 'complex')
    A_grid = A.reshape(plotting_shape)
    
    # Good colormaps: viridis, inferno, plasma, magma
    colormap = 'plasma'
    
    # %%*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Plot the beamforming map
    
    fig1 = plt.figure(figsize=(10, 8))
    map1 = plt.pcolormesh(scan_x, scan_y, 10*np.log10(A_grid/A_grid.max()),
                          cmap=colormap, vmax=0, vmin=-dynamic_range,
                          shading='gouraud')  # Changed to 'gouraud' for smoother appearance
    map1.cmap.set_under('w')
    plt.title(f'Conventional Beamforming (kc = {kc}, f = {freq_hz} Hz)', fontsize=16)
    plt.axis('equal')
    plt.xlim([-scan_sides[0]/2, scan_sides[0]/2])
    plt.ylim([-scan_sides[1]/2, scan_sides[1]/2])
    plt.xlabel(r'$x$ [m]', fontsize=14)
    plt.ylabel(r'$y$ [m]', fontsize=14)
    
    cbar1 = fig1.colorbar(map1)
    cbar1.set_label('Normalised dB', fontsize=14)
    cbar1.ax.tick_params(labelsize=12)
    
    # Indicate the leading edge, trailing edge and sideplates on beamforming plot
    plt.vlines(-b, -d, d, color='k', linewidth=2)
    plt.vlines(b, -d, d, color='k', linewidth=2)
    plt.hlines(-d, -scan_sides[0]/2, scan_sides[0]/2, color='k', linewidth=2)
    plt.hlines(d, -scan_sides[0]/2, scan_sides[0]/2, color='k', linewidth=2)
    plt.text(-b+0.01, d-0.05, 'LE', fontsize=14, color='w', fontweight='bold')
    plt.text(b+0.01, d-0.05, 'TE', fontsize=14, color='w', fontweight='bold')
    
    plt.tight_layout()
    
    if save_fig:
        # Save in multiple formats
        plt.savefig(f"{output_dir}/AirfoilBeamf_kc{kc}_f{freq_hz}Hz.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/AirfoilBeamf_kc{kc}_f{freq_hz}Hz.svg", format='svg', bbox_inches='tight')
        plt.savefig(f"{output_dir}/AirfoilBeamf_kc{kc}_f{freq_hz}Hz.pdf", format='pdf', bbox_inches='tight')
        print(f"Saved beamforming map for kc = {kc}, f = {freq_hz} Hz")

print(f"\nAll plots have been saved to {output_dir}")
