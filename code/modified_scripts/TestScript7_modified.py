"""Airfoil turbulence analysis script."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import os

import amiet_tools as AmT

# Disable LaTeX rendering to avoid dependency issues
plt.rc('text', usetex=False)
plt.close('all')

# Create a new folder for all test case results (with timestamp to avoid overwriting)
import datetime
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
results_dir = f"../all_test_cases_results_{timestamp}"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Flag for saving figures - always true for this modified script
save_fig = True

# Flag for calculating and plotting near-field radiation slices (may take many minutes)
calc_nearfield = False

# Set better plot styles
plt.rcParams['font.size'] = 12
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Get absolute paths to the setup files
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(os.path.dirname(current_dir))
setup_file = os.path.join(root_dir, 'DARP2016_TestSetup.json')
airfoil_file = os.path.join(root_dir, 'DARP2016_AirfoilGeom.json')

# load test setup from file
DARP2016Setup = AmT.loadTestSetup(setup_file)

# export variables to current namespace
(c0, rho0, p_ref, Ux, turb_intensity, length_scale, z_sl, Mach, beta,
 flow_param, dipole_axis) = DARP2016Setup.export_values()

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# define airfoil points over the whole chord

# load airfoil geometry from file
DARP2016Airfoil = AmT.loadAirfoilGeom(airfoil_file)
(b, d, Nx, Ny, XYZ_airfoil, dx, dy) = DARP2016Airfoil.export_values()
XYZ_airfoil_calc = XYZ_airfoil.reshape(3, Nx*Ny)

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# frequency of operation
kc = 5                          # chordwise normalised frequency = k0*(2*b)
f0 = kc*c0/(2*np.pi*(2*b))      # approx 1.8 kHz

FreqVars = AmT.FrequencyVars(f0, DARP2016Setup)
(k0, Kx, Ky_crit) = FreqVars.export_values()

ac_wavelength = 2*np.pi/k0

# Create far field points for directivity
R_farfield = 50     # [m]
M_farfield = 181    # number of far-field mics in arc

theta_farfield = np.linspace(-np.pi/2, np.pi/2, M_farfield)
xy_farfield = R_farfield*np.sin(theta_farfield)
z_farfield = -R_farfield*np.cos(theta_farfield)

XZ_farfield = np.array([xy_farfield, np.zeros(xy_farfield.shape), z_farfield])
YZ_farfield = np.array([np.zeros(xy_farfield.shape), xy_farfield, z_farfield])

# Pressure field generated by the airfoil at the mesh
p_XZ_farfield = np.zeros(xy_farfield.shape, 'complex')
p_YZ_farfield = np.zeros(xy_farfield.shape, 'complex')

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# turbulence/gust amplitude
w0 = 1

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# First create a comprehensive plot showing response function magnitude along chord for multiple k_y values

# Define k_y values as multiples of Ky_crit
ky_multiples = np.array([0.0, 0.35, 0.75, 1.0, 1.25, 2.0])
ky_values = ky_multiples * Ky_crit
num_ky = len(ky_values)

# Define colors and line styles for the plot
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
line_styles = ['-', '--', '-.', '-.', ':', ':']

# Find points along the chord at mid-span (y=0)
mid_span_indices = np.where(np.isclose(XYZ_airfoil[1].flatten(), 0.0))[0]
xs_coords = XYZ_airfoil[0].flatten()[mid_span_indices]
# Sort by x-coordinate
sort_idx = np.argsort(xs_coords)
xs_sorted = xs_coords[sort_idx]

print(f"Number of points along mid-span: {len(xs_sorted)}")

# Create figure for chordwise magnitude plot
fig_chord_multi = plt.figure(figsize=(10, 6))
ax_chord_multi = fig_chord_multi.add_subplot(111)

# Calculate and plot response function magnitude for each k_y value
for i, ky in enumerate(ky_values):
    # Calculate pressure jump for this k_y value
    delta_p = AmT.delta_p(rho0, b, w0, Kx, ky, XYZ_airfoil[0:2], Mach)
    
    # Extract mid-span values and sort
    delta_p_midspan = delta_p.flatten()[mid_span_indices][sort_idx]
    
    # Calculate magnitude
    delta_p_mag = np.abs(delta_p_midspan)
    
    # Plot magnitude vs. normalized chord position
    ax_chord_multi.plot(xs_sorted/b, delta_p_mag, 
                        color=colors[i], linestyle=line_styles[i], linewidth=2,
                        label=f'$k_y = {ky_multiples[i]:.2f}k_y^{{crit}}$')

ax_chord_multi.set_xlim([-1, 1])  # Normalized chord limits
ax_chord_multi.set_ylim([0, 2])   # Limit y-axis to max value of 2
ax_chord_multi.set_xlabel('$x_s/b$', fontsize=14)
ax_chord_multi.set_ylabel('Magnitude', fontsize=14)
ax_chord_multi.set_title(f'$|g(x_s, k_y)|$', fontsize=16)
ax_chord_multi.legend(fontsize=12)
ax_chord_multi.grid(True, alpha=0.3)
plt.tight_layout()

# Save the figure
if save_fig:
    plt.savefig(f"{results_dir}/chordwise_magnitude_k0c_{kc}_multiple_ky.png", dpi=300)
    print(f"Saved chordwise magnitude plot to {results_dir}")

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Now create contour plots of the pressure jump real part for different k_y values

# Reshape airfoil coordinates for proper 2D plotting
x_mesh = XYZ_airfoil[0].reshape(Nx, Ny)
y_mesh = XYZ_airfoil[1].reshape(Nx, Ny)

# Create a figure for 2 rows of 3 contour plots (6 plots total)
fig_contours = plt.figure(figsize=(12, 10))
plt.subplots_adjust(wspace=0.4, hspace=0.3)

# Create a symmetric colormap
custom_cmap = plt.cm.RdBu_r
vmax = 0.1  # Maximum value for color scaling

# IMPORTANT: To make the plots smoother, we'll create a higher resolution grid and interpolate
# Create a higher resolution mesh for smoother plotting
x_interp = np.linspace(-b, b, 500)
y_interp = np.linspace(-d, d, 250)
X_interp, Y_interp = np.meshgrid(x_interp, y_interp)

# Loop through all k_y values to create contour plots
for i, ky in enumerate(ky_values):
    # Calculate pressure jump for this k_y
    delta_p = AmT.delta_p(rho0, b, w0, Kx, ky, XYZ_airfoil[0:2], Mach)
    delta_p_2d = delta_p.reshape(Nx, Ny)
    
    # Determine subplot position (2 rows, 3 columns)
    ax = fig_contours.add_subplot(2, 3, i+1)
    
    # Instead of directly plotting the irregular grid with pcolormesh,
    # we'll use pcolormesh on a regular grid after interpolation
    from scipy.interpolate import griddata
    
    # Convert mesh to points for interpolation input
    points = np.column_stack((x_mesh.flatten(), y_mesh.flatten()))
    values = np.real(delta_p_2d).flatten()
    
    # Interpolate onto regular grid
    interp_data = griddata(points, values, (X_interp, Y_interp), method='cubic', fill_value=0)
    
    # Plot with pcolormesh using the interpolated data
    contour = ax.pcolormesh(X_interp, Y_interp, interp_data, 
                           cmap=custom_cmap, shading='gouraud',
                           vmin=-vmax, vmax=vmax)
    
    ax.set_title(f'$k_y = {ky_multiples[i]:.2f}k_y^{{crit}}$', fontsize=14)
    ax.set_xlabel('$x_s$', fontsize=12)
    if i % 3 == 0:  # First column
        ax.set_ylabel('$y_s$', fontsize=12)
    
    # Set aspect ratio to be equal
    ax.set_aspect('equal')
    
    # Set axis limits more precisely
    ax.set_xlim([-b*1.1, b*1.1])
    ax.set_ylim([-d*1.1, d*1.1])

# Add a colorbar for the entire figure
cbar_ax = fig_contours.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig_contours.colorbar(contour, cax=cbar_ax)
cbar.set_label('Pressure Jump', fontsize=14)

# Set an overall title
fig_contours.suptitle(f'Real Part of Surface Pressure Jump ($k_0c = {kc}$, $M_x \\approx {Mach:.2f}$)', 
                     fontsize=16, y=0.98)

# Save the figure
if save_fig:
    plt.savefig(f"{results_dir}/pressure_jump_contours_k0c_{kc}.png", dpi=300)
    print(f"Saved pressure jump contours to {results_dir}")

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Run all four test cases and save in the results folder
for test_case in range(1, 5):
    if test_case == 1:
        # supercritical, normal incidence gust
        ky = 0
        fig_title = 'Kphi_000'
        desc = 'normal_incidence'

    elif test_case == 2:
        # supercritical gust, oblique incidence
        ky = 0.35*Ky_crit
        fig_title = 'Kphi_035'
        desc = 'oblique_incidence'

    elif test_case == 3:
        # supercritical gust, oblique incidence / close to critical
        ky = 0.75*Ky_crit
        fig_title = 'Kphi_075' 
        desc = 'near_critical'

    elif test_case == 4:
        # subcritical gust
        ky = 1.25*Ky_crit
        fig_title = 'Kphi_125'
        desc = 'subcritical'
    
    print(f"\nRunning test case {test_case}: {desc}")
    print(f"ky/Ky_crit = {ky/Ky_crit:.2f}")

    # Calculate the pressure 'jump' over the airfoil
    delta_p1 = AmT.delta_p(rho0, b, w0, Kx, ky, XYZ_airfoil[0:2], Mach)

    # reshape airfoil ac. source strengths, apply weights by grid area
    delta_p1_calc = (delta_p1*dx).reshape(Nx*Ny)*dy

    # Plot the airfoil source strength distribution
    re_p1_max = np.max(np.real(np.abs(delta_p1)))
    abs_p1_max = np.max(np.abs(delta_p1))

    fig_airfoil = plt.figure(figsize=(6.4, 5))
    ax1_airfoil = plt.subplot(121)
    re_p = ax1_airfoil.pcolormesh(XYZ_airfoil[0], XYZ_airfoil[1],
                                np.real(delta_p1)/re_p1_max, cmap='seismic',
                                shading='gouraud', vmin=-1, vmax=+1)  # Changed to gouraud for smooth appearance
    ax1_airfoil.axis('equal')
    ax1_airfoil.set_title(r'$Re\{\Delta p(x_s, y_s)\}$ [a.u.]', fontsize=18)

    ax1_airfoil.set_yticks([-d, 0, +d])
    ax1_airfoil.set_yticklabels([r'$-d$', r'$0$', r'$+d$'], fontsize=18)
    ax1_airfoil.set_xticks([-b, 0, b])
    ax1_airfoil.set_xticklabels([r'$-b$', r'$0$', r'$+b$'], fontsize=18)

    ax1_airfoil.set_ylim(-1.1*d, 1.1*d)
    ax1_airfoil.set_xlabel(r'$x_s$', fontsize=18)
    ax1_airfoil.set_ylabel(r'$y_s$', fontsize=18)
    plt.colorbar(re_p)

    ax2_airfoil = plt.subplot(122)
    abs_p = ax2_airfoil.pcolormesh(XYZ_airfoil[0], XYZ_airfoil[1],
                                20*np.log10(np.abs(delta_p1)/abs_p1_max),
                                cmap='inferno', shading='gouraud',  # Changed to gouraud for smooth appearance
                                vmax=0, vmin=-30)
    ax2_airfoil.axis('equal')
    ax2_airfoil.set_title(r'$|\Delta p(x_s, y_s)|$ [dB]', fontsize=18)

    ax2_airfoil.set_yticks([-d, 0, +d])
    ax2_airfoil.set_yticklabels([r'$-d$', r'$0$', r'$+d$'], fontsize=18)
    ax2_airfoil.set_xticks([-b, 0, b])
    ax2_airfoil.set_xticklabels([r'$-b$', r'$0$', r'$+b$'], fontsize=18)

    ax2_airfoil.set_ylim(-1.1*d, 1.1*d)
    ax2_airfoil.set_xlabel(r'$x_s$', fontsize=18)
    plt.colorbar(abs_p)
    fig_airfoil.set_tight_layout(True)

    # Save with both test case number and descriptive name for easy reference
    if save_fig:
        plt.savefig(f"{results_dir}/case{test_case}_DeltaP_{desc}_{fig_title}.png", dpi=300)
        print(f"Saved pressure jump case {test_case} to {results_dir}")

    # Calculate matrices of convected dipole Greens functions for each source to observers
    G_ffXZ = AmT.dipole3D(XYZ_airfoil_calc, XZ_farfield, k0, dipole_axis,
                        flow_param)

    G_ffYZ = AmT.dipole3D(XYZ_airfoil_calc, YZ_farfield, k0, dipole_axis,
                        flow_param)

    # Calculate the pressure in the far field
    p_XZ_farfield = G_ffXZ @ delta_p1_calc
    p_YZ_farfield = G_ffYZ @ delta_p1_calc

    # Plot the far field directivities [in dB]
    # normalise with respect to maximum FF pressure for parallel gust
    # (obtained from previous analysis)
    p_ff_max = 0.001139

    p_XZ_ff_norm = p_XZ_farfield/p_ff_max
    p_YZ_ff_norm = p_YZ_farfield/p_ff_max

    fig_dir_XZ = plt.figure(figsize=(6, 4))
    ax_dir_XZ = fig_dir_XZ.add_subplot(111, polar=True)
    plot_dir_XZ = ax_dir_XZ.plot(theta_farfield, 20*np.log10(np.abs(p_XZ_ff_norm)))
    ax_dir_XZ.set_thetamin(-90)
    ax_dir_XZ.set_thetamax(90)
    ax_dir_XZ.set_ylim([-40, 0])
    ax_dir_XZ.set_theta_zero_location('N')
    ax_dir_XZ.set_theta_direction('clockwise')
    ax_dir_XZ.set_thetagrids([-90, -45, 0, 45, 90],
                            labels=[r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$',
                                    r'$\theta = 0$', r'$+\frac{\pi}{4}$',
                                    r'$+\frac{\pi}{2}$'], size=18)
    ax_dir_XZ.set_rgrids([0., -10, -20, -30, -40],
                        labels=['0 dB', '-10', '-20', '-30', '-40'],
                        fontsize=12)

    # compensate axes position for half-circle plot
    ax_dir_XZ.set_position([0.1, -0.55, 0.8, 2])

    title_dir_XZ = ax_dir_XZ.set_title('Directivity: y=0 plane (phi=0)', 
                                      fontsize=18, pad=-55)

    if save_fig:
        fig_dir_XZ.savefig(f"{results_dir}/case{test_case}_dir_XZ_{desc}_{fig_title}.png", dpi=300)
        print(f"Saved XZ directivity case {test_case} to {results_dir}")

    fig_dir_YZ = plt.figure(figsize=(6, 4))
    ax_dir_YZ = fig_dir_YZ.add_subplot(111, polar=True)
    plot_dir_YZ = ax_dir_YZ.plot(theta_farfield, 20*np.log10(np.abs(p_YZ_ff_norm)))
    ax_dir_YZ.set_thetamin(-90)
    ax_dir_YZ.set_thetamax(90)
    ax_dir_YZ.set_ylim([-40, 0])
    ax_dir_YZ.set_theta_zero_location('N')
    ax_dir_YZ.set_theta_direction('clockwise')
    ax_dir_YZ.set_thetagrids([-90, -45, 0, 45, 90],
                            labels=[r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$',
                                    r'$\theta = 0$', r'$+\frac{\pi}{4}$',
                                    r'$+\frac{\pi}{2}$'], size=18)
    ax_dir_YZ.set_rgrids([0., -10, -20, -30, -40],
                        labels=['0 dB', '-10', '-20', '-30', '-40'],
                        fontsize=12)

    # compensate axes position for half-circle plot
    ax_dir_YZ.set_position([0.1, -0.55, 0.8, 2])

    title_dir_YZ = ax_dir_YZ.set_title('Directivity: x=0 plane (phi=pi/2)',
                                    fontsize=18, pad=-55)

    if save_fig:
        fig_dir_YZ.savefig(f"{results_dir}/case{test_case}_dir_YZ_{desc}_{fig_title}.png", dpi=300)
        print(f"Saved YZ directivity case {test_case} to {results_dir}")
        
    # Clear all figures before next iteration
    plt.close('all')

print("\nAll test cases and additional plots completed. Results saved to:", results_dir)

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Calculate the near-field radiation over the 2D cuts (may take many minutes)

if calc_nearfield:
    # create mesh for acoustic field 2D cuts
    mesh_side = 6*ac_wavelength  # [m]
    N_mesh = 201

    coord_vector = np.linspace(-mesh_side/2., mesh_side/2., N_mesh)
    X_mesh1, Z_mesh1 = np.meshgrid(coord_vector, coord_vector)
    Y_mesh1 = np.zeros(X_mesh1.shape)
    XZ_mesh1 = np.array([X_mesh1, Y_mesh1, Z_mesh1])

    Y_mesh2, Z_mesh2 = np.meshgrid(coord_vector, coord_vector)
    X_mesh2 = np.zeros(Y_mesh2.shape)
    YZ_mesh2 = np.array([X_mesh2, Y_mesh2, Z_mesh2])

    # Pressure field generated by the airfoil at the 2D mesh
    pressure_XZ_calc = np.zeros(X_mesh1.shape[0]*X_mesh1.shape[1], 'complex')
    pressure_YZ_calc = np.zeros(X_mesh2.shape[0]*X_mesh2.shape[1], 'complex')

    XZ_mesh1_calc = XZ_mesh1.reshape(3, XZ_mesh1.shape[1]*XZ_mesh1.shape[2])
    YZ_mesh2_calc = YZ_mesh2.reshape(3, YZ_mesh2.shape[1]*YZ_mesh2.shape[2])

    # Too many points in 2D cuts to generate entire G matrices; calculate the ac.
    # field of source grid point individually
    for s in range(delta_p1_calc.shape[0]):
        G_pXZ = AmT.dipole3D(XYZ_airfoil_calc[:, s, np.newaxis], XZ_mesh1_calc, k0,
                             dipole_axis, flow_param)

        G_pYZ = AmT.dipole3D(XYZ_airfoil_calc[:, s, np.newaxis], YZ_mesh2_calc, k0,
                             dipole_axis, flow_param)

        # Calculate the pressure in the near field
        pressure_XZ_calc += delta_p1_calc[s]*G_pXZ[:, 0]
        pressure_YZ_calc += delta_p1_calc[s]*G_pYZ[:, 0]

    # reshape nearfield meshes
    pressure_XZ = pressure_XZ_calc.reshape(XZ_mesh1[0].shape)
    pressure_YZ = pressure_YZ_calc.reshape(YZ_mesh2[0].shape)

    # find max pressure values for setting up color scale
    p_XZ_max = np.max(np.abs(np.real(pressure_XZ)))
    p_YZ_max = np.max(np.abs(np.real(pressure_YZ)))
    p_max = np.max((p_XZ_max, p_YZ_max))

    plt.figure(figsize=(6.3, 5.))
    plt.pcolormesh(X_mesh1/ac_wavelength, Z_mesh1/ac_wavelength,
                   np.real(pressure_XZ)/p_max, cmap='seismic',
                   shading='gouraud', vmin=-1.5, vmax=1.5)  # Changed to gouraud for smooth appearance
    plt.plot((-b/ac_wavelength, b/ac_wavelength), (0, 0), 'k', linewidth=6)
    plt.xlabel(r"$x/\lambda_0$", fontsize=18)
    plt.ylabel(r"$z/\lambda_0$", fontsize=18)
    plt.axis('equal')
    plt.xlim([-3, 3])
    plt.ylim([-3, 3])
    cbar1 = plt.colorbar()
    cbar1.set_label('Ac. Pressure [a.u.]', fontsize=15)
    plt.title(r"Acoustic Field on $y=0$ plane", fontsize=18)

    if save_fig:
        plt.savefig(f"{results_dir}/p_XZ_{fig_title}.png", dpi=300)
        print(f"Saved XZ acoustic field to {results_dir}")

    plt.figure(figsize=(6.3, 5.))
    plt.pcolormesh(Y_mesh2/ac_wavelength, Z_mesh2/ac_wavelength,
                   np.real(pressure_YZ)/p_max, cmap='seismic',
                   shading='gouraud', vmin=-1.5, vmax=1.5)  # Changed to gouraud for smooth appearance
    plt.plot((-d/ac_wavelength, d/ac_wavelength), (0, 0), 'k', linewidth=6)
    plt.xlabel(r"$y/\lambda_0$", fontsize=18)
    plt.ylabel(r"$z/\lambda_0$", fontsize=18)
    plt.axis('equal')
    plt.xlim([-3, 3])
    plt.ylim([-3, 3])

    cbar2 = plt.colorbar()
    cbar2.set_label('Ac. Pressure [a.u.]', fontsize=15)
    plt.title(r"Acoustic Field on $x=0$ plane", fontsize=18)

    if save_fig:
        plt.savefig(f"{results_dir}/p_YZ_{fig_title}.png", dpi=300)
        print(f"Saved YZ acoustic field to {results_dir}") 