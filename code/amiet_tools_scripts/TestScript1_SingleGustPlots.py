"""
amiet_tools - a Python package for turbulence-aerofoil noise prediction.
https://github.com/fchirono/amiet_tools
Copyright (c) 2020, Fabio Casagrande Hirono


TestScript1_SingleGustPlots.py

Test script 1: calculate interaction of a single turbulent gust with a flat
plate aerofoil, and plot:
    a) the surface pressure jump (real part and magnitude);
    b) the radiated acoustic field (real part) near the aerofoil over the
    x=0 and y=0 planes;
    c) the chordwise (y=0) and spanwise (x=0) far-field directivities (in dB).

The flow speed is assumed constant through all space (including far-field) -
i.e. there are no shear layer refraction effects.


This code was used to generate the Figures in Chap. 4, Section 4.1 of the
Authors' PhD thesis [Casagrande Hirono, 2018].


Author:
Fabio Casagrande Hirono
fchirono@gmail.com

"""


import numpy as np
import os

import amiet_tools as AmT

import matplotlib.pyplot as plt
plt.rc('text', usetex=False)
plt.close('all')


# flag for saving figures
save_fig = True
# folder to save figures
fig_folder = '../script1plots/'
# Create directory if it doesn't exist
os.makedirs(os.path.dirname(os.path.abspath(fig_folder)), exist_ok=True)

# Small value to avoid divide by zero in log10
eps = 1e-15

# flag for calculating and plotting near-field radiation slices
# --->>>may take many minutes!
calc_nearfield = False

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Get absolute paths to the setup files
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(current_dir)
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

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# turbulence/gust amplitude
w0 = 1

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# parallel incidence wavenumber component - select one case

test_case = 2       # int between 1 and 4

if test_case == 1:
    # supercritical, normal incidence gust
    ky = 0
    fig_title = 'Kphi_000'

elif test_case == 2:
    # supercritical gust, oblique incidence
    ky = 0.35*Ky_crit
    fig_title = 'Kphi_035'

elif test_case == 3:
    # supercritical gust, oblique incidence / close to critical
    ky = 0.75*Ky_crit
    fig_title = 'Kphi_075'

elif test_case == 4:
    # subcritical gust
    ky = 1.25*Ky_crit
    fig_title = 'Kphi_125'

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
ky_crit = k0/beta
# Calculate the pressure 'jump' over the airfoil
delta_p1 = AmT.delta_p(rho0, b, w0, Kx, ky, XYZ_airfoil[0:2], Mach)
# reshape airfoil ac. source strengths, apply weights by grid area
delta_p1_calc = (delta_p1*dx).reshape(Nx*Ny)*dy

# %%*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Plot the airfoil source strength distribution
re_p1_max = np.max(np.real(np.abs(delta_p1)))
abs_p1_max = np.max(np.abs(delta_p1))
# Calculate pressure jump along the chord at mid-span (y=0)
# Extract mid-span indices
mid_span_indices = np.where(np.isclose(XYZ_airfoil[1], 0))[0]
xs_midspan = XYZ_airfoil[0].flatten()[mid_span_indices]
delta_p1_midspan = delta_p1.flatten()[mid_span_indices]

# Sort by x-coordinate to ensure proper plotting
sort_idx = np.argsort(xs_midspan)
xs_midspan = xs_midspan[sort_idx]
delta_p1_midspan = delta_p1_midspan[sort_idx]

# Normalize x-coordinates by semi-chord
xs_normalized = xs_midspan / b

# Create figure for pressure jump magnitude along chord
fig_chord = plt.figure(figsize=(8, 6))
ax_chord = fig_chord.add_subplot(111)

# Plot current case
ax_chord.plot(xs_normalized, np.abs(delta_p1_midspan), 'b-', linewidth=2, 
              label=f'$|k_y/k_\\psi^{{(crit)}}| = {np.abs(ky/Ky_crit):.2f}$')

# Calculate and plot for additional ky values
# Supercritical case (|ky| < ky_crit)
ky_super = 0.5 * Ky_crit
delta_p_super = AmT.delta_p(rho0, b, w0, Kx, ky_super, XYZ_airfoil[0:2], Mach).flatten()[mid_span_indices][sort_idx]
ax_chord.plot(xs_normalized, np.abs(delta_p_super), 'r--', linewidth=2,
              label=f'$|k_y/k_\\psi^{{(crit)}}| = 0.50$')

# Near critical case (|ky| ≈ ky_crit)
ky_near = 0.95 * Ky_crit
delta_p_near = AmT.delta_p(rho0, b, w0, Kx, ky_near, XYZ_airfoil[0:2], Mach).flatten()[mid_span_indices][sort_idx]
ax_chord.plot(xs_normalized, np.abs(delta_p_near), 'g-.', linewidth=2,
              label=f'$|k_y/k_\\psi^{{(crit)}}| = 0.95$')

# Subcritical case (|ky| > ky_crit)
ky_sub = 1.5 * Ky_crit
delta_p_sub = AmT.delta_p(rho0, b, w0, Kx, ky_sub, XYZ_airfoil[0:2], Mach).flatten()[mid_span_indices][sort_idx]
ax_chord.plot(xs_normalized, np.abs(delta_p_sub), 'm:', linewidth=2,
              label=f'$|k_y/k_\\psi^{{(crit)}}| = 1.50$')

ax_chord.set_xlabel('$x_s/b$', fontsize=18)
ax_chord.set_ylabel('$|\\Delta p(x_s, y_s=0)|$', fontsize=18)
ax_chord.set_title(f'Pressure Jump Magnitude Along Chord ($k_0c = {kc}$)', fontsize=18)
ax_chord.grid(True)
ax_chord.legend(fontsize=14)

if save_fig:
    plt.savefig(fig_folder + 'DeltaP_chord_comparison_' + fig_title + '.png')

# Create a reshaped version of delta_p1 for the 2D plots
delta_p1_2d = delta_p1.reshape(Nx, Ny)

fig_airfoil = plt.figure(figsize=(6.4, 5))
ax1_airfoil = plt.subplot(121)

# Calculate meshgrid for proper pcolormesh plotting
x_mesh = XYZ_airfoil[0].reshape(Nx, Ny)
y_mesh = XYZ_airfoil[1].reshape(Nx, Ny)

re_p = ax1_airfoil.pcolormesh(x_mesh, y_mesh,
                              np.real(delta_p1_2d)/re_p1_max, cmap='seismic',
                              shading='nearest', vmin=-1, vmax=+1)
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
abs_p = ax2_airfoil.pcolormesh(x_mesh, y_mesh,
                               20*np.log10(np.abs(delta_p1_2d)/abs_p1_max + eps),
                               cmap='inferno', shading='nearest',
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

if save_fig:
    plt.savefig(fig_folder + 'DeltaP_real_mag_' + fig_title + '.png')

# %%*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

# %%*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Plot 2D surface map of δp(x,y) for specific gusts (supercritical and subcritical)

# Critical gust spanwise wavenumber
ky_crit = k0/beta

# Define a supercritical and subcritical ky value
ky_super = 0.5 * ky_crit  # Supercritical: |ky| < ky_crit
ky_sub = 1.5 * ky_crit    # Subcritical: |ky| > ky_crit

# Calculate pressure jump for supercritical gust
delta_p_super = AmT.delta_p(rho0, b, w0, Kx, ky_super, XYZ_airfoil[0:2], Mach)
delta_p_super_2d = delta_p_super.reshape(Nx, Ny)

# Calculate pressure jump for subcritical gust
delta_p_sub = AmT.delta_p(rho0, b, w0, Kx, ky_sub, XYZ_airfoil[0:2], Mach)
delta_p_sub_2d = delta_p_sub.reshape(Nx, Ny)

# Normalize by maximum values
abs_p_super_max = np.max(np.abs(delta_p_super_2d))
abs_p_sub_max = np.max(np.abs(delta_p_sub_2d))

# Create figure for the two gusts
fig_gusts = plt.figure(figsize=(12, 5))

# Plot supercritical gust
ax_super = plt.subplot(121)
contour_super = ax_super.pcolormesh(x_mesh, y_mesh,
                                   20*np.log10(np.abs(delta_p_super_2d)/abs_p_super_max + eps),
                                   cmap='inferno', shading='nearest',
                                   vmax=0, vmin=-30)
ax_super.axis('equal')
ax_super.set_title(r'Supercritical gust: $|k_y| < k_y^{crit}$', fontsize=16)
ax_super.set_xlabel(r'$x_s$', fontsize=14)
ax_super.set_ylabel(r'$y_s$', fontsize=14)
ax_super.set_yticks([-d, 0, +d])
ax_super.set_yticklabels([r'$-d$', r'$0$', r'$+d$'], fontsize=14)
ax_super.set_xticks([-b, 0, b])
ax_super.set_xticklabels([r'$-b$', r'$0$', r'$+b$'], fontsize=14)
ax_super.set_ylim(-1.1*d, 1.1*d)
cbar_super = plt.colorbar(contour_super, ax=ax_super)
cbar_super.set_label(r'$|\Delta p(x_s, y_s)|$ [dB]', fontsize=14)

# Plot subcritical gust
ax_sub = plt.subplot(122)
contour_sub = ax_sub.pcolormesh(x_mesh, y_mesh,
                               20*np.log10(np.abs(delta_p_sub_2d)/abs_p_sub_max + eps),
                               cmap='inferno', shading='nearest',
                               vmax=0, vmin=-30)
ax_sub.axis('equal')
ax_sub.set_title(r'Subcritical gust: $|k_y| > k_y^{crit}$', fontsize=16)
ax_sub.set_xlabel(r'$x_s$', fontsize=14)
ax_sub.set_yticks([-d, 0, +d])
ax_sub.set_yticklabels([r'$-d$', r'$0$', r'$+d$'], fontsize=14)
ax_sub.set_xticks([-b, 0, b])
ax_sub.set_xticklabels([r'$-b$', r'$0$', r'$+b$'], fontsize=14)
ax_sub.set_ylim(-1.1*d, 1.1*d)
cbar_sub = plt.colorbar(contour_sub, ax=ax_sub)
cbar_sub.set_label(r'$|\Delta p(x_s, y_s)|$ [dB]', fontsize=14)

fig_gusts.tight_layout()

# Add text showing the actual ky values used
ax_super.text(-0.9*b, -0.9*d, f'$k_y = {ky_super:.2f}$\n$k_y^{{crit}} = {ky_crit:.2f}$', 
              fontsize=12, bbox=dict(facecolor='white', alpha=0.7))
ax_sub.text(-0.9*b, -0.9*d, f'$k_y = {ky_sub:.2f}$\n$k_y^{{crit}} = {ky_crit:.2f}$', 
            fontsize=12, bbox=dict(facecolor='white', alpha=0.7))

if save_fig:
    plt.savefig(fig_folder + 'DeltaP_supercritical_subcritical_gusts.png')

# Plot airfoil directivity in x-plane and y-plane

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


# calculate matrices of convected dipole Greens functions for each source
# to observers
G_ffXZ = AmT.dipole3D(XYZ_airfoil_calc, XZ_farfield, k0, dipole_axis,
                      flow_param)

G_ffYZ = AmT.dipole3D(XYZ_airfoil_calc, YZ_farfield, k0, dipole_axis,
                      flow_param)

# Calculate the pressure in the far field
p_XZ_farfield = G_ffXZ @ delta_p1_calc
p_YZ_farfield = G_ffYZ @ delta_p1_calc

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# plot the far field directivities [in dB]

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

title_dir_XZ = ax_dir_XZ.set_title('Normalised Directivity on $y=0$ plane ($\\phi=0$)',
                                   fontsize=18, pad=-55)

if save_fig:
    plt.savefig(fig_folder + 'FF_PhiX_' + fig_title + '.png')


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

title_dir_YZ = ax_dir_YZ.set_title('Normalised Directivity on $x=0$ plane ($\\phi=\\pi/2$)',
                                   fontsize=18, pad=-55)

if save_fig:
    plt.savefig(fig_folder + 'FF_PhiY_' + fig_title + '.png')


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
                   shading='nearest', vmin=-1.5, vmax=1.5)
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
        plt.savefig(fig_folder + 'NF_real_X0_' + fig_title + '.png')

    plt.figure(figsize=(6.3, 5.))
    plt.pcolormesh(Y_mesh2/ac_wavelength, Z_mesh2/ac_wavelength,
                   np.real(pressure_YZ)/p_max, cmap='seismic',
                   shading='nearest', vmin=-1.5, vmax=1.5)
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
        plt.savefig(fig_folder + 'NF_real_Y0_' + fig_title + '.png')
