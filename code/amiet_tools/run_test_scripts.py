#!/usr/bin/env python3
"""
Script to modify and run test scripts using the DARP2016 files
"""
import os
import subprocess
import shutil

# Create plots directory if it doesn't exist
plots_dir = "plots"
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Path to original scripts
script_paths = [
    "test_scripts/TestScript1_SingleGustPlots.py",
    "test_scripts/TestScript2_MultGustsPlots.py",
    "test_scripts/TestScript3_MultGustsSurfPressure.py"
]

# Directory for temporary scripts
temp_dir = "temp_scripts"
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

# Create modified versions of the scripts
for script_path in script_paths:
    script_name = os.path.basename(script_path)
    temp_script_path = f"{temp_dir}/{script_name}"
    
    # Read the original script
    with open(script_path, 'r') as f:
        content = f.read()
    
    # Disable LaTeX rendering
    content = content.replace("plt.rc('text', usetex=True)", "plt.rc('text', usetex=False)")
    
    # Enable saving figures and redirect to plots directory
    content = content.replace("save_fig = False", "save_fig = True")
    content = content.replace(".savefig('", ".savefig('plots/")
    
    # Change file format from eps to png to avoid LaTeX dependency
    content = content.replace(".eps')", ".png')")
    
    # Fix TestSetup initialization in the script
    if "AmT.loadTestSetup('../DARP2016_TestSetup.txt')" in content:
        content = content.replace(
            "AmT.loadTestSetup('../DARP2016_TestSetup.txt')",
            "AmT.TestSetup(c0=340.0, rho0=1.2, p_ref=20e-6, Ux=60.0, turb_intensity=0.025, length_scale=0.007, z_sl=-0.075)"
        )
    
    if "AmT.loadAirfoilGeom('../DARP2016_AirfoilGeom.txt')" in content:
        content = content.replace(
            "AmT.loadAirfoilGeom('../DARP2016_AirfoilGeom.txt')",
            "AmT.AirfoilGeom(b=0.075, d=0.225, Nx=100, Ny=101)"
        )
    
    # Write the modified script
    with open(temp_script_path, 'w') as f:
        f.write(content)
    
    print(f"Created modified script: {temp_script_path}")

# Run each modified script
for script_path in script_paths:
    script_name = os.path.basename(script_path)
    temp_script_path = f"{temp_dir}/{script_name}"
    
    print(f"Running {script_name}...")
    try:
        result = subprocess.run(['python', temp_script_path], check=True)
        print(f"Successfully completed {script_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")

# Clean up temp directory
shutil.rmtree(temp_dir)
print("All scripts completed. Check the plots directory for output files.") 