#!/usr/bin/env python3
"""
Script to run TestScript4_BeamformingSingleFreq.py and save results to whyresults folder
"""
import os
import subprocess
import shutil

# Create whyresults directory if it doesn't exist
results_dir = "whyresults"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Get current working directory
current_dir = os.getcwd()

# Path to original script
script_path = "test_scripts/TestScript4_BeamformingSingleFreq.py"
script_name = os.path.basename(script_path)

# Read the original script
with open(script_path, 'r') as f:
    content = f.read()

# Disable LaTeX rendering to avoid formatting issues
content = content.replace("plt.rc('text', usetex=True)", "plt.rc('text', usetex=False)")

# Enable saving figures
content = content.replace("save_fig = False", "save_fig = True")

# Specify the output directory for the figure
f0_regex = "f0 = kc*c0/(2*np.pi*(2*b))"
if f0_regex in content:
    # Add a line to create a string representation of f0 for filename
    content = content.replace(f0_regex, f0_regex + "\nf0_str = str(int(f0))")
    # Change the savefig line to include the path and use f0_str
    content = content.replace(
        "plt.savefig('AirfoilBeamf_' + f0 + 'Hz.png', dpi=200)",
        f"plt.savefig('{current_dir}/{results_dir}/AirfoilBeamf_' + f0_str + 'Hz.png', dpi=200)"
    )
else:
    # If regex doesn't match, use a more general replacement
    content = content.replace(
        "plt.savefig('AirfoilBeamf_' + f0 + 'Hz.png', dpi=200)",
        f"plt.savefig('{current_dir}/{results_dir}/AirfoilBeamf_' + str(int(f0)) + 'Hz.png', dpi=200)"
    )

# Fix the path to the DARP2016 files to use absolute paths
content = content.replace("'DARP2016_TestSetup.txt'", f"'{current_dir}/DARP2016_TestSetup.json'")
content = content.replace("'DARP2016_AirfoilGeom.txt'", f"'{current_dir}/DARP2016_AirfoilGeom.json'")

# Create a temporary modified script
temp_script_path = os.path.join(os.path.dirname(script_path), "temp_" + script_name)

# Write the modified script
with open(temp_script_path, 'w') as f:
    f.write(content)

print(f"Created modified script: {temp_script_path}")

# Run the modified script
print(f"\nRunning {script_name}...")
try:
    result = subprocess.run(['python', temp_script_path], 
                           check=True)
    
    print(f"Successfully completed {script_name}")
except subprocess.CalledProcessError as e:
    print(f"Error running {script_name}: {e}")

# Clean up the temporary script
os.remove(temp_script_path)

print(f"\nScript completed. Check the {results_dir} directory for output files.") 