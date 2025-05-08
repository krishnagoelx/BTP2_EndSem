#!/usr/bin/env python3
"""
Script to remove author information and copyright notices from all Python files in the codebase.
"""

import os
import re
import glob

def clean_file(file_path):
    """Remove author information and copyright notices from a file."""
    print(f"Cleaning {file_path}")
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Remove docstring with copyright and author information
    content = re.sub(r'""".*?Copyright.*?"""', '"""Airfoil turbulence analysis script."""', content, flags=re.DOTALL)
    
    # Remove lines containing author, copyright, etc.
    lines = content.split('\n')
    cleaned_lines = []
    skip_mode = False
    
    for line in lines:
        if any(term in line.lower() for term in ["author:", "copyright", "@author", "fabio", "casagrande", "hirono", "fchirono"]):
            continue
        cleaned_lines.append(line)
    
    cleaned_content = '\n'.join(cleaned_lines)
    
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(cleaned_content)

def process_directory(directory):
    """Process all Python files in a directory and its subdirectories."""
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                clean_file(file_path)

# Process all Python files in the code directory
process_directory('code')

print("All Python files have been processed and copyright/author information has been removed.") 