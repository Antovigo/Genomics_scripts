#!/usr/bin/env python3

from pathlib import Path
import shutil
from tqdm import tqdm
import sys

# Define input and output paths
input_folder = Path("/home/anv888/Scratch/Genomics/Breseq_output_LTEE_ecoli")
output_folder = Path("/home/anv888/Scratch/Genomics/Data_LTEE_ecoli")

# Create required subfolders
subfolders = ["GD", "VCF", "HTML", "Coverage"]
for subfolder in subfolders:
    (output_folder / subfolder).mkdir(parents=True, exist_ok=True)

# Get all subfolders in input directory (sorted alphabetically)
input_subfolders = sorted([f for f in input_folder.iterdir() if f.is_dir()], 
                         key=lambda x: x.name)

# Define file mapping (source to destination)
def get_file_mapping(name):
    return [
        (input_folder / name / "output/output.vcf", output_folder / "VCF" / f"{name.name}.vcf"),
        (input_folder / name / "data/annotated.gd", output_folder / "GD" / f"{name.name}.gd"),
        (input_folder / name / "output/index.html", output_folder / "HTML" / f"{name.name}.html"),
        (input_folder / name / "08_mutation_identification/U00096.coverage.tab", output_folder / "Coverage" / f"{name.name}_U00096.coverage.tab"), # for e coli
        (input_folder / name / "08_mutation_identification/JX871397.coverage.tab", output_folder / "Coverage" / f"{name.name}_JX871397.coverage.tab") # for phage phi80
        # (input_folder / name / "08_mutation_identification/NZ_CP009977.coverage.tab", output_folder / "Coverage" / f"{name.name}_chr1.coverage.tab"), # for v natriegens
        # (input_folder / name / "08_mutation_identification/NZ_CP009978.coverage.tab", output_folder / "Coverage" / f"{name.name}_chr2.coverage.tab") # for v natriegens
    ]

# First, verify all files exist
print("Checking for required files...")
for name in input_subfolders:
    for src, _ in get_file_mapping(name):
        if not src.exists():
            print(f"Error: Required file not found: {src}")
            sys.exit(1)

# If we get here, all files exist. Proceed with copying
print("All required files found. Starting copy operation...")

# Calculate total number of files for progress bar
total_files = len(input_subfolders) * 4

# Copy files with progress bar
with tqdm(total=total_files, desc="Copying files") as pbar:
    for name in input_subfolders:
        pbar.set_description(f"Copying {name.name}")
        for src, dst in get_file_mapping(name):
            shutil.copy2(src, dst)
            pbar.update(1)

print("File copying completed successfully!")
