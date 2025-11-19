#!/usr/bin/env python3
"""
Complete pipeline for kinase conformational state analysis.
Steps:
1. Add chain IDs to all PDB files
2. Run kincore analysis on each file (parallelized)
3. Process output to generate conformational state distribution
4. Create structure mapping file for docking analysis

Caution: have to activate kincore
"""

import os
import subprocess
import sys
import re
import multiprocessing as mp
from collections import defaultdict
from pathlib import Path

# Configuration
RECEPTOR_BASE_DIR = "/data/work/dock/virtual_screening/input/receptors"
CHAIN_ID = "A"
OUTPUT_DIR = "/data/work/dock/analisis/data"
KINCORE_OUTPUT_FILE = "kincore_raw_output.txt"
FINAL_COUNTS_FILE = "conf_state.txt"

# Multiprocessing configuration
#NUM_CORES = mp.cpu_count() - 1  # Leave one core for system (recommended)
# Uncomment and modify to set custom number of cores:
NUM_CORES = 8  # Use 8 cores

def extract_structure_id_from_filename(filename):
    """Extract structure_id from receptor_XXX_chainA.pdb format"""
    match = re.search(r'receptor_(\d+)_chainA\.pdb', filename)
    if match:
        return int(match.group(1))
    return None

def add_chain_id_to_pdb(input_file, chain_id="A"):
    """Add chain ID to PDB file and save as new file."""
    if not os.path.exists(input_file):
        print(f"File {input_file} does not exist.")
        return None
    
    # Create output filename
    base_name = input_file.replace('.pdb', '')
    output_file = f"{base_name}_chainA.pdb"
    
    # Skip if already exists
    if os.path.exists(output_file):
        print(f"Using existing: {output_file}")
        return output_file
    
    # Read the input file
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
    
    # Write to new file with chain ID added
    modified_lines = 0
    with open(output_file, 'w') as outfile:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Check if chain identifier is already present (22nd column)
                if len(line) > 21 and line[21] != chain_id:
                    # Insert the chain identifier at the 22nd column
                    modified_line = line[:21] + chain_id + line[22:]
                    outfile.write(modified_line)
                    modified_lines += 1
                else:
                    # If the chain ID is already correct, write the line as is
                    outfile.write(line)
            else:
                # Write non-atom lines as is
                outfile.write(line)
    
    print(f"Created: {output_file}")
    return output_file

def kincore_worker(file_info):
    """Worker function for multiprocessing kincore analysis."""
    chainA_file, kinase_name, structure_id = file_info
    
    try:
        # Run kincore command using the actual Python script
        result = subprocess.run(['python', '/data/work/Kincore-standalone2/kinase_state.py', chainA_file], 
                              capture_output=True, 
                              text=True, 
                              check=True)
        
        # Format output with kinase name and structure_id prefix
        output_lines = []
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if line.strip():
                output_lines.append(f"KINASE:{kinase_name}\tSTRUCTURE_ID:{structure_id}\t{line}")
        
        return {
            'success': True,
            'file': chainA_file,
            'kinase': kinase_name,
            'structure_id': structure_id,
            'output': '\n'.join(output_lines) + '\n' if output_lines else '',
            'error': None
        }
        
    except subprocess.CalledProcessError as e:
        return {
            'success': False,
            'file': chainA_file,
            'kinase': kinase_name,
            'structure_id': structure_id,
            'output': '',
            'error': f"CalledProcessError: {e.stderr}"
        }
    except Exception as e:
        return {
            'success': False,
            'file': chainA_file,
            'kinase': kinase_name,
            'structure_id': structure_id,
            'output': '',
            'error': f"Exception: {str(e)}"
        }

def process_conformational_states(filename, output_file="conformational_state_counts.txt"):
    """Process kincore output to count conformational states overall and per kinase."""
    # Define expected combinations in desired order
    order = [
        "DFGin-ABAminus", "DFGin-BLAminus", "DFGin-BLAplus",
        "DFGin-BLBminus", "DFGin-BLBplus", "DFGin-BLBtrans",
        "DFGin-Unassigned", "DFGinter-BABtrans", "DFGinter-Unassigned",
        "DFGout-BBAminus", "DFGout-Unassigned", "Unassigned-Unassigned"
    ]
    
    # Initialize counts and file tracking
    overall_counts = defaultdict(int)
    kinase_counts = defaultdict(lambda: defaultdict(int))
    
    # Track which files belong to each state
    state_files = defaultdict(list)  # Overall: state -> [files]
    kinase_state_files = defaultdict(lambda: defaultdict(list))  # Per kinase: kinase -> state -> [files]
    
    # NEW: Track structure mapping for docking analysis
    structure_mapping = []  # List of (kinase, structure_id, filename, conformational_state)
    
    # Track processed files to avoid duplicates (since kincore outputs 5 lines per structure)
    processed_files = set()
    
    for item in order:
        overall_counts[item] = 0
    
    # Read lines
    with open(filename, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    # Process lines with kinase information
    unique_structures = 0
    validation_errors = []
    
    for line in lines:
        if line.startswith("KINASE:"):
            # Parse kinase, structure_id, and kincore output
            parts = line.split('\t')
            if len(parts) >= 3:
                kinase_info = parts[0]
                structure_id_info = parts[1]
                kincore_line = parts[2]
                
                # Extract kinase name and structure_id
                kinase_name = kinase_info.replace("KINASE:", "")
                structure_id = int(structure_id_info.replace("STRUCTURE_ID:", ""))
                
                # Process kincore output
                kincore_parts = kincore_line.split()
                if len(kincore_parts) >= 4:
                    pdb_filename = kincore_parts[0]  # First column is PDB filename
                    
                    # Create unique identifier for this structure
                    structure_key = f"{kinase_name}_{structure_id}"
                    
                    # Skip if we've already processed this structure
                    if structure_key in processed_files:
                        continue
                    
                    processed_files.add(structure_key)
                    
                    # Validate structure_id matches filename
                    expected_filename = f"receptor_{structure_id:03d}_chainA.pdb"
                    if pdb_filename != expected_filename:
                        validation_errors.append(
                            f"WARNING: {kinase_name} structure_id {structure_id} "
                            f"filename mismatch: expected {expected_filename}, got {pdb_filename}"
                        )
                    
                    fourth_column = kincore_parts[3]
                    split_forth = fourth_column.split('_')
                    
                    if len(split_forth) >= 3:
                        second_part = split_forth[1]  # DFG state
                        third_part = split_forth[2]   # Conformation
                        
                        # Handle None cases
                        if second_part == "None":
                            result = "Unassigned-Unassigned"
                        elif third_part == "None":
                            result = f"{second_part}-Unassigned"
                        else:
                            result = f"{second_part}-{third_part}"
                    else:
                        result = "Unassigned-Unassigned"
                        
                    # Count for overall and per kinase
                    if result in overall_counts:
                        overall_counts[result] += 1
                        kinase_counts[kinase_name][result] += 1
                    else:
                        overall_counts["Unassigned-Unassigned"] += 1
                        kinase_counts[kinase_name]["Unassigned-Unassigned"] += 1
                        result = "Unassigned-Unassigned"
                    
                    # Track files for each state
                    state_files[result].append(f"{kinase_name}/{pdb_filename}")
                    kinase_state_files[kinase_name][result].append(pdb_filename)
                    
                    # NEW: Add to structure mapping
                    structure_mapping.append((kinase_name, structure_id, pdb_filename, result))
                        
                    # Initialize kinase counts for all states
                    for item in order:
                        if item not in kinase_counts[kinase_name]:
                            kinase_counts[kinase_name][item] = 0
                            
                    unique_structures += 1
    
    # Print validation errors
    if validation_errors:
        print("\n=== VALIDATION WARNINGS ===")
        for error in validation_errors:
            print(error)
        print("=== END WARNINGS ===\n")
    
    # Save overall results
    with open(output_file, 'w', encoding='utf-8') as f_out:
        f_out.write(f"Overall Kinase Conformational State Distribution\n")
        f_out.write(f"Total unique structures analyzed: {unique_structures}\n")
        f_out.write(f"="*50 + "\n\n")
        
        total_count = sum(overall_counts.values())
        for key in order:
            percentage = (overall_counts[key] / total_count * 100) if total_count > 0 else 0
            f_out.write(f"{key}: {overall_counts[key]} ({percentage:.1f}%)\n")
        
        f_out.write(f"\nTotal: {total_count} structures\n")
    
    # Save per-kinase results
    per_kinase_file = output_file.replace('.txt', '_per_kinase.txt')
    with open(per_kinase_file, 'w', encoding='utf-8') as f_out:
        f_out.write(f"Per-Kinase Conformational State Distribution\n")
        f_out.write(f"="*60 + "\n\n")
        
        for kinase_name in sorted(kinase_counts.keys()):
            f_out.write(f"KINASE: {kinase_name.upper()}\n")
            f_out.write(f"-" * 40 + "\n")
            
            kinase_total = sum(kinase_counts[kinase_name].values())
            for key in order:
                count = kinase_counts[kinase_name][key]
                percentage = (count / kinase_total * 100) if kinase_total > 0 else 0
                f_out.write(f"{key}: {count} ({percentage:.1f}%)\n")
            
            f_out.write(f"Total: {kinase_total} structures\n\n")
    
    # NEW: Save structure mapping file - THIS IS THE KEY OUTPUT FOR DOCKING ANALYSIS
    mapping_file = output_file.replace('.txt', '_structure_mapping.csv')
    with open(mapping_file, 'w', encoding='utf-8') as f_map:
        f_map.write("kinase,structure_id,filename,conformational_state\n")
        
        # Sort by kinase name and structure_id for consistency
        structure_mapping.sort(key=lambda x: (x[0], x[1]))
        
        for kinase_name, structure_id, filename, state in structure_mapping:
            # Convert kinase name to uppercase for consistency with docking data
            f_map.write(f"{kinase_name.upper()},{structure_id},{filename},{state}\n")
    
    # Validation: Check for missing structure_ids per kinase
    validation_file = output_file.replace('.txt', '_validation_report.txt')
    with open(validation_file, 'w', encoding='utf-8') as f_val:
        f_val.write("Structure ID Validation Report\n")
        f_val.write("="*40 + "\n\n")
        
        kinase_structure_ids = defaultdict(set)
        for kinase_name, structure_id, filename, state in structure_mapping:
            kinase_structure_ids[kinase_name].add(structure_id)
        
        for kinase_name in sorted(kinase_structure_ids.keys()):
            structure_ids = kinase_structure_ids[kinase_name]
            expected_ids = set(range(0, 49))  # Assuming structures 1-49
            missing_ids = expected_ids - structure_ids
            extra_ids = structure_ids - expected_ids
            
            f_val.write(f"KINASE: {kinase_name}\n")
            f_val.write(f"Total structures: {len(structure_ids)}\n")
            f_val.write(f"Expected range: 0-49 ({len(expected_ids)} structures)\n")
            
            if missing_ids:
                f_val.write(f"MISSING structure_ids: {sorted(missing_ids)}\n")
            else:
                f_val.write("All expected structure_ids present\n")
            
            if extra_ids:
                f_val.write(f"EXTRA structure_ids: {sorted(extra_ids)}\n")
            
            f_val.write(f"\n")
    
    # Save files by conformational state - Overall
    state_files_file = output_file.replace('.txt', '_files_by_state.txt')
    with open(state_files_file, 'w', encoding='utf-8') as f_out:
        f_out.write(f"PDB Files by Conformational State (Overall)\n")
        f_out.write(f"="*60 + "\n\n")
        
        for state in order:
            if overall_counts[state] > 0:
                f_out.write(f"STATE: {state} ({overall_counts[state]} files)\n")
                f_out.write(f"-" * 50 + "\n")
                
                # Sort files by kinase name for better organization
                sorted_files = sorted(state_files[state])
                for file_path in sorted_files:
                    f_out.write(f"  {file_path}\n")
                f_out.write(f"\n")
    
    # Save files by conformational state - Per kinase (reorganized format)
    kinase_state_files_file = output_file.replace('.txt', '_files_by_state_per_kinase.txt')
    with open(kinase_state_files_file, 'w', encoding='utf-8') as f_out:
        f_out.write(f"PDB Files by Conformational State (Organized by Kinase)\n")
        f_out.write(f"="*60 + "\n\n")
        
        for kinase_name in sorted(kinase_counts.keys()):
            f_out.write(f"KINASE: {kinase_name.upper()}\n")
            f_out.write(f"=" * 40 + "\n")
            
            # Only show states that have files for this kinase
            kinase_has_files = False
            for state in order:
                if kinase_counts[kinase_name][state] > 0:
                    kinase_has_files = True
                    f_out.write(f"\n  {state}: ({kinase_counts[kinase_name][state]} files)\n")
                    f_out.write(f"  " + "-" * 45 + "\n")
                    
                    # Sort files for better organization
                    sorted_files = sorted(kinase_state_files[kinase_name][state])
                    for filename in sorted_files:
                        f_out.write(f"    {filename}\n")
            
            if not kinase_has_files:
                f_out.write(f"  No files found\n")
            f_out.write(f"\n")
    
    # Also create a files-by-state overview (for quick lookup of which kinases have each state)
    state_overview_file = output_file.replace('.txt', '_state_overview.txt')
    with open(state_overview_file, 'w', encoding='utf-8') as f_out:
        f_out.write(f"Conformational State Overview (Which Kinases Have Each State)\n")
        f_out.write(f"="*60 + "\n\n")
        
        for state in order:
            if overall_counts[state] > 0:
                f_out.write(f"STATE: {state} (Total: {overall_counts[state]} files)\n")
                f_out.write(f"-" * 50 + "\n")
                
                # Group by kinase and count files per kinase for this state
                kinase_file_counts = defaultdict(int)
                for file_path in state_files[state]:
                    kinase = file_path.split('/')[0]
                    kinase_file_counts[kinase] += 1
                
                # Sort kinases by number of files in this state (descending)
                sorted_kinases = sorted(kinase_file_counts.items(), key=lambda x: x[1], reverse=True)
                
                for kinase, count in sorted_kinases:
                    percentage = (count / overall_counts[state] * 100) if overall_counts[state] > 0 else 0
                    f_out.write(f"  {kinase.upper()}: {count} files ({percentage:.1f}% of this state)\n")
                f_out.write(f"\n")
    
    print(f"Overall conformational state analysis saved to {output_file}")
    print(f"Per-kinase conformational state analysis saved to {per_kinase_file}")
    print(f"Files by conformational state (overall) saved to {state_files_file}")
    print(f"Files by kinase and state saved to {kinase_state_files_file}")
    print(f"State overview (which kinases have each state) saved to {state_overview_file}")
    print(f"*** STRUCTURE MAPPING FOR DOCKING ANALYSIS saved to {mapping_file} ***")
    print(f"Validation report saved to {validation_file}")
    
    return overall_counts, kinase_counts, state_files, kinase_state_files, structure_mapping

def main():
    """Main pipeline execution."""
    print("=== Kinase Conformational State Analysis Pipeline ===")
    print(f"Using {NUM_CORES} CPU cores for parallel processing")
    print(f"Total CPU cores available: {mp.cpu_count()}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Step 1: Find all PDB files and add chain IDs
    print("\nStep 1: Adding chain IDs to PDB files...")
    
    all_pdb_files = []
    chainA_files = []
    
    for kinase_dir in sorted(os.listdir(RECEPTOR_BASE_DIR)):
        kinase_path = os.path.join(RECEPTOR_BASE_DIR, kinase_dir)
        if os.path.isdir(kinase_path):
            print(f"Processing kinase: {kinase_dir}")
            
            # Find all PDB files in this kinase directory
            pdb_files = [f for f in os.listdir(kinase_path) if f.endswith('.pdb') and not f.endswith('_chainA.pdb')]
            
            for pdb_file in sorted(pdb_files):
                pdb_path = os.path.join(kinase_path, pdb_file)
                all_pdb_files.append(pdb_path)
                
                # Extract structure_id from filename
                structure_id = extract_structure_id_from_filename(pdb_file.replace('.pdb', '_chainA.pdb'))
                if structure_id is None:
                    print(f"WARNING: Could not extract structure_id from {pdb_file}")
                    continue
                
                # Add chain ID
                chainA_file = add_chain_id_to_pdb(pdb_path, CHAIN_ID)
                if chainA_file:
                    chainA_files.append((chainA_file, kinase_dir, structure_id))  # Store kinase name and structure_id with file
    
    print(f"Found {len(chainA_files)} files with chain IDs")
    
    # Step 2: Run kincore on all chainA files (parallel processing)
    print(f"\nStep 2: Running kincore analysis on {len(chainA_files)} files using {NUM_CORES} CPU cores...")
    print(f"Estimated time reduction: ~{mp.cpu_count()-1}x faster than single-core processing")
    
    kincore_output_path = os.path.join(OUTPUT_DIR, KINCORE_OUTPUT_FILE)
    successful_analyses = 0
    
    # Use multiprocessing to parallelize kincore analysis
    print("Starting parallel kincore analysis...")
    with mp.Pool(NUM_CORES) as pool:
        # Submit all jobs to the pool
        results = pool.map(kincore_worker, chainA_files)
    
    print("Parallel processing complete. Writing results...")
    
    # Write results to output file and track success/failures
    with open(kincore_output_path, 'w') as output_file:
        for i, result in enumerate(results, 1):
            if i % 50 == 0 or i == len(results):  # Progress every 50 files
                print(f"Writing results: {i}/{len(chainA_files)}")
            
            if result['success']:
                output_file.write(result['output'])
                successful_analyses += 1
            else:
                print(f"  ERROR - {result['kinase']}/structure_{result['structure_id']}: {result['error']}")
    
    print(f"Successfully analyzed {successful_analyses}/{len(chainA_files)} files")
    
    # Step 3: Process conformational states
    print(f"\nStep 3: Processing conformational states...")
    
    final_counts_path = os.path.join(OUTPUT_DIR, FINAL_COUNTS_FILE)
    overall_counts, kinase_counts, state_files, kinase_state_files, structure_mapping = process_conformational_states(kincore_output_path, final_counts_path)
    
    # Print summary to console
    print(f"\n=== FINAL RESULTS ===")
    print(f"Total unique structures analyzed: {sum(overall_counts.values())}")
    print(f"Overall results saved to: {final_counts_path}")
    print(f"Per-kinase results saved to: {final_counts_path.replace('.txt', '_per_kinase.txt')}")
    print(f"Files by state (overall) saved to: {final_counts_path.replace('.txt', '_files_by_state.txt')}")
    print(f"Files by kinase and state saved to: {final_counts_path.replace('.txt', '_files_by_state_per_kinase.txt')}")
    print(f"State overview saved to: {final_counts_path.replace('.txt', '_state_overview.txt')}")
    print(f"Raw kincore output saved to: {kincore_output_path}")
    print(f"\n*** KEY OUTPUT FOR DOCKING ANALYSIS: {final_counts_path.replace('.txt', '_structure_mapping.csv')} ***")
    print(f"Validation report: {final_counts_path.replace('.txt', '_validation_report.txt')}")
    
    # Print top conformational states
    sorted_counts = sorted(overall_counts.items(), key=lambda x: x[1], reverse=True)
    print(f"\nTop overall conformational states:")
    for state, count in sorted_counts[:5]:
        if count > 0:
            percentage = count / sum(overall_counts.values()) * 100
            print(f"  {state}: {count} ({percentage:.1f}%) - {len(state_files[state])} files")
    
    # Print summary per kinase
    print(f"\nPer-kinase summary:")
    for kinase in sorted(kinase_counts.keys()):
        total = sum(kinase_counts[kinase].values())
        top_state = max(kinase_counts[kinase].items(), key=lambda x: x[1])
        print(f"  {kinase}: {total} structures, top state: {top_state[0]} ({top_state[1]} structures)")

if __name__ == "__main__":
    # Ensure proper multiprocessing behavior
    main()
