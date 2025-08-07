#!/usr/bin/env python3
"""
Full DLG Processor - Process all kinases
"""

import os
import re
import sys
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool, cpu_count
import json

# Force unbuffered output
import sys
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

class AdaDLGProcessor:
    """DLG processor for Ada cluster docking output"""
    
    def __init__(self, n_processes=None):
        self.n_processes = n_processes or max(1, cpu_count() - 1)
        print(f"Initialized with {self.n_processes} processes", flush=True)
        
    def parse_dlg_file(self, dlg_path):
        """Parse a single DLG file"""
        try:
            with open(dlg_path, 'r') as f:
                content = f.read()
            
            # Extract ligand name from filename pattern
            filename = os.path.basename(dlg_path)
            chembl_match = re.search(r'REMARKName=([^.]+)', filename)
            if chembl_match:
                ligand_name = chembl_match.group(1)
            else:
                ligand_match = re.search(r'REMARK\s+Name\s*=\s*(\S+)', content)
                if ligand_match:
                    ligand_name = ligand_match.group(1)
                else:
                    ligand_name = Path(dlg_path).stem
            
            # Find all binding energies
            energies = re.findall(r'Estimated Free Energy of Binding\s*=\s*([-+]?\d*\.?\d+)\s*kcal/mol', content)
            
            if energies:
                energies = [float(e) for e in energies]
                return {
                    'ligand_name': ligand_name,
                    'best_score': min(energies),
                    'mean_score': np.mean(energies),
                    'std_score': np.std(energies) if len(energies) > 1 else 0.0,
                    'n_poses': len(energies)
                }
            else:
                return None
                
        except Exception as e:
            return None
    
    def determine_compound_type(self, ligand_name, file_path):
        """Determine if compound is active or decoy"""
        basename = os.path.basename(file_path)
        
        if basename.startswith('actives_'):
            return 'active'
        elif basename.startswith('decoys_') or basename.startswith('decoy_'):
            return 'decoy'
        else:
            return 'decoy'
    
    def process_receptor_directory(self, args):
        """Process all DLG files in a receptor directory"""
        receptor_path, kinase_name, receptor_id = args
        results = []
        
        # Find all DLG files
        dlg_files = glob.glob(os.path.join(receptor_path, "*.dlg"))
        
        if not dlg_files:
            return results
        
        # Process each file
        for dlg_file in dlg_files:
            result = self.parse_dlg_file(dlg_file)
            if result:
                result['kinase'] = kinase_name.upper()
                result['structure_id'] = receptor_id
                result['file_path'] = dlg_file
                result['compound_type'] = self.determine_compound_type(
                    result['ligand_name'], dlg_file
                )
                results.append(result)
        
        # Print progress for this receptor
        print(f"  Processed {kinase_name}/receptor_{receptor_id:03d}: {len(results)} compounds", flush=True)
        
        return results
    
    def process_all_kinases(self, base_directory):
        """Process all kinases in the directory"""
        print(f"\nFULL MODE: Processing all kinases", flush=True)
        print(f"Base directory: {base_directory}", flush=True)
        start_time = datetime.now()
        
        # Get all kinase directories
        kinase_dirs = sorted([d for d in os.listdir(base_directory) 
                            if os.path.isdir(os.path.join(base_directory, d))])
        
        print(f"Found {len(kinase_dirs)} kinases: {', '.join(kinase_dirs)}", flush=True)
        
        # Collect all tasks
        all_tasks = []
        kinase_receptor_counts = {}
        
        for kinase_name in kinase_dirs:
            kinase_path = os.path.join(base_directory, kinase_name)
            
            # List receptor directories
            receptor_dirs = sorted([d for d in os.listdir(kinase_path) 
                                  if os.path.isdir(os.path.join(kinase_path, d)) 
                                  and d.startswith('receptor_')])
            
            kinase_receptor_counts[kinase_name] = len(receptor_dirs)
            
            # Add tasks for this kinase
            for receptor_dir in receptor_dirs:
                receptor_path = os.path.join(kinase_path, receptor_dir)
                receptor_match = re.search(r'receptor_(\d+)', receptor_dir)
                receptor_id = int(receptor_match.group(1)) if receptor_match else 1
                all_tasks.append((receptor_path, kinase_name, receptor_id))
        
        # Print summary
        total_receptors = sum(kinase_receptor_counts.values())
        print(f"\nTotal receptor directories to process: {total_receptors}", flush=True)
        for kinase, count in kinase_receptor_counts.items():
            print(f"  {kinase}: {count} receptors", flush=True)
        
        # Process all tasks
        print(f"\nProcessing {len(all_tasks)} receptor directories...", flush=True)
        all_results = []
        
        with Pool(self.n_processes) as pool:
            for i, results in enumerate(pool.imap_unordered(self.process_receptor_directory, all_tasks)):
                all_results.extend(results)
                if (i + 1) % 50 == 0:
                    print(f"Progress: {i + 1}/{len(all_tasks)} receptors completed", flush=True)
        
        print(f"\nTotal compounds processed: {len(all_results)}", flush=True)
        
        if not all_results:
            print("ERROR: No results found!", flush=True)
            return None
        
        # Convert to DataFrame
        df = pd.DataFrame(all_results)
        
        # Create clean dataset
        clean_df = df[['kinase', 'structure_id', 'ligand_name', 'compound_type', 
                      'best_score', 'mean_score', 'std_score', 'n_poses']].copy()
        clean_df.rename(columns={
            'ligand_name': 'compound_id',
            'best_score': 'docking_score'
        }, inplace=True)
        
        # Save outputs
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # Save clean CSV
        output_csv = f"docking_results_clean_{timestamp}.csv"
        clean_df.to_csv(output_csv, index=False)
        print(f"\nSaved clean results to: {output_csv}", flush=True)
        
        # Generate and save summary
        summary = self.generate_summary(df, base_directory)
        summary_file = f"docking_results_summary_{timestamp}.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"Saved summary to: {summary_file}", flush=True)
        
        # Print final summary
        elapsed = datetime.now() - start_time
        print(f"\n{'='*60}", flush=True)
        print(f"Processing complete in {elapsed}", flush=True)
        print(f"Total kinases: {df['kinase'].nunique()}", flush=True)
        print(f"Total receptors: {len(df.groupby(['kinase', 'structure_id']))}", flush=True)
        print(f"Total unique compounds: {df['ligand_name'].nunique()}", flush=True)
        print(f"Total docking results: {len(df)}", flush=True)
        print(f"{'='*60}", flush=True)
        
        # Show per-kinase summary
        print("\nPer-kinase summary:", flush=True)
        for kinase in sorted(df['kinase'].unique()):
            kinase_df = df[df['kinase'] == kinase]
            type_counts = kinase_df['compound_type'].value_counts()
            print(f"  {kinase}: {type_counts.get('active', 0)} actives, {type_counts.get('decoy', 0)} decoys", flush=True)
        
        return clean_df
    
    def generate_summary(self, df, base_directory):
        """Generate summary statistics"""
        summary = {
            'processing_date': datetime.now().isoformat(),
            'base_directory': base_directory,
            'total_kinases': df['kinase'].nunique(),
            'total_structures': len(df.groupby(['kinase', 'structure_id'])),
            'total_compounds': df['ligand_name'].nunique(),
            'total_docking_runs': len(df),
            'kinase_summary': {}
        }
        
        for kinase in sorted(df['kinase'].unique()):
            kinase_df = df[df['kinase'] == kinase]
            
            # Calculate statistics
            active_scores = kinase_df[kinase_df['compound_type'] == 'active']['best_score']
            decoy_scores = kinase_df[kinase_df['compound_type'] == 'decoy']['best_score']
            
            summary['kinase_summary'][kinase] = {
                'n_structures': kinase_df['structure_id'].nunique(),
                'n_compounds': kinase_df['ligand_name'].nunique(),
                'n_actives': len(kinase_df[kinase_df['compound_type'] == 'active']),
                'n_decoys': len(kinase_df[kinase_df['compound_type'] == 'decoy']),
                'avg_active_score': float(active_scores.mean()) if len(active_scores) > 0 else None,
                'avg_decoy_score': float(decoy_scores.mean()) if len(decoy_scores) > 0 else None,
                'score_separation': float(decoy_scores.mean() - active_scores.mean()) 
                                   if len(active_scores) > 0 and len(decoy_scores) > 0 else None
            }
        
        return summary
    
    def process_test_kinase(self, base_directory, kinase_name):
        """Process just one kinase for testing"""
        print(f"\nTEST MODE: Processing only {kinase_name}", flush=True)
        kinase_path = os.path.join(base_directory, kinase_name)
        
        if not os.path.exists(kinase_path):
            print(f"Error: Kinase directory not found: {kinase_path}", flush=True)
            return None
        
        # List receptor directories
        receptor_dirs = sorted([d for d in os.listdir(kinase_path) 
                               if os.path.isdir(os.path.join(kinase_path, d)) 
                               and d.startswith('receptor_')])
        
        print(f"Found {len(receptor_dirs)} receptor directories", flush=True)
        
        # Prepare tasks
        tasks = []
        for receptor_dir in receptor_dirs:
            receptor_path = os.path.join(kinase_path, receptor_dir)
            receptor_match = re.search(r'receptor_(\d+)', receptor_dir)
            receptor_id = int(receptor_match.group(1)) if receptor_match else 1
            tasks.append((receptor_path, kinase_name, receptor_id))
        
        # Check first receptor
        if tasks:
            first_receptor = tasks[0][0]
            dlg_count = len(glob.glob(os.path.join(first_receptor, "*.dlg")))
            print(f"Sample: {first_receptor} contains {dlg_count} DLG files", flush=True)
        
        # Process with progress
        print(f"\nProcessing {len(tasks)} receptors...", flush=True)
        all_results = []
        
        with Pool(self.n_processes) as pool:
            for results in pool.imap(self.process_receptor_directory, tasks):
                all_results.extend(results)
        
        print(f"\nTotal compounds processed: {len(all_results)}", flush=True)
        
        if not all_results:
            print("ERROR: No results found!", flush=True)
            return None
        
        # Convert to DataFrame
        df = pd.DataFrame(all_results)
        
        # Create clean dataset
        clean_df = df[['kinase', 'structure_id', 'ligand_name', 'compound_type', 
                      'best_score', 'mean_score', 'std_score', 'n_poses']].copy()
        clean_df.rename(columns={
            'ligand_name': 'compound_id',
            'best_score': 'docking_score'
        }, inplace=True)
        
        # Save output
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_csv = f"test_{kinase_name}_results_{timestamp}.csv"
        clean_df.to_csv(output_csv, index=False)
        print(f"\nSaved results to: {output_csv}", flush=True)
        
        # Show summary
        type_counts = df['compound_type'].value_counts()
        print(f"\nSummary for {kinase_name}:", flush=True)
        print(f"  Total compounds: {len(df)}", flush=True)
        print(f"  Unique compounds: {df['ligand_name'].nunique()}", flush=True)
        print(f"  Actives: {type_counts.get('active', 0)}", flush=True)
        print(f"  Decoys: {type_counts.get('decoy', 0)}", flush=True)
        print(f"  Receptors: {df['structure_id'].nunique()}", flush=True)
        
        return clean_df


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Process DLG files on Ada cluster')
    parser.add_argument('--base-dir', 
                       default='/store/jaeohshin/work/dock/virtual_screening/docking_output',
                       help='Base directory containing kinase folders')
    parser.add_argument('--processes', type=int, default=20,
                       help='Number of parallel processes')
    parser.add_argument('--test-kinase', 
                       help='Process only one kinase for testing')
    parser.add_argument('--all', action='store_true',
                       help='Process all kinases')
    
    args = parser.parse_args()
    
    print(f"Starting processing at {datetime.now()}", flush=True)
    print(f"Base directory: {args.base_dir}", flush=True)
    print(f"Using {args.processes} processes", flush=True)
    
    processor = AdaDLGProcessor(n_processes=args.processes)
    
    if args.all:
        # Process all kinases
        processor.process_all_kinases(args.base_dir)
    elif args.test_kinase:
        # Process test kinase
        processor.process_test_kinase(args.base_dir, args.test_kinase)
    else:
        print("Please specify either --test-kinase KINASE_NAME or --all")
        parser.print_help()
    
    print(f"\nCompleted at {datetime.now()}", flush=True)


if __name__ == '__main__':
    main()