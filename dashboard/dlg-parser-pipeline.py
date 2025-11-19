import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import glob
from tqdm import tqdm

class DLGParser:
    """Parser for AutoDock DLG (Docking Log) files"""
    
    def __init__(self):
        self.results = []
        
    def parse_dlg_file(self, dlg_path):
        """
        Parse a single DLG file and extract docking results
        
        Returns:
            dict: Parsed results including scores, poses, and metadata
        """
        results = {
            'file_path': dlg_path,
            'runs': []
        }
        
        with open(dlg_path, 'r') as f:
            content = f.read()
        
        # Extract ligand name (usually in the REMARK section)
        ligand_match = re.search(r'REMARK\s+Name\s*=\s*(\S+)', content)
        if ligand_match:
            results['ligand_name'] = ligand_match.group(1)
        else:
            # Try to get from filename
            results['ligand_name'] = Path(dlg_path).stem
        
        # Find all docking runs
        run_pattern = r'DOCKED: USER\s+Run\s*=\s*(\d+).*?DOCKED: USER\s+Estimated Free Energy of Binding\s*=\s*([-\d.]+)\s*kcal/mol.*?DOCKED: USER\s+\(1\)\s+Final Intermolecular Energy\s*=\s*([-\d.]+)\s*kcal/mol.*?DOCKED: USER\s+\(2\)\s+Final Total Internal Energy\s*=\s*([-\d.]+)\s*kcal/mol.*?DOCKED: USER\s+\(3\)\s+Torsional Free Energy\s*=\s*([-\d.]+)\s*kcal/mol'
        
        runs = re.finditer(run_pattern, content, re.DOTALL)
        
        for run in runs:
            run_data = {
                'run_number': int(run.group(1)),
                'binding_energy': float(run.group(2)),
                'intermolecular_energy': float(run.group(3)),
                'internal_energy': float(run.group(4)),
                'torsional_energy': float(run.group(5))
            }
            results['runs'].append(run_data)
        
        # Find clustering information if available
        cluster_pattern = r'CLUSTERING HISTOGRAM.*?RANKING(.*?)(?:_____|$)'
        cluster_match = re.search(cluster_pattern, content, re.DOTALL)
        
        if cluster_match:
            ranking_text = cluster_match.group(1)
            rank_pattern = r'(\d+)\s+([-\d.]+)\s+(\d+)\s+(\d+)'
            rankings = re.findall(rank_pattern, ranking_text)
            
            results['clusters'] = []
            for rank in rankings:
                cluster_data = {
                    'rank': int(rank[0]),
                    'lowest_binding_energy': float(rank[1]),
                    'run_number': int(rank[2]),
                    'cluster_size': int(rank[3])
                }
                results['clusters'].append(cluster_data)
        
        return results
    
    def parse_directory(self, directory_path, pattern="*.dlg"):
        """
        Parse all DLG files in a directory
        
        Args:
            directory_path: Path to directory containing DLG files
            pattern: Glob pattern for DLG files
            
        Returns:
            pd.DataFrame: Parsed results from all files
        """
        dlg_files = glob.glob(os.path.join(directory_path, pattern))
        
        all_results = []
        
        print(f"Found {len(dlg_files)} DLG files to parse...")
        
        for dlg_file in tqdm(dlg_files, desc="Parsing DLG files"):
            try:
                result = self.parse_dlg_file(dlg_file)
                
                # Extract best docking score (lowest energy)
                if result['runs']:
                    best_run = min(result['runs'], key=lambda x: x['binding_energy'])
                    
                    # Flatten the data for DataFrame
                    flat_result = {
                        'file_path': dlg_file,
                        'ligand_name': result['ligand_name'],
                        'best_binding_energy': best_run['binding_energy'],
                        'best_run_number': best_run['run_number'],
                        'num_runs': len(result['runs']),
                        'all_binding_energies': [r['binding_energy'] for r in result['runs']]
                    }
                    
                    # Add cluster information if available
                    if 'clusters' in result and result['clusters']:
                        best_cluster = result['clusters'][0]  # Already ranked
                        flat_result['best_cluster_energy'] = best_cluster['lowest_binding_energy']
                        flat_result['best_cluster_size'] = best_cluster['cluster_size']
                    
                    all_results.append(flat_result)
                    
            except Exception as e:
                print(f"Error parsing {dlg_file}: {str(e)}")
                continue
        
        return pd.DataFrame(all_results)


def organize_docking_data(base_directory, kinase_structure_mapping=None):
    """
    Organize docking data from DLG files into a structured format
    
    Expected directory structure:
    base_directory/
    ├── KINASE1/
    │   ├── structure_1/
    │   │   ├── active_compound1.dlg
    │   │   ├── active_compound2.dlg
    │   │   ├── decoy_compound1.dlg
    │   │   └── ...
    │   ├── structure_2/
    │   └── ...
    ├── KINASE2/
    └── ...
    
    Args:
        base_directory: Root directory containing all docking results
        kinase_structure_mapping: Optional dict mapping directory names to kinase names
        
    Returns:
        pd.DataFrame: Organized docking results
    """
    parser = DLGParser()
    all_data = []
    
    # Walk through directory structure
    for kinase_dir in os.listdir(base_directory):
        kinase_path = os.path.join(base_directory, kinase_dir)
        
        if not os.path.isdir(kinase_path):
            continue
            
        kinase_name = kinase_structure_mapping.get(kinase_dir, kinase_dir) if kinase_structure_mapping else kinase_dir
        
        print(f"\nProcessing kinase: {kinase_name}")
        
        for structure_dir in os.listdir(kinase_path):
            structure_path = os.path.join(kinase_path, structure_dir)
            
            if not os.path.isdir(structure_path):
                continue
            
            # Extract structure ID from directory name
            structure_id_match = re.search(r'structure_?(\d+)', structure_dir, re.IGNORECASE)
            structure_id = int(structure_id_match.group(1)) if structure_id_match else structure_dir
            
            print(f"  Processing structure: {structure_dir}")
            
            # Parse all DLG files in this structure directory
            df_structure = parser.parse_directory(structure_path)
            
            if not df_structure.empty:
                # Add kinase and structure information
                df_structure['kinase'] = kinase_name
                df_structure['structure_id'] = structure_id
                
                # Determine compound type based on filename or directory structure
                # Adjust this logic based on your naming convention
                df_structure['compound_type'] = df_structure['ligand_name'].apply(
                    lambda x: 'active' if 'active' in x.lower() or 'actives' in x.lower() 
                    else 'decoy'
                )
                
                # You might need to adjust this based on your compound naming
                df_structure['compound_id'] = df_structure['ligand_name']
                
                all_data.append(df_structure)
    
    # Combine all data
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        
        # Rename columns to match dashboard expectations
        final_df = final_df.rename(columns={
            'best_binding_energy': 'docking_score'
        })
        
        return final_df
    else:
        return pd.DataFrame()


def prepare_dude_data(docking_results_df, dude_actives_file=None, dude_decoys_file=None):
    """
    If you have separate files listing DUD-E actives and decoys,
    use this function to properly label compound types
    
    Args:
        docking_results_df: DataFrame from organize_docking_data()
        dude_actives_file: Path to file listing active compounds
        dude_decoys_file: Path to file listing decoy compounds
        
    Returns:
        pd.DataFrame: Updated DataFrame with correct compound types
    """
    if dude_actives_file:
        with open(dude_actives_file, 'r') as f:
            actives = set(line.strip() for line in f)
        
        # Update compound types based on active list
        docking_results_df['compound_type'] = docking_results_df['compound_id'].apply(
            lambda x: 'active' if any(active in x for active in actives) else 'decoy'
        )
    
    return docking_results_df


def generate_summary_statistics(df):
    """
    Generate summary statistics for the docking results
    
    Args:
        df: DataFrame with parsed docking results
        
    Returns:
        dict: Summary statistics
    """
    summary = {
        'total_kinases': df['kinase'].nunique(),
        'total_structures': len(df.groupby(['kinase', 'structure_id'])),
        'total_compounds': df['compound_id'].nunique(),
        'total_docking_runs': len(df),
        'kinase_stats': []
    }
    
    # Per-kinase statistics
    for kinase in df['kinase'].unique():
        kinase_df = df[df['kinase'] == kinase]
        
        kinase_stat = {
            'kinase': kinase,
            'n_structures': kinase_df['structure_id'].nunique(),
            'n_actives': len(kinase_df[kinase_df['compound_type'] == 'active']['compound_id'].unique()),
            'n_decoys': len(kinase_df[kinase_df['compound_type'] == 'decoy']['compound_id'].unique()),
            'avg_active_score': kinase_df[kinase_df['compound_type'] == 'active']['docking_score'].mean(),
            'avg_decoy_score': kinase_df[kinase_df['compound_type'] == 'decoy']['docking_score'].mean(),
            'score_separation': (
                kinase_df[kinase_df['compound_type'] == 'decoy']['docking_score'].mean() -
                kinase_df[kinase_df['compound_type'] == 'active']['docking_score'].mean()
            )
        }
        
        summary['kinase_stats'].append(kinase_stat)
    
    return summary


# Example usage:
if __name__ == "__main__":
    # Example 1: Parse a single DLG file
    parser = DLGParser()
    single_result = parser.parse_dlg_file("path/to/your/file.dlg")
    print("Single file results:", single_result)
    
    # Example 2: Parse all DLG files in a directory
    df_results = parser.parse_directory("path/to/dlg/directory/")
    print(f"Parsed {len(df_results)} files")
    print(df_results.head())
    
    # Example 3: Organize complete docking dataset
    # Assuming directory structure: base_dir/KINASE_NAME/structure_X/*.dlg
    base_dir = "/path/to/your/docking/results"
    
    # If your directory names don't match kinase names exactly:
    kinase_mapping = {
        'abl1': 'ABL1',
        'akt1': 'AKT1',
        'akt2': 'AKT2',
        # ... add more mappings as needed
    }
    
    # Process all data
    full_df = organize_docking_data(base_dir, kinase_mapping)
    
    # Save to CSV for use with the dashboard
    full_df.to_csv('docking_results_processed.csv', index=False)
    
    # Generate summary
    summary = generate_summary_statistics(full_df)
    print("\nSummary Statistics:")
    print(f"Total kinases: {summary['total_kinases']}")
    print(f"Total structures: {summary['total_structures']}")
    print(f"Total compounds: {summary['total_compounds']}")
    print(f"Total docking runs: {summary['total_docking_runs']}")
    
    # Show per-kinase statistics
    kinase_stats_df = pd.DataFrame(summary['kinase_stats'])
    print("\nPer-kinase statistics:")
    print(kinase_stats_df.to_string())