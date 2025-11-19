"""
Integrated DLG Parser and Interactive Dashboard for Kinase Virtual Screening
This combines DLG file parsing with the interactive dashboard
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import glob
from tqdm import tqdm

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from sklearn.metrics import roc_curve, auc
import plotly.figure_factory as ff

# ============================================================================
# PART 1: DLG PARSER
# ============================================================================

class DLGParser:
    """Parser for AutoDock DLG (Docking Log) files"""
    
    def __init__(self):
        self.results = []
        
    def parse_dlg_file(self, dlg_path):
        """Parse a single DLG file and extract docking results"""
        results = {
            'file_path': dlg_path,
            'runs': []
        }
        
        with open(dlg_path, 'r') as f:
            content = f.read()
        
        # Extract ligand name
        ligand_match = re.search(r'REMARK\s+Name\s*=\s*(\S+)', content)
        if ligand_match:
            results['ligand_name'] = ligand_match.group(1)
        else:
            results['ligand_name'] = Path(dlg_path).stem
        
        # Find all docking runs - simplified pattern
        run_pattern = r'Run\s*=\s*(\d+).*?Estimated Free Energy of Binding\s*=\s*([-\d.]+)\s*kcal/mol'
        runs = re.finditer(run_pattern, content, re.DOTALL)
        
        for run in runs:
            run_data = {
                'run_number': int(run.group(1)),
                'binding_energy': float(run.group(2))
            }
            results['runs'].append(run_data)
        
        return results
    
    def parse_directory(self, directory_path, pattern="*.dlg"):
        """Parse all DLG files in a directory"""
        dlg_files = glob.glob(os.path.join(directory_path, pattern))
        all_results = []
        
        print(f"Found {len(dlg_files)} DLG files in {directory_path}")
        
        for dlg_file in dlg_files:
            try:
                result = self.parse_dlg_file(dlg_file)
                
                if result['runs']:
                    best_run = min(result['runs'], key=lambda x: x['binding_energy'])
                    
                    flat_result = {
                        'file_path': dlg_file,
                        'ligand_name': result['ligand_name'],
                        'best_binding_energy': best_run['binding_energy'],
                        'best_run_number': best_run['run_number'],
                        'num_runs': len(result['runs']),
                        'all_binding_energies': [r['binding_energy'] for r in result['runs']]
                    }
                    
                    all_results.append(flat_result)
                    
            except Exception as e:
                print(f"Error parsing {dlg_file}: {str(e)}")
                continue
        
        return pd.DataFrame(all_results)


# ============================================================================
# PART 2: DATA PROCESSING FUNCTIONS
# ============================================================================

def process_docking_data(base_directory):
    """
    Process docking data from directory structure.
    
    Expected structure:
    base_directory/
    ├── KINASE1/
    │   ├── structure_1/
    │   │   ├── compound1.dlg
    │   │   └── compound2.dlg
    │   └── structure_2/
    └── KINASE2/
    """
    parser = DLGParser()
    all_data = []
    
    print(f"Processing docking data from: {base_directory}")
    
    # Get all kinase directories
    kinase_dirs = [d for d in os.listdir(base_directory) 
                   if os.path.isdir(os.path.join(base_directory, d))]
    
    for kinase_dir in tqdm(kinase_dirs, desc="Processing kinases"):
        kinase_path = os.path.join(base_directory, kinase_dir)
        kinase_name = kinase_dir.upper()  # Standardize kinase names
        
        # Get all structure directories
        structure_dirs = [d for d in os.listdir(kinase_path) 
                         if os.path.isdir(os.path.join(kinase_path, d))]
        
        for structure_dir in structure_dirs:
            structure_path = os.path.join(kinase_path, structure_dir)
            
            # Extract structure ID
            structure_match = re.search(r'(\d+)', structure_dir)
            structure_id = int(structure_match.group(1)) if structure_match else 1
            
            # Parse all DLG files
            df_structure = parser.parse_directory(structure_path)
            
            if not df_structure.empty:
                df_structure['kinase'] = kinase_name
                df_structure['structure_id'] = structure_id
                
                # Determine compound type from filename or path
                # Adjust this based on your naming convention
                df_structure['compound_type'] = df_structure['ligand_name'].apply(
                    lambda x: 'active' if 'active' in x.lower() else 'decoy'
                )
                
                df_structure['compound_id'] = df_structure['ligand_name']
                df_structure['docking_score'] = df_structure['best_binding_energy']
                
                all_data.append(df_structure)
    
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()


def load_or_process_data(base_directory=None, csv_file=None):
    """
    Load data from CSV if available, otherwise process from DLG files
    """
    if csv_file and os.path.exists(csv_file):
        print(f"Loading processed data from {csv_file}")
        return pd.read_csv(csv_file)
    elif base_directory:
        print(f"Processing DLG files from {base_directory}")
        df = process_docking_data(base_directory)
        
        # Save for future use
        output_file = 'processed_docking_results.csv'
        df.to_csv(output_file, index=False)
        print(f"Saved processed data to {output_file}")
        
        return df
    else:
        print("No data source provided. Using sample data.")
        return generate_sample_data()


def generate_sample_data():
    """Generate sample data for demo purposes"""
    kinases = [f'KINASE_{i}' for i in range(1, 26)]
    data = []
    
    np.random.seed(42)
    for kinase in kinases:
        performance_offset = np.random.uniform(-2, 2)
        
        for struct_id in range(1, 51):
            # Actives
            for compound_id in range(50):
                score = np.random.normal(-8 + performance_offset, 1.5)
                data.append({
                    'kinase': kinase,
                    'structure_id': struct_id,
                    'compound_id': f'ACTIVE_{compound_id}',
                    'compound_type': 'active',
                    'docking_score': score
                })
            
            # Decoys
            for compound_id in range(500):
                score = np.random.normal(-5 + performance_offset, 2)
                data.append({
                    'kinase': kinase,
                    'structure_id': struct_id,
                    'compound_id': f'DECOY_{compound_id}',
                    'compound_type': 'decoy',
                    'docking_score': score
                })
    
    return pd.DataFrame(data)


# ============================================================================
# PART 3: DASHBOARD APPLICATION
# ============================================================================

# Initialize the Dash app
app = dash.Dash(__name__)

# Configuration - MODIFY THESE PATHS
BASE_DIRECTORY = "path/to/your/docking/results"  # Change this!
CSV_FILE = "processed_docking_results.csv"  # Will be created after first run

# Load data
print("Loading docking data...")
df = load_or_process_data(base_directory=BASE_DIRECTORY, csv_file=CSV_FILE)
print(f"Loaded {len(df)} docking results")
print(f"Kinases: {df['kinase'].nunique()}")
print(f"Unique compounds: {df['compound_id'].nunique()}")

# Calculate ensemble metrics
def calculate_ensemble_metrics(df):
    """Calculate AUC and enrichment factors using ensemble average"""
    results = []
    
    for kinase in df['kinase'].unique():
        kinase_data = df[df['kinase'] == kinase]
        
        # Calculate ensemble average score for each compound
        ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
        
        # Calculate ROC AUC
        actives = ensemble_scores[ensemble_scores['compound_type'] == 'active']['docking_score'].values
        decoys = ensemble_scores[ensemble_scores['compound_type'] == 'decoy']['docking_score'].values
        
        if len(actives) > 0 and len(decoys) > 0:
            y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
            y_scores = -np.concatenate([actives, decoys])  # Negative because more negative = better
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            roc_auc = auc(fpr, tpr)
            
            # Calculate enrichment factors
            all_compounds = ensemble_scores.sort_values('docking_score')
            n_actives = len(actives)
            n_total = len(all_compounds)
            
            # EF at 1%
            top_1_percent = max(1, int(0.01 * n_total))
            actives_in_top_1 = (all_compounds.head(top_1_percent)['compound_type'] == 'active').sum()
            ef_1 = (actives_in_top_1 / top_1_percent) / (n_actives / n_total) if n_actives > 0 else 0
            
            # EF at 5%
            top_5_percent = max(1, int(0.05 * n_total))
            actives_in_top_5 = (all_compounds.head(top_5_percent)['compound_type'] == 'active').sum()
            ef_5 = (actives_in_top_5 / top_5_percent) / (n_actives / n_total) if n_actives > 0 else 0
            
            results.append({
                'kinase': kinase,
                'AUC': roc_auc,
                'EF1%': ef_1,
                'EF5%': ef_5,
                'n_actives': n_actives,
                'n_decoys': len(decoys),
                'n_structures': kinase_data['structure_id'].nunique()
            })
    
    return pd.DataFrame(results)

metrics_df = calculate_ensemble_metrics(df)

# Dashboard layout
app.layout = html.Div([
    html.H1("Kinase Virtual Screening Dashboard", style={'text-align': 'center'}),
    
    # Data summary
    html.Div([
        html.P(f"Total kinases: {df['kinase'].nunique()} | "
               f"Total structures: {df.groupby(['kinase', 'structure_id']).ngroups} | "
               f"Total compounds: {df['compound_id'].nunique()} | "
               f"Total docking runs: {len(df)}",
               style={'text-align': 'center', 'font-size': '16px'})
    ]),
    
    # Top row: Overview metrics
    html.Div([
        html.Div([
            html.H3("Performance Overview"),
            dcc.Graph(id='heatmap-overview')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H3("Kinase Ranking by AUC"),
            dcc.Graph(id='ranking-bar')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Kinase selector
    html.Div([
        html.H3("Select Kinase for Detailed Analysis:"),
        dcc.Dropdown(
            id='kinase-dropdown',
            options=[{'label': k, 'value': k} for k in sorted(df['kinase'].unique())],
            value=sorted(df['kinase'].unique())[0]
        )
    ], style={'width': '48%', 'margin': 'auto', 'margin-top': '20px'}),
    
    # Selected kinase info
    html.Div(id='kinase-info', style={'text-align': 'center', 'margin': '10px'}),
    
    # Detailed analysis row
    html.Div([
        html.Div([
            html.H4("ROC Curve Analysis"),
            dcc.Graph(id='roc-curve')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H4("Score Distribution"),
            dcc.Graph(id='score-distribution')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Structure variability analysis
    html.Div([
        html.Div([
            html.H4("Structure Performance Variability"),
            dcc.Graph(id='structure-performance')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H4("Top Scoring Compounds"),
            html.Div(id='top-compounds-table')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Export button
    html.Div([
        html.Button("Export Current Kinase Data", id="export-button"),
        dcc.Download(id="download-dataframe-csv")
    ], style={'text-align': 'center', 'margin': '20px'})
])

# Callbacks
@app.callback(
    Output('heatmap-overview', 'figure'),
    Input('heatmap-overview', 'id')
)
def update_heatmap(_):
    # Prepare data for heatmap
    heatmap_data = metrics_df.set_index('kinase')[['AUC', 'EF1%', 'EF5%']].T
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdBu',
        text=np.round(heatmap_data.values, 2),
        texttemplate='%{text}',
        textfont={"size": 10},
        colorbar=dict(title="Value")
    ))
    
    fig.update_layout(
        title="Performance Metrics Across All Kinases",
        xaxis_title="Kinase",
        yaxis_title="Metric",
        height=400,
        xaxis={'tickangle': -45}
    )
    
    return fig

@app.callback(
    Output('ranking-bar', 'figure'),
    Input('ranking-bar', 'id')
)
def update_ranking(_):
    sorted_df = metrics_df.sort_values('AUC', ascending=True).tail(20)  # Top 20
    
    colors = ['red' if auc < 0.6 else 'yellow' if auc < 0.8 else 'green' 
              for auc in sorted_df['AUC']]
    
    fig = go.Figure(go.Bar(
        x=sorted_df['AUC'],
        y=sorted_df['kinase'],
        orientation='h',
        marker_color=colors,
        text=sorted_df['AUC'].round(3),
        textposition='auto'
    ))
    
    fig.add_vline(x=0.5, line_dash="dash", line_color="gray", 
                  annotation_text="Random")
    