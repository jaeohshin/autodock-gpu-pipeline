#!/usr/bin/env python3
"""
Updated Dashboard with EF10% and separate heatmaps for AUC and EF metrics
"""

import os
import sys
import pandas as pd
import numpy as np

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from sklearn.metrics import roc_curve, auc

# ============================================================================
# CONFIGURATION
# ============================================================================

CSV_FILE = "docking_results_clean_20250731_202701.csv"  # Update with your filename

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading docking results...")

if not os.path.exists(CSV_FILE):
    print(f"ERROR: CSV file '{CSV_FILE}' not found!")
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv')]
    print("\nAvailable CSV files:")
    for f in csv_files:
        print(f"  - {f}")
    sys.exit(1)

df = pd.read_csv(CSV_FILE)

print(f"Loaded {len(df)} docking results")
print(f"Kinases: {df['kinase'].nunique()}")

# ============================================================================
# CALCULATE ENSEMBLE METRICS WITH EF10%
# ============================================================================

def calculate_ensemble_metrics(df):
    """Calculate AUC and enrichment factors including EF10%"""
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
            y_scores = -np.concatenate([actives, decoys])
            
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
            
            # EF at 10%
            top_10_percent = max(1, int(0.10 * n_total))
            actives_in_top_10 = (all_compounds.head(top_10_percent)['compound_type'] == 'active').sum()
            ef_10 = (actives_in_top_10 / top_10_percent) / (n_actives / n_total) if n_actives > 0 else 0
            
            results.append({
                'kinase': kinase,
                'AUC': roc_auc,
                'EF1%': ef_1,
                'EF5%': ef_5,
                'EF10%': ef_10,
                'n_actives': n_actives,
                'n_decoys': len(decoys),
                'n_structures': kinase_data['structure_id'].nunique()
            })
    
    return pd.DataFrame(results)

print("\nCalculating ensemble metrics...")
metrics_df = calculate_ensemble_metrics(df)
print("Metrics calculated successfully")

# ============================================================================
# DASHBOARD LAYOUT
# ============================================================================

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("Kinase Virtual Screening Dashboard", style={'text-align': 'center'}),
    
    # Data info
    html.Div([
        html.P(f"Data source: {CSV_FILE}", style={'text-align': 'center', 'font-style': 'italic'}),
        html.P(f"Total kinases: {df['kinase'].nunique()} | "
               f"Total structures: {df.groupby(['kinase', 'structure_id']).ngroups} | "
               f"Total compounds: {df['compound_id'].nunique()}",
               style={'text-align': 'center', 'font-size': '16px'})
    ]),
    
    # Top row: Split heatmaps for AUC and EF metrics
    html.Div([
        html.Div([
            html.H3("AUC Performance"),
            dcc.Graph(id='heatmap-auc')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H3("Enrichment Factors"),
            dcc.Graph(id='heatmap-ef')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Second row: Ranking
    html.Div([
        html.H3("Kinase Ranking by AUC", style={'text-align': 'center'}),
        dcc.Graph(id='ranking-bar')
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

# ============================================================================
# CALLBACKS
# ============================================================================

@app.callback(
    Output('heatmap-auc', 'figure'),
    Input('heatmap-auc', 'id')
)
def update_heatmap_auc(_):
    # Create heatmap for AUC only
    auc_data = metrics_df.set_index('kinase')[['AUC']].T
    
    fig = go.Figure(data=go.Heatmap(
        z=auc_data.values,
        x=auc_data.columns,
        y=['AUC'],
        colorscale='RdBu',
        text=np.round(auc_data.values, 3),
        texttemplate='%{text}',
        textfont={"size": 10},
        colorbar=dict(title="AUC"),
        zmin=0,
        zmax=1
    ))
    
    fig.update_layout(
        title="AUC Values Across Kinases",
        xaxis_title="Kinase",
        yaxis_title="",
        height=200,
        xaxis={'tickangle': -45}
    )
    
    return fig

@app.callback(
    Output('heatmap-ef', 'figure'),
    Input('heatmap-ef', 'id')
)
def update_heatmap_ef(_):
    # Create heatmap for EF metrics
    ef_data = metrics_df.set_index('kinase')[['EF1%', 'EF5%', 'EF10%']].T
    
    fig = go.Figure(data=go.Heatmap(
        z=ef_data.values,
        x=ef_data.columns,
        y=ef_data.index,
        colorscale='Viridis',
        text=np.round(ef_data.values, 1),
        texttemplate='%{text}',
        textfont={"size": 10},
        colorbar=dict(title="EF")
    ))
    
    fig.update_layout(
        title="Enrichment Factors Across Kinases",
        xaxis_title="Kinase",
        yaxis_title="",
        height=300,
        xaxis={'tickangle': -45}
    )
    
    return fig

@app.callback(
    Output('ranking-bar', 'figure'),
    Input('ranking-bar', 'id')
)
def update_ranking(_):
    sorted_df = metrics_df.sort_values('AUC', ascending=True)
    
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
    
    fig.update_layout(
        title="Kinases Ranked by AUC",
        xaxis_title="AUC",
        yaxis_title="Kinase",
        height=600,
        xaxis_range=[0, 1]
    )
    
    return fig

@app.callback(
    Output('kinase-info', 'children'),
    Input('kinase-dropdown', 'value')
)
def update_kinase_info(selected_kinase):
    kinase_metrics = metrics_df[metrics_df['kinase'] == selected_kinase].iloc[0]
    
    return html.Div([
        html.P(
            f"Selected: {selected_kinase}",
            style={'font-weight': 'bold', 'font-size': '18px'}
        ),
        html.P(
            f"Structures: {kinase_metrics['n_structures']} | "
            f"Actives: {kinase_metrics['n_actives']} | "
            f"Decoys: {kinase_metrics['n_decoys']}",
            style={'font-size': '14px'}
        ),
        html.P(
            f"AUC: {kinase_metrics['AUC']:.3f} | "
            f"EF1%: {kinase_metrics['EF1%']:.1f} | "
            f"EF5%: {kinase_metrics['EF5%']:.1f} | "
            f"EF10%: {kinase_metrics['EF10%']:.1f}",
            style={'font-weight': 'bold', 'font-size': '16px', 'color': 'blue'}
        )
    ])

# Keep all other callbacks the same...
@app.callback(
    Output('roc-curve', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_roc_curve(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    fig = go.Figure()
    
    # Calculate ensemble ROC
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
    actives = ensemble_scores[ensemble_scores['compound_type'] == 'active']['docking_score'].values
    decoys = ensemble_scores[ensemble_scores['compound_type'] == 'decoy']['docking_score'].values
    
    y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
    y_scores = -np.concatenate([actives, decoys])
    
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    
    # Ensemble ROC
    fig.add_trace(go.Scatter(
        x=fpr, y=tpr,
        mode='lines',
        name=f'Ensemble (AUC={auc(fpr, tpr):.3f})',
        line=dict(width=3, color='red')
    ))
    
    # Add ROCs for a sample of individual structures
    structure_ids = sorted(kinase_data['structure_id'].unique())
    sample_structures = structure_ids[::max(1, len(structure_ids)//5)][:5]
    
    for struct_id in sample_structures:
        struct_data = kinase_data[kinase_data['structure_id'] == struct_id]
        
        actives = struct_data[struct_data['compound_type'] == 'active']['docking_score'].values
        decoys = struct_data[struct_data['compound_type'] == 'decoy']['docking_score'].values
        
        if len(actives) > 0 and len(decoys) > 0:
            y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
            y_scores = -np.concatenate([actives, decoys])
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            
            fig.add_trace(go.Scatter(
                x=fpr, y=tpr,
                mode='lines',
                name=f'Structure {struct_id} (AUC={auc(fpr, tpr):.3f})',
                line=dict(dash='dash'),
                opacity=0.6
            ))
    
    # Random classifier line
    fig.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1],
        mode='lines',
        name='Random',
        line=dict(dash='dot', color='gray')
    ))
    
    fig.update_layout(
        title=f'ROC Curves for {selected_kinase}',
        xaxis_title='False Positive Rate',
        yaxis_title='True Positive Rate',
        height=400,
        xaxis=dict(constrain='domain'),
        yaxis=dict(scaleanchor="x", scaleratio=1, constrain='domain')
    )
    
    return fig

@app.callback(
    Output('score-distribution', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_score_distribution(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Get ensemble scores
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].agg(['mean', 'std']).reset_index()
    
    fig = go.Figure()
    
    # Box plots for actives and decoys
    for compound_type, color in [('active', 'green'), ('decoy', 'red')]:
        type_data = ensemble_scores[ensemble_scores['compound_type'] == compound_type]
        
        fig.add_trace(go.Box(
            y=type_data['mean'],
            name=f'{compound_type.capitalize()} (n={len(type_data)})',
            marker_color=color,
            boxpoints='outliers'
        ))
    
    fig.update_layout(
        title=f'Ensemble Score Distribution for {selected_kinase}',
        yaxis_title='Docking Score (kcal/mol)',
        showlegend=True,
        height=400
    )
    
    return fig
@app.callback(
    Output('structure-performance', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_structure_performance(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Calculate AUC and EF1% for each structure
    structure_metrics = []
    
    for struct_id in sorted(kinase_data['structure_id'].unique()):
        struct_data = kinase_data[kinase_data['structure_id'] == struct_id]
        
        actives = struct_data[struct_data['compound_type'] == 'active']['docking_score'].values
        decoys = struct_data[struct_data['compound_type'] == 'decoy']['docking_score'].values
        
        if len(actives) > 0 and len(decoys) > 0:
            # Calculate AUC
            y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
            y_scores = -np.concatenate([actives, decoys])
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            auc_score = auc(fpr, tpr)
            
            # Calculate EF1%
            struct_scores = struct_data[['compound_id', 'compound_type', 'docking_score']].copy()
            struct_scores = struct_scores.sort_values('docking_score')
            
            n_actives = len(actives)
            n_total = len(struct_scores)
            top_1_percent = max(1, int(0.01 * n_total))
            actives_in_top_1 = (struct_scores.head(top_1_percent)['compound_type'] == 'active').sum()
            ef_1 = (actives_in_top_1 / top_1_percent) / (n_actives / n_total) if n_actives > 0 else 0
            
            structure_metrics.append({
                'structure_id': struct_id,
                'AUC': auc_score,
                'EF1%': ef_1,
                'n_compounds': len(struct_data)
            })
    
    struct_df = pd.DataFrame(structure_metrics)
    
    # Create subplots with secondary y-axis
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.1,
        subplot_titles=('AUC by Structure', 'EF1% by Structure'),
        row_heights=[0.5, 0.5]
    )
    
    # Plot AUC
    fig.add_trace(
        go.Scatter(
            x=struct_df['structure_id'],
            y=struct_df['AUC'],
            mode='markers+lines',
            name='Structure AUC',
            marker=dict(size=8, color='blue'),
            line=dict(width=1),
            showlegend=True
        ),
        row=1, col=1
    )
    
    # Add ensemble AUC line
    ensemble_auc = metrics_df[metrics_df['kinase'] == selected_kinase]['AUC'].values[0]
    fig.add_hline(
        y=ensemble_auc, 
        line_dash="dash", 
        line_color="red",
        annotation_text=f"Ensemble: {ensemble_auc:.3f}",
        annotation_position="right",
        row=1, col=1
    )
    
    # Add AUC standard deviation band
    mean_auc = struct_df['AUC'].mean()
    std_auc = struct_df['AUC'].std()
    
    fig.add_hrect(
        y0=mean_auc-std_auc, 
        y1=mean_auc+std_auc,
        fillcolor="lightgray", 
        opacity=0.3,
        annotation_text="±1 SD",
        annotation_position="top left",
        row=1, col=1
    )
    
    # Plot EF1%
    fig.add_trace(
        go.Scatter(
            x=struct_df['structure_id'],
            y=struct_df['EF1%'],
            mode='markers+lines',
            name='Structure EF1%',
            marker=dict(size=8, color='green'),
            line=dict(width=1),
            showlegend=True
        ),
        row=2, col=1
    )
    
    # Add ensemble EF1% line
    ensemble_ef1 = metrics_df[metrics_df['kinase'] == selected_kinase]['EF1%'].values[0]
    fig.add_hline(
        y=ensemble_ef1, 
        line_dash="dash", 
        line_color="darkgreen",
        annotation_text=f"Ensemble: {ensemble_ef1:.1f}",
        annotation_position="right",
        row=2, col=1
    )
    
    # Add EF1% standard deviation band
    mean_ef1 = struct_df['EF1%'].mean()
    std_ef1 = struct_df['EF1%'].std()
    
    fig.add_hrect(
        y0=max(0, mean_ef1-std_ef1),  # EF can't be negative
        y1=mean_ef1+std_ef1,
        fillcolor="lightgreen", 
        opacity=0.3,
        annotation_text="±1 SD",
        annotation_position="top left",
        row=2, col=1
    )
    
    # Update layout
    fig.update_xaxes(title_text="Structure ID", row=2, col=1)
    fig.update_yaxes(title_text="AUC", row=1, col=1)
    fig.update_yaxes(title_text="EF1%", row=2, col=1)
    
    fig.update_layout(
        title=f'Structure Performance Variability for {selected_kinase}',
        height=600,
        hovermode='x unified',
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )
    
    # Add annotations with statistics
    fig.add_annotation(
        text=f"AUC: μ={mean_auc:.3f}, σ={std_auc:.3f}",
        xref="paper", yref="paper",
        x=0.02, y=0.95,
        showarrow=False,
        font=dict(size=12),
        bgcolor="white",
        bordercolor="black",
        borderwidth=1
    )
    
    fig.add_annotation(
        text=f"EF1%: μ={mean_ef1:.1f}, σ={std_ef1:.1f}",
        xref="paper", yref="paper",
        x=0.02, y=0.45,
        showarrow=False,
        font=dict(size=12),
        bgcolor="white",
        bordercolor="black",
        borderwidth=1
    )
    
    return fig

@app.callback(
    Output('top-compounds-table', 'children'),
    Input('kinase-dropdown', 'value')
)
def update_top_compounds(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Get ensemble scores with statistics
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
        'docking_score': ['mean', 'std', 'count']
    }).reset_index()
    
    ensemble_scores.columns = ['compound_id', 'compound_type', 'mean_score', 'std_score', 'n_structures']
    ensemble_scores = ensemble_scores.sort_values('mean_score').head(20)
    
    # Create table
    table_header = [
        html.Thead([
            html.Tr([
                html.Th("Rank"),
                html.Th("Compound ID"),
                html.Th("Type"),
                html.Th("Mean Score"),
                html.Th("Std Dev"),
                html.Th("N Structures")
            ])
        ])
    ]
    
    rows = []
    for i, row in ensemble_scores.iterrows():
        row_style = {'backgroundColor': '#90EE90' if row['compound_type'] == 'active' else '#FFB6C1'}
        
        rows.append(
            html.Tr([
                html.Td(i+1),
                html.Td(row['compound_id'][:20] + '...' if len(row['compound_id']) > 20 else row['compound_id']),
                html.Td(row['compound_type']),
                html.Td(f"{row['mean_score']:.2f}"),
                html.Td(f"{row['std_score']:.2f}"),
                html.Td(row['n_structures'])
            ], style=row_style)
        )
    
    table_body = [html.Tbody(rows)]
    
    return html.Table(table_header + table_body, style={
        'width': '100%',
        'text-align': 'center',
        'border': '1px solid black',
        'font-size': '12px'
    })

@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("export-button", "n_clicks"),
    State('kinase-dropdown', 'value'),
    prevent_initial_call=True,
)
def export_kinase_data(n_clicks, selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Get ensemble scores
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
        'docking_score': ['mean', 'std', 'min', 'max', 'count']
    }).reset_index()
    
    ensemble_scores.columns = ['compound_id', 'compound_type', 'mean_score', 'std_score', 
                               'best_score', 'worst_score', 'n_structures']
    
    return dcc.send_data_frame(ensemble_scores.to_csv, f"{selected_kinase}_ensemble_scores.csv")

# ============================================================================
# RUN THE APP
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*50)
    print("Starting dashboard...")
    print("Local access: http://localhost:8050/")
    print("Network access: http://10.10.1.32:8050/")
    print("\nShare this URL with others: http://10.10.1.32:8050/")
    print("Press Ctrl+C to stop")
    print("="*50 + "\n")
    
    # Use 0.0.0.0 to make it accessible from any network interface
    app.run(debug=False, host='0.0.0.0', port=8050)