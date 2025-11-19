#!/usr/bin/env python3
"""
Updated Dashboard with three selection strategies:
1. Conformation-based (original): EF calculated across all conformation-compound pairs
2. Best score: Select best conformation per compound  
3. Consensus: Average score per compound
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

# Update with your filename
CSV_FILE = "docking_results_clean_20250917_232813.csv"
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
# CALCULATE METRICS WITH THREE STRATEGIES
# ============================================================================



def calculate_ligand_centric_metrics(df, strategy='best_score'):
    """Calculate AUC and enrichment factors using ligand-centric approach"""
    results = []
    
    for kinase in df['kinase'].unique():
        kinase_data = df[df['kinase'] == kinase]
        
        # Select representative score per compound
        if strategy == 'best_score':
            # Select the best (lowest) scoring conformation for each compound
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].min().reset_index()
        elif strategy == 'consensus':
            # Use average score across conformations
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
        
        # Calculate ROC AUC
        actives = compound_scores[compound_scores['compound_type'] == 'active']['docking_score'].values
        decoys = compound_scores[compound_scores['compound_type'] == 'decoy']['docking_score'].values
        
        if len(actives) > 0 and len(decoys) > 0:
            y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
            y_scores = -np.concatenate([actives, decoys])
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            roc_auc = auc(fpr, tpr)
            
            # Calculate enrichment factors
            all_compounds = compound_scores.sort_values('docking_score')
            n_actives = len(actives)
            n_total = len(all_compounds)
            
            # EF calculations
            ef_results = {}
            for pct in [0.01, 0.05, 0.10]:
                top_n = max(1, int(pct * n_total))
                actives_in_top = (all_compounds.head(top_n)['compound_type'] == 'active').sum()
                ef_results[f'EF{int(pct*100)}%'] = (actives_in_top / top_n) / (n_actives / n_total) if n_actives > 0 else 0
            
            results.append({
                'kinase': kinase,
                'strategy': strategy,
                'AUC': roc_auc,
                **ef_results,
                'n_actives': n_actives,
                'n_decoys': len(decoys),
                'n_structures': kinase_data['structure_id'].nunique(),
                'n_compounds': n_total
            })
    
    return pd.DataFrame(results)

print("\nCalculating metrics with two strategies...")

# Calculate both realistic strategies
best_metrics = calculate_ligand_centric_metrics(df, 'best_score')
consensus_metrics = calculate_ligand_centric_metrics(df, 'consensus')

# Combine metrics
all_metrics = pd.concat([best_metrics, consensus_metrics], ignore_index=True)

print("Metrics calculated successfully")
print(f"Strategies: {all_metrics['strategy'].unique()}")

# ============================================================================
# DASHBOARD LAYOUT
# ============================================================================

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("Kinase Virtual Screening Dashboard - Strategy Comparison", style={'text-align': 'center'}),
    
    # Data info
    html.Div([
        html.P(f"Data source: {CSV_FILE}", style={'text-align': 'center', 'font-style': 'italic'}),
        html.P(f"Total kinases: {df['kinase'].nunique()} | "
               f"Total structures: {df.groupby(['kinase', 'structure_id']).ngroups} | "
               f"Total compounds: {df['compound_id'].nunique()}",
               style={'text-align': 'center', 'font-size': '16px'})
    ]),
    
    # Strategy selector
    html.Div([
        html.H3("Select Analysis Strategy:"),
        dcc.RadioItems(
            id='strategy-selector',
            options=[
                {'label': 'Best Score (Standard VS)', 'value': 'best_score'},
                {'label': 'Consensus (Ensemble)', 'value': 'consensus'}
            ],
            value='best_score',
            labelStyle={'display': 'inline-block', 'margin': '10px'},
            style={'text-align': 'center', 'font-size': '16px'}
        ),
        html.Div(id='strategy-description', style={'text-align': 'center', 'margin': '10px', 'font-style': 'italic'})
    ]),
    
    # Strategy comparison section
    html.Div([
        html.H3("Strategy Comparison", style={'text-align': 'center'}),
        html.Div([
            html.Div([
                dcc.Graph(id='strategy-comparison-auc')
            ], style={'width': '48%', 'display': 'inline-block'}),
            
            html.Div([
                dcc.Graph(id='strategy-comparison-ef')
            ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
        ])
    ]),
    
    # Top row: Split heatmaps for AUC and EF metrics (for selected strategy)
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
    Output('strategy-description', 'children'),
    Input('strategy-selector', 'value')
)
def update_strategy_description(strategy):
    descriptions = {
        'best_score': 'Selects best scoring conformation per compound (standard virtual screening approach)',
        'consensus': 'Uses average score across conformations per compound (ensemble docking approach)'
    }
    return descriptions.get(strategy, '')

@app.callback(
    Output('strategy-comparison-auc', 'figure'),
    Input('strategy-comparison-auc', 'id')
)
def update_strategy_comparison_auc(_):
    # Create comparison plot for AUC across strategies
    fig = go.Figure()
    
    strategies = ['best_score', 'consensus']
    colors = ['red', 'green']
    
    for strategy, color in zip(strategies, colors):
        strategy_data = all_metrics[all_metrics['strategy'] == strategy]
        fig.add_trace(go.Box(
            y=strategy_data['AUC'],
            name=strategy.replace('_', ' ').title(),
            marker_color=color,
            boxpoints='all'
        ))
    
    fig.update_layout(
        title="AUC Distribution by Strategy",
        yaxis_title="AUC",
        height=400
    )
    
    return fig

@app.callback(
    Output('strategy-comparison-ef', 'figure'),
    Input('strategy-comparison-ef', 'id')
)
def update_strategy_comparison_ef(_):
    # Create comparison plot for EF1% across strategies
    fig = go.Figure()
    
    strategies = ['best_score', 'consensus']
    colors = ['red', 'green']
    
    for strategy, color in zip(strategies, colors):
        strategy_data = all_metrics[all_metrics['strategy'] == strategy]
        fig.add_trace(go.Box(
            y=strategy_data['EF1%'],
            name=strategy.replace('_', ' ').title(),
            marker_color=color,
            boxpoints='all'
        ))
    
    fig.update_layout(
        title="EF1% Distribution by Strategy",
        yaxis_title="EF1%",
        height=400
    )
    
    return fig

@app.callback(
    Output('heatmap-auc', 'figure'),
    Input('strategy-selector', 'value')
)
def update_heatmap_auc(selected_strategy):
    # Create heatmap for AUC only for selected strategy
    strategy_data = all_metrics[all_metrics['strategy'] == selected_strategy]
    auc_data = strategy_data.set_index('kinase')[['AUC']].T
    
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
        title=f"AUC Values - {selected_strategy.replace('_', ' ').title()} Strategy",
        xaxis_title="Kinase",
        yaxis_title="",
        height=200,
        xaxis={'tickangle': -45}
    )
    
    return fig

@app.callback(
    Output('heatmap-ef', 'figure'),
    Input('strategy-selector', 'value')
)
def update_heatmap_ef(selected_strategy):
    # Create heatmap for EF metrics for selected strategy
    strategy_data = all_metrics[all_metrics['strategy'] == selected_strategy]
    ef_data = strategy_data.set_index('kinase')[['EF1%', 'EF5%', 'EF10%']].T
    
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
        title=f"Enrichment Factors - {selected_strategy.replace('_', ' ').title()} Strategy",
        xaxis_title="Kinase",
        yaxis_title="",
        height=300,
        xaxis={'tickangle': -45}
    )
    
    return fig

@app.callback(
    Output('ranking-bar', 'figure'),
    Input('strategy-selector', 'value')
)
def update_ranking(selected_strategy):
    strategy_data = all_metrics[all_metrics['strategy'] == selected_strategy]
    sorted_df = strategy_data.sort_values('AUC', ascending=True)
    
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
        title=f"Kinases Ranked by AUC - {selected_strategy.replace('_', ' ').title()} Strategy",
        xaxis_title="AUC",
        yaxis_title="Kinase",
        height=600,
        xaxis_range=[0, 1]
    )
    
    return fig

@app.callback(
    Output('kinase-info', 'children'),
    [Input('kinase-dropdown', 'value'),
     Input('strategy-selector', 'value')]
)
def update_kinase_info(selected_kinase, selected_strategy):
    strategy_data = all_metrics[all_metrics['strategy'] == selected_strategy]
    kinase_metrics = strategy_data[strategy_data['kinase'] == selected_kinase]
    
    if len(kinase_metrics) == 0:
        return html.P("No data available for this kinase-strategy combination")
    
    kinase_metrics = kinase_metrics.iloc[0]
    
    return html.Div([
        html.P(
            f"Selected: {selected_kinase} | Strategy: {selected_strategy.replace('_', ' ').title()}",
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

# Keep existing detailed analysis callbacks but update them to use selected strategy
@app.callback(
    Output('roc-curve', 'figure'),
    [Input('kinase-dropdown', 'value'),
     Input('strategy-selector', 'value')]
)
def update_roc_curve(selected_kinase, selected_strategy):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    fig = go.Figure()
    
    # Calculate ROC for selected strategy
    if selected_strategy == 'conformation_based':
        # Use all conformation-compound pairs
        actives = kinase_data[kinase_data['compound_type'] == 'active']['docking_score'].values
        decoys = kinase_data[kinase_data['compound_type'] == 'decoy']['docking_score'].values
    elif selected_strategy == 'best_score':
        # Use best score per compound
        compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].min().reset_index()
        actives = compound_scores[compound_scores['compound_type'] == 'active']['docking_score'].values
        decoys = compound_scores[compound_scores['compound_type'] == 'decoy']['docking_score'].values
    else:  # consensus
        # Use average score per compound
        compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
        actives = compound_scores[compound_scores['compound_type'] == 'active']['docking_score'].values
        decoys = compound_scores[compound_scores['compound_type'] == 'decoy']['docking_score'].values
    
    y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
    y_scores = -np.concatenate([actives, decoys])
    
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    
    # Main ROC curve
    fig.add_trace(go.Scatter(
        x=fpr, y=tpr,
        mode='lines',
        name=f'{selected_strategy.replace("_", " ").title()} (AUC={auc(fpr, tpr):.3f})',
        line=dict(width=3, color='red')
    ))
    
    # Add comparison with other strategy
    other_strategy = 'consensus' if selected_strategy == 'best_score' else 'best_score'
    other_color = 'green' if other_strategy == 'consensus' else 'blue'
    
    if other_strategy == 'best_score':
        comp_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].min().reset_index()
    else:
        comp_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
    
    comp_actives = comp_scores[comp_scores['compound_type'] == 'active']['docking_score'].values
    comp_decoys = comp_scores[comp_scores['compound_type'] == 'decoy']['docking_score'].values
    
    comp_y_true = np.concatenate([np.ones(len(comp_actives)), np.zeros(len(comp_decoys))])
    comp_y_scores = -np.concatenate([comp_actives, comp_decoys])
    
    comp_fpr, comp_tpr, _ = roc_curve(comp_y_true, comp_y_scores)
    
    fig.add_trace(go.Scatter(
        x=comp_fpr, y=comp_tpr,
        mode='lines',
        name=f'{other_strategy.replace("_", " ").title()} (AUC={auc(comp_fpr, comp_tpr):.3f})',
        line=dict(dash='dash', color=other_color),
        opacity=0.7
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
    [Input('kinase-dropdown', 'value'),
     Input('strategy-selector', 'value')]
)
def update_score_distribution(selected_kinase, selected_strategy):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    fig = go.Figure()
    
    if selected_strategy == 'conformation_based':
        # Show all conformation-compound pairs
        for compound_type, color in [('active', 'green'), ('decoy', 'red')]:
            type_data = kinase_data[kinase_data['compound_type'] == compound_type]
            fig.add_trace(go.Box(
                y=type_data['docking_score'],
                name=f'{compound_type.capitalize()} (n={len(type_data)})',
                marker_color=color,
                boxpoints='outliers'
            ))
    else:
        # Show compound-level scores
        if selected_strategy == 'best_score':
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].min().reset_index()
        else:  # consensus
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
        
        for compound_type, color in [('active', 'green'), ('decoy', 'red')]:
            type_data = compound_scores[compound_scores['compound_type'] == compound_type]
            fig.add_trace(go.Box(
                y=type_data['docking_score'],
                name=f'{compound_type.capitalize()} (n={len(type_data)})',
                marker_color=color,
                boxpoints='outliers'
            ))
    
    fig.update_layout(
        title=f'Score Distribution for {selected_kinase} - {selected_strategy.replace("_", " ").title()}',
        yaxis_title='Docking Score (kcal/mol)',
        showlegend=True,
        height=400
    )
    
    return fig

# Keep the remaining callbacks (structure-performance, top-compounds-table, export) the same as they work with the current selected strategy

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
    
    # Add ensemble AUC line (using consensus strategy)
    consensus_data = all_metrics[(all_metrics['kinase'] == selected_kinase) & 
                                (all_metrics['strategy'] == 'consensus')]
    if len(consensus_data) > 0:
        ensemble_auc = consensus_data['AUC'].values[0]
        fig.add_hline(
            y=ensemble_auc, 
            line_dash="dash", 
            line_color="red",
            annotation_text=f"Consensus: {ensemble_auc:.3f}",
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
    
    # Add ensemble EF1% line (using consensus strategy)
    if len(consensus_data) > 0:
        ensemble_ef1 = consensus_data['EF1%'].values[0]
        fig.add_hline(
            y=ensemble_ef1, 
            line_dash="dash", 
            line_color="darkgreen",
            annotation_text=f"Consensus: {ensemble_ef1:.1f}",
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
    [Input('kinase-dropdown', 'value'),
     Input('strategy-selector', 'value')]
)
def update_top_compounds(selected_kinase, selected_strategy):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Get scores based on selected strategy
    if selected_strategy == 'conformation_based':
        # Show top scoring conformation-compound pairs
        top_pairs = kinase_data[['compound_id', 'compound_type', 'docking_score', 'structure_id']].copy()
        top_pairs = top_pairs.sort_values('docking_score').head(20)
        
        # Create table for conformation pairs
        table_header = [
            html.Thead([
                html.Tr([
                    html.Th("Rank"),
                    html.Th("Compound ID"),
                    html.Th("Type"),
                    html.Th("Score"),
                    html.Th("Structure ID")
                ])
            ])
        ]
        
        rows = []
        for i, (_, row) in enumerate(top_pairs.iterrows()):
            row_style = {'backgroundColor': '#90EE90' if row['compound_type'] == 'active' else '#FFB6C1'}
            
            rows.append(
                html.Tr([
                    html.Td(i+1),
                    html.Td(row['compound_id'][:15] + '...' if len(row['compound_id']) > 15 else row['compound_id']),
                    html.Td(row['compound_type']),
                    html.Td(f"{row['docking_score']:.2f}"),
                    html.Td(str(row['structure_id'])[:10])
                ], style=row_style)
            )
        
    else:
        # Get compound-level scores with statistics
        if selected_strategy == 'best_score':
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
                'docking_score': ['min', 'std', 'count']
            }).reset_index()
            score_col = 'min'
        else:  # consensus
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
                'docking_score': ['mean', 'std', 'count']
            }).reset_index()
            score_col = 'mean'
        
        compound_scores.columns = ['compound_id', 'compound_type', 'score', 'std_score', 'n_structures']
        compound_scores = compound_scores.sort_values('score').head(20)
        
        # Create table for compounds
        table_header = [
            html.Thead([
                html.Tr([
                    html.Th("Rank"),
                    html.Th("Compound ID"),
                    html.Th("Type"),
                    html.Th(f"{score_col.title()} Score"),
                    html.Th("Std Dev"),
                    html.Th("N Structures")
                ])
            ])
        ]
        
        rows = []
        for i, (_, row) in enumerate(compound_scores.iterrows()):
            row_style = {'backgroundColor': '#90EE90' if row['compound_type'] == 'active' else '#FFB6C1'}
            
            rows.append(
                html.Tr([
                    html.Td(i+1),
                    html.Td(row['compound_id'][:15] + '...' if len(row['compound_id']) > 15 else row['compound_id']),
                    html.Td(row['compound_type']),
                    html.Td(f"{row['score']:.2f}"),
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
    [Input("export-button", "n_clicks"),
     State('kinase-dropdown', 'value'),
     State('strategy-selector', 'value')],
    prevent_initial_call=True,
)
def export_kinase_data(n_clicks, selected_kinase, selected_strategy):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    if selected_strategy == 'conformation_based':
        # Export all conformation-compound pairs
        export_data = kinase_data[['compound_id', 'compound_type', 'docking_score', 'structure_id']].copy()
        filename = f"{selected_kinase}_{selected_strategy}_data.csv"
    else:
        # Export compound-level scores
        if selected_strategy == 'best_score':
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
                'docking_score': ['min', 'std', 'max', 'count']
            }).reset_index()
            compound_scores.columns = ['compound_id', 'compound_type', 'best_score', 'std_score', 'worst_score', 'n_structures']
        else:  # consensus
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type']).agg({
                'docking_score': ['mean', 'std', 'min', 'max', 'count']
            }).reset_index()
            compound_scores.columns = ['compound_id', 'compound_type', 'mean_score', 'std_score', 'best_score', 'worst_score', 'n_structures']
        
        export_data = compound_scores
        filename = f"{selected_kinase}_{selected_strategy}_compounds.csv"
    
    return dcc.send_data_frame(export_data.to_csv, filename)

# ============================================================================
# RUN THE APP
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*50)
    print("Starting enhanced dashboard with strategy comparison...")
    print("Available strategies:")
    print("  1. Conformation-based: EF across all pairs (not realistic for VS)")
    print("  2. Best score: Best conformation per compound (standard VS)")
    print("  3. Consensus: Average score per compound (ensemble approach)")
    print("\nLocal access: http://localhost:8050/")
    print("Network access: http://10.10.1.32:8050/")
    print("\nShare this URL with others: http://10.10.1.32:8050/")
    print("Press Ctrl+C to stop")
    print("="*50 + "\n")
    
    # Use 0.0.0.0 to make it accessible from any network interface
    app.run(debug=True, host='0.0.0.0', port=8050)
