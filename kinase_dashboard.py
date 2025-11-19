import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
import plotly.figure_factory as ff

# Initialize the Dash app
app = dash.Dash(__name__)

# Generate sample data (replace with your actual data)
def generate_sample_data():
    """
    Replace this with loading your actual docking results
    Expected format:
    - kinase: kinase name
    - structure_id: 1-50 for each kinase
    - compound_id: compound identifier
    - compound_type: 'active' or 'decoy'
    - docking_score: docking score (more negative = better)
    """
    kinases = [f'KINASE_{i}' for i in range(1, 26)]
    data = []
    
    np.random.seed(42)
    for kinase in kinases:
        # Different kinases have different performance
        performance_offset = np.random.uniform(-2, 2)
        
        for struct_id in range(1, 51):
            # Actives - generally better scores
            for compound_id in range(100):
                score = np.random.normal(-8 + performance_offset, 1.5)
                data.append({
                    'kinase': kinase,
                    'structure_id': struct_id,
                    'compound_id': f'ACTIVE_{compound_id}',
                    'compound_type': 'active',
                    'docking_score': score
                })
            
            # Decoys - generally worse scores
            for compound_id in range(1000):
                score = np.random.normal(-5 + performance_offset, 2)
                data.append({
                    'kinase': kinase,
                    'structure_id': struct_id,
                    'compound_id': f'DECOY_{compound_id}',
                    'compound_type': 'decoy',
                    'docking_score': score
                })
    
    return pd.DataFrame(data)

# Load or generate data
df = generate_sample_data()

# Calculate ensemble metrics for each kinase
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
        
        y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
        y_scores = -np.concatenate([actives, decoys])  # Negative because more negative = better
        
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # Calculate enrichment factors
        all_compounds = ensemble_scores.sort_values('docking_score')
        n_actives = len(actives)
        n_total = len(all_compounds)
        
        # EF at 1%
        top_1_percent = int(0.01 * n_total)
        actives_in_top_1 = (all_compounds.head(top_1_percent)['compound_type'] == 'active').sum()
        ef_1 = (actives_in_top_1 / top_1_percent) / (n_actives / n_total)
        
        # EF at 5%
        top_5_percent = int(0.05 * n_total)
        actives_in_top_5 = (all_compounds.head(top_5_percent)['compound_type'] == 'active').sum()
        ef_5 = (actives_in_top_5 / top_5_percent) / (n_actives / n_total)
        
        results.append({
            'kinase': kinase,
            'AUC': roc_auc,
            'EF1%': ef_1,
            'EF5%': ef_5,
            'n_actives': n_actives,
            'n_decoys': len(decoys)
        })
    
    return pd.DataFrame(results)

metrics_df = calculate_ensemble_metrics(df)

# Define the layout
app.layout = html.Div([
    html.H1("Kinase Virtual Screening Dashboard", style={'text-align': 'center'}),
    
    # Top row: Overview metrics
    html.Div([
        html.Div([
            html.H3("Performance Overview"),
            dcc.Graph(id='heatmap-overview')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H3("Kinase Ranking"),
            dcc.Graph(id='ranking-bar')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Kinase selector
    html.Div([
        html.H3("Select Kinase for Detailed Analysis:"),
        dcc.Dropdown(
            id='kinase-dropdown',
            options=[{'label': k, 'value': k} for k in sorted(df['kinase'].unique())],
            value=df['kinase'].unique()[0]
        )
    ], style={'width': '48%', 'margin': 'auto', 'margin-top': '20px'}),
    
    # Detailed analysis row
    html.Div([
        html.Div([
            html.H4("ROC Curve"),
            dcc.Graph(id='roc-curve')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H4("Score Distribution Across Structures"),
            dcc.Graph(id='score-distribution')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Structure variability analysis
    html.Div([
        html.Div([
            html.H4("Structure Ensemble Performance"),
            dcc.Graph(id='structure-performance')
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.H4("Top Compounds"),
            html.Div(id='top-compounds-table')
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    
    # Cross-kinase analysis
    html.Div([
        html.H3("Cross-Kinase Compound Analysis"),
        html.P("Select a compound to see its scores across all kinases:"),
        dcc.Dropdown(
            id='compound-dropdown',
            options=[],  # Will be populated dynamically
            value=None
        ),
        dcc.Graph(id='compound-cross-kinase')
    ], style={'margin-top': '20px'})
])

# Callback for overview heatmap
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
        textfont={"size": 10}
    ))
    
    fig.update_layout(
        title="Performance Metrics Heatmap",
        xaxis_title="Kinase",
        yaxis_title="Metric",
        height=300
    )
    
    return fig

# Callback for ranking bar chart
@app.callback(
    Output('ranking-bar', 'figure'),
    Input('ranking-bar', 'id')
)
def update_ranking(_):
    sorted_df = metrics_df.sort_values('AUC', ascending=True)
    
    fig = go.Figure(go.Bar(
        x=sorted_df['AUC'],
        y=sorted_df['kinase'],
        orientation='h',
        marker_color='lightblue',
        text=sorted_df['AUC'].round(3),
        textposition='auto'
    ))
    
    fig.update_layout(
        title="Kinases Ranked by AUC",
        xaxis_title="AUC",
        yaxis_title="Kinase",
        height=300
    )
    
    return fig

# Callback for ROC curve
@app.callback(
    Output('roc-curve', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_roc_curve(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Calculate ROC for ensemble and individual structures
    fig = go.Figure()
    
    # Ensemble ROC
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
    actives = ensemble_scores[ensemble_scores['compound_type'] == 'active']['docking_score'].values
    decoys = ensemble_scores[ensemble_scores['compound_type'] == 'decoy']['docking_score'].values
    
    y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
    y_scores = -np.concatenate([actives, decoys])
    
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    
    fig.add_trace(go.Scatter(
        x=fpr, y=tpr,
        mode='lines',
        name=f'Ensemble (AUC={auc(fpr, tpr):.3f})',
        line=dict(width=3, color='red')
    ))
    
    # Add a few individual structure ROCs for comparison
    for struct_id in [1, 10, 25, 50]:
        struct_data = kinase_data[kinase_data['structure_id'] == struct_id]
        actives = struct_data[struct_data['compound_type'] == 'active']['docking_score'].values
        decoys = struct_data[struct_data['compound_type'] == 'decoy']['docking_score'].values
        
        y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
        y_scores = -np.concatenate([actives, decoys])
        
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        
        fig.add_trace(go.Scatter(
            x=fpr, y=tpr,
            mode='lines',
            name=f'Structure {struct_id}',
            line=dict(dash='dash'),
            opacity=0.5
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
        height=400
    )
    
    return fig

# Callback for score distribution
@app.callback(
    Output('score-distribution', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_score_distribution(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    fig = go.Figure()
    
    # Box plot for actives and decoys
    for compound_type, color in [('active', 'green'), ('decoy', 'red')]:
        type_data = kinase_data[kinase_data['compound_type'] == compound_type]
        
        fig.add_trace(go.Box(
            y=type_data['docking_score'],
            name=compound_type.capitalize(),
            marker_color=color,
            boxpoints='outliers'
        ))
    
    fig.update_layout(
        title=f'Score Distribution for {selected_kinase}',
        yaxis_title='Docking Score',
        showlegend=True,
        height=400
    )
    
    return fig

# Callback for structure performance
@app.callback(
    Output('structure-performance', 'figure'),
    Input('kinase-dropdown', 'value')
)
def update_structure_performance(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Calculate AUC for each structure
    structure_aucs = []
    
    for struct_id in range(1, 51):
        struct_data = kinase_data[kinase_data['structure_id'] == struct_id]
        actives = struct_data[struct_data['compound_type'] == 'active']['docking_score'].values
        decoys = struct_data[struct_data['compound_type'] == 'decoy']['docking_score'].values
        
        if len(actives) > 0 and len(decoys) > 0:
            y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
            y_scores = -np.concatenate([actives, decoys])
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            structure_aucs.append({'structure_id': struct_id, 'AUC': auc(fpr, tpr)})
    
    struct_df = pd.DataFrame(structure_aucs)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=struct_df['structure_id'],
        y=struct_df['AUC'],
        mode='markers+lines',
        marker=dict(size=8),
        line=dict(width=1)
    ))
    
    # Add ensemble AUC line
    ensemble_auc = metrics_df[metrics_df['kinase'] == selected_kinase]['AUC'].values[0]
    fig.add_hline(y=ensemble_auc, line_dash="dash", line_color="red",
                  annotation_text=f"Ensemble AUC: {ensemble_auc:.3f}")
    
    fig.update_layout(
        title=f'Individual Structure Performance for {selected_kinase}',
        xaxis_title='Structure ID',
        yaxis_title='AUC',
        height=400
    )
    
    return fig

# Callback for top compounds table
@app.callback(
    [Output('top-compounds-table', 'children'),
     Output('compound-dropdown', 'options')],
    Input('kinase-dropdown', 'value')
)
def update_top_compounds(selected_kinase):
    kinase_data = df[df['kinase'] == selected_kinase]
    
    # Get ensemble scores
    ensemble_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].agg(['mean', 'std']).reset_index()
    ensemble_scores = ensemble_scores.sort_values('mean').head(20)
    
    # Create table
    table_header = [
        html.Thead([
            html.Tr([
                html.Th("Rank"),
                html.Th("Compound ID"),
                html.Th("Type"),
                html.Th("Avg Score"),
                html.Th("Std Dev")
            ])
        ])
    ]
    
    table_body = [
        html.Tbody([
            html.Tr([
                html.Td(i+1),
                html.Td(row['compound_id']),
                html.Td(row['compound_type']),
                html.Td(f"{row['mean']:.2f}"),
                html.Td(f"{row['std']:.2f}")
            ]) for i, row in ensemble_scores.iterrows()
        ])
    ]
    
    table = html.Table(table_header + table_body, style={
        'width': '100%',
        'text-align': 'center',
        'border': '1px solid black'
    })
    
    # Update compound dropdown options
    compound_options = [{'label': cid, 'value': cid} for cid in ensemble_scores['compound_id'].head(10)]
    
    return table, compound_options

# Callback for cross-kinase compound analysis
@app.callback(
    Output('compound-cross-kinase', 'figure'),
    Input('compound-dropdown', 'value')
)
def update_compound_cross_kinase(selected_compound):
    if not selected_compound:
        return go.Figure()
    
    # Get scores for this compound across all kinases
    compound_scores = []
    
    for kinase in df['kinase'].unique():
        kinase_data = df[(df['kinase'] == kinase) & (df['compound_id'] == selected_compound)]
        if not kinase_data.empty:
            avg_score = kinase_data['docking_score'].mean()
            compound_scores.append({'kinase': kinase, 'score': avg_score})
    
    scores_df = pd.DataFrame(compound_scores).sort_values('score')
    
    fig = go.Figure(go.Bar(
        x=scores_df['kinase'],
        y=scores_df['score'],
        marker_color='lightcoral'
    ))
    
    fig.update_layout(
        title=f'Docking Scores for {selected_compound} Across All Kinases',
        xaxis_title='Kinase',
        yaxis_title='Average Docking Score',
        height=300
    )
    
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
