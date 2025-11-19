#!/usr/bin/env python3
"""
Quick start script to process your DLG files and launch the dashboard
"""

import os
import sys
import argparse

def check_dependencies():
    """Check if required packages are installed"""
    required_packages = {
        'pandas': 'pandas',
        'numpy': 'numpy',
        'dash': 'dash',
        'plotly': 'plotly',
        'sklearn': 'scikit-learn',
        'tqdm': 'tqdm'
    }
    
    missing = []
    for import_name, package_name in required_packages.items():
        try:
            __import__(import_name)
        except ImportError:
            missing.append(package_name)
    
    if missing:
        print("Missing required packages:")
        print(f"Please install: pip install {' '.join(missing)}")
        return False
    return True

def setup_directory_structure():
    """
    Print expected directory structure and check if it exists
    """
    print("""
Expected directory structure:
base_directory/
├── KINASE1/           (e.g., ABL1, AKT1, BRAF, CDK2, etc.)
│   ├── structure_1/   (or any naming with numbers)
│   │   ├── active_compound1.dlg
│   │   ├── active_compound2.dlg
│   │   ├── decoy_compound1.dlg
│   │   └── ...
│   ├── structure_2/
│   └── ...
├── KINASE2/
└── ...

Your compound files should have 'active' or 'decoy' in the filename,
or you'll need to modify the compound_type detection logic.
""")

def create_run_script(base_directory):
    """Create a customized run script"""
    script_content = f'''#!/usr/bin/env python3
"""
Auto-generated script to run your kinase screening dashboard
Generated for: {base_directory}
"""

import os
import sys

# Add the directory containing the dashboard script to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import and run the dashboard
from integrated_dlg_dashboard import app, load_or_process_data, calculate_ensemble_metrics
import integrated_dlg_dashboard as dashboard

# Set your base directory
dashboard.BASE_DIRECTORY = r"{base_directory}"
dashboard.CSV_FILE = "processed_docking_results.csv"

# Load or process data
print("Loading docking data...")
dashboard.df = load_or_process_data(
    base_directory=dashboard.BASE_DIRECTORY, 
    csv_file=dashboard.CSV_FILE
)

# Recalculate metrics
print("Calculating ensemble metrics...")
dashboard.metrics_df = calculate_ensemble_metrics(dashboard.df)

# Print summary
print(f"\\nData Summary:")
print(f"Total kinases: {{dashboard.df['kinase'].nunique()}}")
print(f"Total compounds: {{dashboard.df['compound_id'].nunique()}}")
print(f"Total docking results: {{len(dashboard.df)}}")

# Run the dashboard
print("\\nStarting dashboard at http://127.0.0.1:8050/")
print("Press Ctrl+C to stop")

if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)
'''
    
    with open('run_dashboard.py', 'w') as f:
        f.write(script_content)
    
    print(f"\nCreated 'run_dashboard.py' configured for: {base_directory}")
    print("Run it with: python run_dashboard.py")

def analyze_dlg_directory(base_directory):
    """Quick analysis of DLG directory structure"""
    print(f"\nAnalyzing directory: {base_directory}")
    
    kinase_count = 0
    structure_count = 0
    dlg_count = 0
    
    for kinase_dir in os.listdir(base_directory):
        kinase_path = os.path.join(base_directory, kinase_dir)
        if os.path.isdir(kinase_path):
            kinase_count += 1
            print(f"\nKinase: {kinase_dir}")
            
            structures = []
            for item in os.listdir(kinase_path):
                item_path = os.path.join(kinase_path, item)
                if os.path.isdir(item_path):
                    structures.append(item)
                    structure_count += 1
                    
                    # Count DLG files
                    dlg_files = [f for f in os.listdir(item_path) if f.endswith('.dlg')]
                    dlg_count += len(dlg_files)
            
            print(f"  Structures: {len(structures)}")
            if structures:
                print(f"  Example structures: {structures[:3]}")
    
    print(f"\nSummary:")
    print(f"Total kinases: {kinase_count}")
    print(f"Total structures: {structure_count}")
    print(f"Total DLG files: {dlg_count}")
    
    return kinase_count > 0 and dlg_count > 0

def main():
    parser = argparse.ArgumentParser(description='Process DLG files and launch dashboard')
    parser.add_argument('base_directory', help='Base directory containing kinase folders')
    parser.add_argument('--analyze-only', action='store_true', 
                       help='Only analyze directory structure without processing')
    
    args = parser.parse_args()
    
    # Check if directory exists
    if not os.path.exists(args.base_directory):
        print(f"Error: Directory '{args.base_directory}' does not exist!")
        return 1
    
    # Check dependencies
    if not check_dependencies():
        return 1
    
    # Show expected structure
    setup_directory_structure()
    
    # Analyze directory
    if not analyze_dlg_directory(args.base_directory):
        print("\nError: No valid DLG files found in the expected structure!")
        return 1
    
    if args.analyze_only:
        return 0
    
    # Create run script
    create_run_script(os.path.abspath(args.base_directory))
    
    # Copy the integrated dashboard script if needed
    if not os.path.exists('integrated_dlg_dashboard.py'):
        print("\nNote: You need to save the 'integrated_dlg_dashboard.py' file from the previous artifact!")
    
    print("\n" + "="*50)
    print("Setup complete! Next steps:")
    print("1. Make sure 'integrated_dlg_dashboard.py' is in the current directory")
    print("2. Run: python run_dashboard.py")
    print("3. Open your browser to http://127.0.0.1:8050/")
    print("="*50)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())