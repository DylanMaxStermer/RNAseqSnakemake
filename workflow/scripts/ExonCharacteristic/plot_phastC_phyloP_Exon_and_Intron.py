#!/usr/bin/env python3
"""
Plot PhastCons or PhyloP conservation for multiple datasets.
Continuous profile across all regions - simple and clean.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def detect_metric(df):
    """
    Detect which conservation metric is available in the dataframe.
    
    Returns:
        tuple: (metric_name, mean_col, median_col, sd_col) or None if not found
    """
    if 'phastcons_mean' in df.columns:
        return ('PhastCons', 'phastcons_mean', 'phastcons_median', 'phastcons_sd')
    elif 'phylop_mean' in df.columns:
        return ('PhyloP', 'phylop_mean', 'phylop_median', 'phylop_sd')
    elif 'mean' in df.columns:
        # Backward compatibility with old format
        return ('PhastCons', 'mean', 'median', 'sd')
    else:
        return None

def plot_conservation_comparison(summary_files, labels, output_prefix, metric='auto', output_format='png', statistic='mean'):
    """
    Plot conservation profiles for multiple datasets on the same plot.
    
    Args:
        summary_files: List of paths to _summary.tsv files
        labels: List of labels for each dataset
        output_prefix: Prefix for output plot files
        metric: 'auto', 'phastcons', or 'phylop'
        output_format: 'png', 'pdf', or 'svg' (default: png)
        statistic: 'mean' or 'median' (default: mean)
    """
    
    if len(summary_files) != len(labels):
        print("ERROR: Number of input files must match number of labels")
        return
    
    # Colors for different datasets
    dataset_colors = ['#E63946', '#2A9D8F', '#F4A261', '#264653', '#E76F51', '#8338EC']
    
    # Load all datasets
    datasets = []
    metric_name = None
    value_col = None
    
    for i, (file, label) in enumerate(zip(summary_files, labels)):
        df = pd.read_csv(file, sep='\t', comment='#')
        
        # Detect or validate metric
        if metric == 'auto':
            detected = detect_metric(df)
            if detected is None:
                print(f"ERROR: Could not detect conservation metric in {file}")
                return
            if metric_name is None:
                metric_name, mean_col, median_col, sd_col = detected
                print(f"Auto-detected metric: {metric_name}")
            elif detected[0] != metric_name:
                print(f"WARNING: Mixed metrics detected. Using {metric_name} for all files.")
        elif metric.lower() == 'phastcons':
            if 'phastcons_mean' in df.columns:
                metric_name, mean_col, median_col, sd_col = 'PhastCons', 'phastcons_mean', 'phastcons_median', 'phastcons_sd'
            elif 'mean' in df.columns:
                metric_name, mean_col, median_col, sd_col = 'PhastCons', 'mean', 'median', 'sd'
            else:
                print(f"ERROR: PhastCons columns not found in {file}")
                return
        elif metric.lower() == 'phylop':
            if 'phylop_mean' in df.columns:
                metric_name, mean_col, median_col, sd_col = 'PhyloP', 'phylop_mean', 'phylop_median', 'phylop_sd'
            else:
                print(f"ERROR: PhyloP columns not found in {file}")
                return
        
        # Select which column to use based on statistic parameter
        if statistic == 'mean':
            value_col = mean_col
        elif statistic == 'median':
            value_col = median_col
        else:
            print(f"ERROR: Invalid statistic '{statistic}'. Must be 'mean' or 'median'")
            return
        
        datasets.append({
            'df': df,
            'label': label,
            'color': dataset_colors[i % len(dataset_colors)],
            'value_col': value_col,
            'mean_col': mean_col,
            'median_col': median_col
        })
        print(f"Loaded {label}: {len(df)} positions, {statistic}={df[value_col].mean():.4f}")
    
    # Check that all datasets have same number of positions
    n_positions = len(datasets[0]['df'])
    for dataset in datasets[1:]:
        if len(dataset['df']) != n_positions:
            print(f"WARNING: {dataset['label']} has {len(dataset['df'])} positions, expected {n_positions}")
    
    # Determine max length for x-axis
    max_positions = max([len(d['df']) for d in datasets])
    
    # ====================================================================
    # Create main comparison plot
    # ====================================================================
    fig, ax = plt.subplots(figsize=(20, 8))
    
    # Add region background shading first (if available)
    df_ref = datasets[0]['df']
    if 'region' in df_ref.columns:
        region_colors_bg = {
            "5'_intron": '#E8F4F8',
            'first_exon_edge': '#FFE5E5',
            'last_exon_edge': '#FFF5E5',
            "3'_intron": '#E5F5F0'
        }
        
        # Add vertical spans for each region
        current_region = df_ref.iloc[0]['region']
        region_start = 0
        
        for i in range(1, len(df_ref)):
            if df_ref.iloc[i]['region'] != current_region:
                # Region changed, draw the previous region
                region_end = i
                if current_region in region_colors_bg:
                    ax.axvspan(region_start, region_end, 
                             alpha=0.1, color=region_colors_bg[current_region], 
                             zorder=0)
                
                # Add vertical line at boundary
                ax.axvline(x=i, color='gray', linestyle=':', linewidth=1.5, 
                          alpha=0.5, zorder=1)
                
                region_start = i
                current_region = df_ref.iloc[i]['region']
        
        # Draw the last region
        if current_region in region_colors_bg:
            ax.axvspan(region_start, len(df_ref), 
                     alpha=0.1, color=region_colors_bg[current_region], 
                     zorder=0)
    
    # Plot each dataset with its own x-positions  
    for dataset in datasets:
        df = dataset['df']
        color = dataset['color']
        label = dataset['label']
        value_col = dataset['value_col']
        
        # Create x-positions for THIS dataset
        x_positions = np.arange(len(df))
        
        # Plot line (mean or median based on user selection)
        ax.plot(x_positions, df[value_col], linewidth=2.5, label=label,
                color=color, zorder=3, alpha=0.9)
    
    # Add region labels at top (if available)
    if 'region' in datasets[0]['df'].columns:
        y_max = ax.get_ylim()[1]
        
        for region_name in ["5'_intron", 'first_exon_edge', 'last_exon_edge', "3'_intron"]:
            region_data = df_ref[df_ref['region'] == region_name]
            if len(region_data) > 0:
                region_indices = region_data.index
                mid_point = (region_indices[0] + region_indices[-1]) / 2
                
                # Clean up label
                label_text = region_name.replace('_', ' ').replace("'", "'")
                
                ax.text(mid_point, y_max * 0.97, label_text, 
                       ha='center', va='top', fontsize=11, fontweight='bold',
                       bbox=dict(boxstyle='round', facecolor='white', 
                                alpha=0.9, edgecolor='gray', linewidth=1))
    
    # Formatting
    ax.set_xlabel('Position Index', fontsize=13, fontweight='bold')
    statistic_label = statistic.capitalize()
    ax.set_ylabel(f'{metric_name} Conservation Score ({statistic_label})', fontsize=13, fontweight='bold')
    ax.set_title(f'{metric_name} Conservation ({statistic_label}): Comparison Across Datasets', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=12, framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    
    # Set y-limits based on metric
    if metric_name == 'PhyloP':
        # PhyloP can be negative (negative selection)
        max_val = max([d['df'][d['value_col']].max() for d in datasets])
        min_val = min([d['df'][d['value_col']].min() for d in datasets])
        y_range = max_val - min_val
        ax.set_ylim(bottom=min_val - y_range*0.05, top=max_val + y_range*0.05)
    else:
        # PhastCons is 0-1
        ax.set_ylim(bottom=0, top=min(1.0, max([d['df'][d['value_col']].max() for d in datasets]) * 1.05))
    
    ax.set_xlim(left=0, right=max_positions) 
    
    plt.tight_layout()
    
    # Save in specified format
    output_file = f"{output_prefix}_comparison_{statistic}.{output_format}"
    
    if output_format == 'png':
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    elif output_format == 'pdf':
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
    elif output_format == 'svg':
        plt.savefig(output_file, format='svg', bbox_inches='tight')
    
    print(f"\nSaved comparison plot: {output_file}")
    
    plt.close()
    
    # ====================================================================
    # Print comparison statistics
    # ====================================================================
    print("\n" + "="*70)
    print(f"COMPARISON STATISTICS ({metric_name} - {statistic_label})")
    print("="*70)
    
    # Overall statistics
    print(f"\nOverall {statistic} conservation:")
    for dataset in datasets:
        overall_val = dataset['df'][dataset['value_col']].mean()
        print(f"  {dataset['label']}: {overall_val:.4f}")
    
    # By region statistics
    if 'region' in datasets[0]['df'].columns:
        print(f"\n{statistic_label} conservation by region:")
        for region_name in ["5'_intron", 'first_exon_edge', 'last_exon_edge', "3'_intron"]:
            print(f"\n  {region_name}:")
            for dataset in datasets:
                region_data = dataset['df'][dataset['df']['region'] == region_name]
                if len(region_data) > 0:
                    val = region_data[dataset['value_col']].mean()
                    print(f"    {dataset['label']}: {val:.4f}")
    
    # Differences
    if len(datasets) == 2:
        print("\n" + "="*70)
        print("PAIRWISE DIFFERENCE")
        print("="*70)
        
        df1 = datasets[0]['df']
        df2 = datasets[1]['df']
        value_col1 = datasets[0]['value_col']
        value_col2 = datasets[1]['value_col']
        
        if len(df1) == len(df2):
            diff = df1[value_col1] - df2[value_col2]
            print(f"\n{datasets[0]['label']} - {datasets[1]['label']}:")
            print(f"  Mean difference: {diff.mean():.4f}")
            print(f"  Max difference: {diff.max():.4f} at position {diff.idxmax()}")
            print(f"  Min difference: {diff.min():.4f} at position {diff.idxmin()}")
            
            # Statistical test
            from scipy import stats as scipy_stats
            t_stat, p_value = scipy_stats.ttest_rel(df1[value_col1], df2[value_col2])
            print(f"\n  Paired t-test:")
            print(f"    t-statistic: {t_stat:.4f}")
            print(f"    p-value: {p_value:.4e}")
            print(f"    Significant: {'Yes' if p_value < 0.05 else 'No'} (α=0.05)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Compare conservation profiles across multiple datasets (PhastCons or PhyloP)',
        epilog="""
Example usage:
  # Auto-detect metric, plot mean values, save as PNG (default)
  python plot_exon_intron.py \\
      -i output/constitutive_summary.tsv output/alternative_summary.tsv \\
      -l "Constitutive" "Alternative" \\
      -o output/plots/comparison
  
  # Plot median values instead of mean
  python plot_exon_intron.py \\
      -i file1.tsv file2.tsv \\
      -l "Dataset1" "Dataset2" \\
      -o output/plots/comparison \\
      -s median
  
  # Save as PDF with median values
  python plot_exon_intron.py \\
      -i file1.tsv file2.tsv file3.tsv \\
      -l "Dataset1" "Dataset2" "Dataset3" \\
      -o output/plots/comparison \\
      -s median -f pdf
  
  # Explicitly specify PhyloP with median, save as SVG
  python plot_exon_intron.py \\
      -i file1.tsv file2.tsv \\
      -l "Dataset1" "Dataset2" \\
      -o output/plots/comparison \\
      -m phylop -s median -f svg
        """
    )
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input _summary.tsv files (space-separated)')
    parser.add_argument('-l', '--labels', nargs='+', required=True,
                       help='Labels for each dataset (space-separated, same order as inputs)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for plot files')
    parser.add_argument('-m', '--metric', default='auto', choices=['auto', 'phastcons', 'phylop'],
                       help='Conservation metric to plot (default: auto-detect)')
    parser.add_argument('-s', '--statistic', default='mean', choices=['mean', 'median'],
                       help='Statistic to plot: mean or median (default: mean)')
    parser.add_argument('-f', '--format', default='png', choices=['png', 'pdf', 'svg'],
                       help='Output file format (default: png)')
    
    args = parser.parse_args()
    
    print("="*70)
    print("Conservation Profile Comparison")
    print("="*70)
    print(f"Number of datasets: {len(args.input)}")
    for i, (file, label) in enumerate(zip(args.input, args.labels)):
        print(f"  {i+1}. {label}: {file}")
    print(f"Output prefix: {args.output}")
    print(f"Metric: {args.metric}")
    print(f"Statistic: {args.statistic}")
    print(f"Format: {args.format}")
    print("="*70 + "\n")
    
    plot_conservation_comparison(args.input, args.labels, args.output, args.metric, args.format, args.statistic)
    
    print("\n" + "="*70)
    print("COMPARISON COMPLETE")
    print("="*70)

# with mean
#################################################################################################################################################################

#!/usr/bin/env python3
# """
# Plot PhastCons or PhyloP conservation for multiple datasets.
# Continuous profile across all regions - simple and clean.
# """

# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import argparse

# def detect_metric(df):
#     """
#     Detect which conservation metric is available in the dataframe.
    
#     Returns:
#         tuple: (metric_name, mean_col, median_col, sd_col) or None if not found
#     """
#     if 'phastcons_mean' in df.columns:
#         return ('PhastCons', 'phastcons_mean', 'phastcons_median', 'phastcons_sd')
#     elif 'phylop_mean' in df.columns:
#         return ('PhyloP', 'phylop_mean', 'phylop_median', 'phylop_sd')
#     elif 'mean' in df.columns:
#         # Backward compatibility with old format
#         return ('PhastCons', 'mean', 'median', 'sd')
#     else:
#         return None

# def plot_conservation_comparison(summary_files, labels, output_prefix, metric='auto', output_format='png'):
#     """
#     Plot conservation profiles for multiple datasets on the same plot.
    
#     Args:
#         summary_files: List of paths to _summary.tsv files
#         labels: List of labels for each dataset
#         output_prefix: Prefix for output plot files
#         metric: 'auto', 'phastcons', or 'phylop'
#         output_format: 'png', 'pdf', or 'svg' (default: png)
#     """
    
#     if len(summary_files) != len(labels):
#         print("ERROR: Number of input files must match number of labels")
#         return
    
#     # Colors for different datasets
#     dataset_colors = ['#E63946', '#2A9D8F', '#F4A261', '#264653', '#E76F51', '#8338EC']
    
#     # Load all datasets
#     datasets = []
#     metric_name = None
#     mean_col = None
    
#     for i, (file, label) in enumerate(zip(summary_files, labels)):
#         df = pd.read_csv(file, sep='\t', comment='#')
        
#         # Detect or validate metric
#         if metric == 'auto':
#             detected = detect_metric(df)
#             if detected is None:
#                 print(f"ERROR: Could not detect conservation metric in {file}")
#                 return
#             if metric_name is None:
#                 metric_name, mean_col, median_col, sd_col = detected
#                 print(f"Auto-detected metric: {metric_name}")
#             elif detected[0] != metric_name:
#                 print(f"WARNING: Mixed metrics detected. Using {metric_name} for all files.")
#         elif metric.lower() == 'phastcons':
#             if 'phastcons_mean' in df.columns:
#                 metric_name, mean_col, median_col, sd_col = 'PhastCons', 'phastcons_mean', 'phastcons_median', 'phastcons_sd'
#             elif 'mean' in df.columns:
#                 metric_name, mean_col, median_col, sd_col = 'PhastCons', 'mean', 'median', 'sd'
#             else:
#                 print(f"ERROR: PhastCons columns not found in {file}")
#                 return
#         elif metric.lower() == 'phylop':
#             if 'phylop_mean' in df.columns:
#                 metric_name, mean_col, median_col, sd_col = 'PhyloP', 'phylop_mean', 'phylop_median', 'phylop_sd'
#             else:
#                 print(f"ERROR: PhyloP columns not found in {file}")
#                 return
        
#         datasets.append({
#             'df': df,
#             'label': label,
#             'color': dataset_colors[i % len(dataset_colors)],
#             'mean_col': mean_col
#         })
#         print(f"Loaded {label}: {len(df)} positions, mean={df[mean_col].mean():.4f}")
    
#     # Check that all datasets have same number of positions
#     n_positions = len(datasets[0]['df'])
#     for dataset in datasets[1:]:
#         if len(dataset['df']) != n_positions:
#             print(f"WARNING: {dataset['label']} has {len(dataset['df'])} positions, expected {n_positions}")
    
#     # Determine max length for x-axis
#     max_positions = max([len(d['df']) for d in datasets])
    
#     # ====================================================================
#     # Create main comparison plot
#     # ====================================================================
#     fig, ax = plt.subplots(figsize=(20, 8))
    
#     # Add region background shading first (if available)
#     df_ref = datasets[0]['df']
#     if 'region' in df_ref.columns:
#         region_colors_bg = {
#             "5'_intron": '#E8F4F8',
#             'first_exon_edge': '#FFE5E5',
#             'last_exon_edge': '#FFF5E5',
#             "3'_intron": '#E5F5F0'
#         }
        
#         # Add vertical spans for each region
#         current_region = df_ref.iloc[0]['region']
#         region_start = 0
        
#         for i in range(1, len(df_ref)):
#             if df_ref.iloc[i]['region'] != current_region:
#                 # Region changed, draw the previous region
#                 region_end = i
#                 if current_region in region_colors_bg:
#                     ax.axvspan(region_start, region_end, 
#                              alpha=0.1, color=region_colors_bg[current_region], 
#                              zorder=0)
                
#                 # Add vertical line at boundary
#                 ax.axvline(x=i, color='gray', linestyle=':', linewidth=1.5, 
#                           alpha=0.5, zorder=1)
                
#                 region_start = i
#                 current_region = df_ref.iloc[i]['region']
        
#         # Draw the last region
#         if current_region in region_colors_bg:
#             ax.axvspan(region_start, len(df_ref), 
#                      alpha=0.1, color=region_colors_bg[current_region], 
#                      zorder=0)
    
#     # Plot each dataset
#     # Plot each dataset with its own x-positions  
#     for dataset in datasets:
#         df = dataset['df']
#         color = dataset['color']
#         label = dataset['label']
#         mean_col = dataset['mean_col']
        
#         # Create x-positions for THIS dataset
#         x_positions = np.arange(len(df))
        
#         # Plot mean line
#         ax.plot(x_positions, df[mean_col], linewidth=2.5, label=label,
#                 color=color, zorder=3, alpha=0.9)
    
#     # Add region labels at top (if available)
#     if 'region' in datasets[0]['df'].columns:
#         y_max = ax.get_ylim()[1]
        
#         for region_name in ["5'_intron", 'first_exon_edge', 'last_exon_edge', "3'_intron"]:
#             region_data = df_ref[df_ref['region'] == region_name]
#             if len(region_data) > 0:
#                 region_indices = region_data.index
#                 mid_point = (region_indices[0] + region_indices[-1]) / 2
                
#                 # Clean up label
#                 label_text = region_name.replace('_', ' ').replace("'", "'")
                
#                 ax.text(mid_point, y_max * 0.97, label_text, 
#                        ha='center', va='top', fontsize=11, fontweight='bold',
#                        bbox=dict(boxstyle='round', facecolor='white', 
#                                 alpha=0.9, edgecolor='gray', linewidth=1))
    
#     # Formatting
#     ax.set_xlabel('Position Index', fontsize=13, fontweight='bold')
#     ax.set_ylabel(f'{metric_name} Conservation Score', fontsize=13, fontweight='bold')
#     ax.set_title(f'{metric_name} Conservation: Comparison Across Datasets', 
#                  fontsize=14, fontweight='bold')
#     ax.legend(loc='upper right', fontsize=12, framealpha=0.95)
#     ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    
#     # Set y-limits based on metric
#     if metric_name == 'PhyloP':
#         # PhyloP can be negative (negative selection)
#         max_val = max([d['df'][d['mean_col']].max() for d in datasets])
#         min_val = min([d['df'][d['mean_col']].min() for d in datasets])
#         y_range = max_val - min_val
#         ax.set_ylim(bottom=min_val - y_range*0.05, top=max_val + y_range*0.05)
#     else:
#         # PhastCons is 0-1
#         ax.set_ylim(bottom=0, top=min(1.0, max([d['df'][d['mean_col']].max() for d in datasets]) * 1.05))
    
#     ax.set_xlim(left=0, right=max_positions) 
    
#     plt.tight_layout()
    
#     # Save in specified format
#     output_file = f"{output_prefix}_comparison.{output_format}"
    
#     if output_format == 'png':
#         plt.savefig(output_file, dpi=300, bbox_inches='tight')
#     elif output_format == 'pdf':
#         plt.savefig(output_file, format='pdf', bbox_inches='tight')
#     elif output_format == 'svg':
#         plt.savefig(output_file, format='svg', bbox_inches='tight')
    
#     print(f"\nSaved comparison plot: {output_file}")
    
#     plt.close()
    
#     # ====================================================================
#     # Print comparison statistics
#     # ====================================================================
#     print("\n" + "="*70)
#     print(f"COMPARISON STATISTICS ({metric_name})")
#     print("="*70)
    
#     # Overall statistics
#     print("\nOverall mean conservation:")
#     for dataset in datasets:
#         overall_mean = dataset['df'][dataset['mean_col']].mean()
#         print(f"  {dataset['label']}: {overall_mean:.4f}")
    
#     # By region statistics
#     if 'region' in datasets[0]['df'].columns:
#         print("\nMean conservation by region:")
#         for region_name in ["5'_intron", 'first_exon_edge', 'last_exon_edge', "3'_intron"]:
#             print(f"\n  {region_name}:")
#             for dataset in datasets:
#                 region_data = dataset['df'][dataset['df']['region'] == region_name]
#                 if len(region_data) > 0:
#                     mean_val = region_data[dataset['mean_col']].mean()
#                     print(f"    {dataset['label']}: {mean_val:.4f}")
    
#     # Differences
#     if len(datasets) == 2:
#         print("\n" + "="*70)
#         print("PAIRWISE DIFFERENCE")
#         print("="*70)
        
#         df1 = datasets[0]['df']
#         df2 = datasets[1]['df']
#         mean_col1 = datasets[0]['mean_col']
#         mean_col2 = datasets[1]['mean_col']
        
#         if len(df1) == len(df2):
#             diff = df1[mean_col1] - df2[mean_col2]
#             print(f"\n{datasets[0]['label']} - {datasets[1]['label']}:")
#             print(f"  Mean difference: {diff.mean():.4f}")
#             print(f"  Max difference: {diff.max():.4f} at position {diff.idxmax()}")
#             print(f"  Min difference: {diff.min():.4f} at position {diff.idxmin()}")
            
#             # Statistical test
#             from scipy import stats as scipy_stats
#             t_stat, p_value = scipy_stats.ttest_rel(df1[mean_col1], df2[mean_col2])
#             print(f"\n  Paired t-test:")
#             print(f"    t-statistic: {t_stat:.4f}")
#             print(f"    p-value: {p_value:.4e}")
#             print(f"    Significant: {'Yes' if p_value < 0.05 else 'No'} (α=0.05)")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description='Compare conservation profiles across multiple datasets (PhastCons or PhyloP)',
#         epilog="""
# Example usage:
#   # Auto-detect metric, save as PNG (default)
#   python plot_exon_intron.py \\
#       -i output/constitutive_summary.tsv output/alternative_summary.tsv \\
#       -l "Constitutive" "Alternative" \\
#       -o output/plots/comparison
  
#   # Save as PDF
#   python plot_exon_intron.py \\
#       -i file1.tsv file2.tsv file3.tsv \\
#       -l "Dataset1" "Dataset2" "Dataset3" \\
#       -o output/plots/comparison \\
#       -f pdf
  
#   # Explicitly specify PhyloP, save as SVG
#   python plot_exon_intron.py \\
#       -i file1.tsv file2.tsv \\
#       -l "Dataset1" "Dataset2" \\
#       -o output/plots/comparison \\
#       -m phylop -f svg
#         """
#     )
#     parser.add_argument('-i', '--input', nargs='+', required=True,
#                        help='Input _summary.tsv files (space-separated)')
#     parser.add_argument('-l', '--labels', nargs='+', required=True,
#                        help='Labels for each dataset (space-separated, same order as inputs)')
#     parser.add_argument('-o', '--output', required=True,
#                        help='Output prefix for plot files')
#     parser.add_argument('-m', '--metric', default='auto', choices=['auto', 'phastcons', 'phylop'],
#                        help='Conservation metric to plot (default: auto-detect)')
#     parser.add_argument('-f', '--format', default='png', choices=['png', 'pdf', 'svg'],
#                        help='Output file format (default: png)')
    
#     args = parser.parse_args()
    
#     print("="*70)
#     print("Conservation Profile Comparison")
#     print("="*70)
#     print(f"Number of datasets: {len(args.input)}")
#     for i, (file, label) in enumerate(zip(args.input, args.labels)):
#         print(f"  {i+1}. {label}: {file}")
#     print(f"Output prefix: {args.output}")
#     print(f"Metric: {args.metric}")
#     print(f"Format: {args.format}")
#     print("="*70 + "\n")
    
#     plot_conservation_comparison(args.input, args.labels, args.output, args.metric, args.format)
    
#     print("\n" + "="*70)
#     print("COMPARISON COMPLETE")
#     print("="*70)