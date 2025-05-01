#!/usr/bin/env python3
import os
import argparse
import pandas as pd
os.environ['QT_QPA_PLATFORM']='offscreen' #Fix for remote ssh tree render <https://github.com/etetoolkit/ete/issues/387>
import ete3  # Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>

def load_sample_data(samplesheet_path):
    """Load sample information from CSV file"""
    df = pd.read_csv(samplesheet_path)
    samples = {row['sample_id']: {'group': row['group'], 'color': row['color']} 
              for _, row in df.iterrows()}
    return samples

def setup_tree_style(show_bootstrap=False):
    """Configure tree style settings"""
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = show_bootstrap
    return ts

def style_tree_nodes(tree, sample_data, color_background=False):
    """Apply styling to tree nodes based on sample data
    
    Args:
        tree: The phylogenetic tree
        sample_data: Dictionary of sample information
        color_background: If True, applies color to background instead of node
    """
    for node in tree.traverse():
        nstyle = ete3.NodeStyle()
        
        if node.is_leaf():
            # Skip nodes not in sample data (like outgroups)
            if node.name not in sample_data:
                continue
                
            color = sample_data[node.name]['color']
            
            if color_background:
                nstyle["fgcolor"] = 'black'
                nstyle["size"] = 15  # Smaller node
                nstyle["bgcolor"] = color
                node.set_style(nstyle)
                
                # Create text with colored background
                name_face = ete3.TextFace(' ' + node.name + ' ', fgcolor='black', fsize=10)
                node.add_face(name_face, column=0, position='branch-right')
            else:
                # Default: apply color to node circle
                nstyle["fgcolor"] = color
                nstyle["size"] = 15
                nstyle["bgcolor"] = 'white'
                node.set_style(nstyle)
                
                # Plain text without background color
                name_face = ete3.TextFace(' ' + node.name, fgcolor='black', fsize=10)
                node.add_face(name_face, column=0, position='branch-right')
        else:
            nstyle["size"] = 0
            node.set_style(nstyle)
            
    return tree

def render_tree(tree, prefix, tree_style):
    """Render tree to PDF and PNG files"""
    tree.render(f"{prefix}.viewPhylo.pdf", tree_style=tree_style, dpi=180)
    tree.render(f"{prefix}.viewPhylo.png", tree_style=tree_style, dpi=180)

def viewPhylo(samplesheet, tree_file, prefix, bootstrap, color_background=False):
    """
    Visualize phylogenetic tree with colored nodes based on sample information
    """
    # Load data
    sample_data = load_sample_data(samplesheet)
    tree = ete3.Tree(tree_file)
    
    # Setup and apply styling
    tree_style = setup_tree_style(bootstrap)
    styled_tree = style_tree_nodes(tree, sample_data, color_background)
    
    # Render output files
    render_tree(styled_tree, prefix, tree_style)

def main():
    parser = argparse.ArgumentParser(description='Visualize phylogenetic tree using ete3')
    parser.add_argument('--samplesheet', required=True, help='Path to sample information CSV file')
    parser.add_argument('--tree_file', required=True, help='Path to newick tree file')
    parser.add_argument('--prefix', required=True, help='Output prefix for visualization files')
    parser.add_argument('--bootstrap', action='store_true', help='Show bootstrap values')
    parser.add_argument('--color_background', action='store_true', 
                       help='Apply colors to text background instead of node circles')
    
    args = parser.parse_args()
    viewPhylo(args.samplesheet, args.tree_file, args.prefix, args.bootstrap, args.color_background)

if __name__ == '__main__':
    main() 