#!/usr/bin/env python3
import os
import argparse
import pandas as pd
os.environ['QT_QPA_PLATFORM']='offscreen' #Fix for remote ssh tree render <https://github.com/etetoolkit/ete/issues/387>
import ete3  

def load_sample_data(samplesheet_path):
    """Load sample information from CSV file"""
    df = pd.read_csv(samplesheet_path)
    samples = {}
    
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        sample_info = {}
        
        # Handle optional group field
        sample_info['group'] = row['group'] if 'group' in row and not pd.isna(row['group']) else 'unknown'
        
        # Handle optional color field with default grey
        sample_info['color'] = row['color'] if 'color' in row and not pd.isna(row['color']) else 'grey'
        
        samples[sample_id] = sample_info
    
    return samples

def setup_tree_style(show_bootstrap=False, circular=False):
    """Configure tree style settings"""
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    #ts.show_branch_support = show_bootstrap
    ts.show_scale = False  # Disable the scale bar in corner
    ts.show_branch_support = False # set this to False becasue we are customizing the support text (see below)
    ts.mode = 'c' if circular else 'r'  # 'c' for circular, 'r' for rectangular
    #ts.branch_vertical_margin = 5
    
    return ts

def style_tree_nodes(tree, sample_data, color_background=False, bootstrap=False):
    """Apply styling to tree nodes based on sample data
    
    Args:
        tree: The phylogenetic tree
        sample_data: Dictionary of sample information
        color_background: If True, applies color to background instead of node
    """
    # Set larger text size for better readability
    text_size = 18
    support_text_size = 18
    
    for node in tree.traverse():
        nstyle = ete3.NodeStyle()
        
        nstyle["hz_line_width"] = 2  # Thickness for horizontal lines
        nstyle["vt_line_width"] = 2  # Thickness for vertical lines
        
        if node.is_leaf():
            # Skip nodes not in sample data (like outgroups)
            if node.name not in sample_data:
                print(f"Node {node.name} not in sample data")
                print(f"Sample data keys: {sample_data.keys()}")
                continue
                
            color = sample_data[node.name]['color']
            
            if color_background:
                nstyle["fgcolor"] = 'black'
                nstyle["size"] = 15  # Smaller node
                nstyle["bgcolor"] = color
                node.set_style(nstyle)
                
                # Create text with colored background
                name_face = ete3.TextFace(' ' + node.name + ' ', fgcolor='black', fsize=text_size)
                node.add_face(name_face, column=0, position='branch-right')
            else:
                # Default: apply color to node circle
                nstyle["fgcolor"] = color
                nstyle["size"] = 15
                nstyle["bgcolor"] = 'white'
                node.set_style(nstyle)
                
                # Plain text without background color
                name_face = ete3.TextFace(' ' + node.name, fgcolor='black', fsize=text_size)
                node.add_face(name_face, column=0, position='branch-right')
            
        else:
            nstyle["size"] = 0
            
            # Customize support text to be larger
            if bootstrap and node.support is not None:  # Check if support value exists
                support_text = f"{node.support:.2f}  "  # add space after support text for better visualization
                support_face = ete3.TextFace(support_text, fgcolor='#8B0000', fsize=support_text_size)
                node.add_face(support_face, column=2, position='branch-bottom')  
            node.set_style(nstyle)
            
    return tree

def render_tree(tree, prefix, tree_style):
    """Render tree to PDF and PNG files"""
    tree.render(f"{prefix}.viewPhylo.pdf", tree_style=tree_style, dpi=500)
    tree.render(f"{prefix}.viewPhylo.png", tree_style=tree_style, dpi=500)
    

def viewPhylo(samplesheet, tree_file, prefix, bootstrap, color_background=False, circular=False):
    """
    Visualize phylogenetic tree with colored nodes based on sample information
    """
    # Load data
    sample_data = load_sample_data(samplesheet)
    tree = ete3.Tree(tree_file)
    
    # Setup and apply styling
    tree_style = setup_tree_style(bootstrap, circular)
    styled_tree = style_tree_nodes(tree, sample_data, color_background, bootstrap)
    
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
    parser.add_argument('--circular', action='store_true', 
                       help='Draw tree in circular layout instead of rectangular')
    
    args = parser.parse_args()
    viewPhylo(args.samplesheet, args.tree_file, args.prefix, args.bootstrap, 
             args.color_background, args.circular)

if __name__ == '__main__':
    main() 