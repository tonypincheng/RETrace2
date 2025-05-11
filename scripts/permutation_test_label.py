from Bio import Phylo
import numpy as np
import random
from collections import Counter
import matplotlib.pyplot as plt
import argparse
import os
import sys
import time
from tqdm import tqdm  # Import tqdm for progress bars

# Set a random seed for reproducibility
random.seed(42)  # You can change 42 to any integer for a different but reproducible sequence

def read_tissue_mapping(file_path):
    """
    Read the samplesheet file and extract sample IDs and group information.
    Expects columns named 'sample_id' and 'group'.
    
    Parameters:
    - file_path: Path to the samplesheet file
    
    Returns:
    - group_dc: Dictionary mapping sample IDs to group/tissue labels
    """
    group_dc = {}
    try:
        with open(file_path, 'r') as file:
            header = None
            sample_idx = None
            group_idx = None
            
            for line in file:
                # Skip comments
                if line.startswith('#'):
                    continue
                    
                columns = line.strip().split(',')  # Assuming CSV format
                
                # Process header
                if header is None:
                    header = columns
                    try:
                        sample_idx = header.index('sample_id')
                        group_idx = header.index('group')
                    except ValueError:
                        print("Error: Could not find 'sample_id' and 'group' columns in samplesheet.")
                        print(f"Header found: {header}")
                        return {}
                    continue
                
                # Process data rows
                if len(columns) > max(sample_idx, group_idx):
                    sample_id = columns[sample_idx].strip()
                    group = columns[group_idx].strip()
                    if sample_id and group:
                        group_dc[sample_id] = group
        
        if not group_dc:
            print("Warning: No group mappings found in samplesheet.")
        return group_dc
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return {}


def calculate_cophenetic_distance(tree, tip1_name, tip2_name):
    """
    Calculate the cophenetic distance between two tips in a phylogenetic tree.
    Returns the number of edges (unweighted) between the tips via their MRCA.
    
    Parameters:
    - tree: Bio.Phylo tree object
    - tip1_name: Name of the first tip
    - tip2_name: Name of the second tip
    
    Returns:
    - distance: Number of edges (cophenetic distance) or -1 if tips not found
    """
    # Find the tips in the tree
    tip1 = [t for t in tree.get_terminals() if t.name == tip1_name]
    tip2 = [t for t in tree.get_terminals() if t.name == tip2_name]
    
    if not tip1 or not tip2:
        return -1  # Tips not found
    
    tip1 = tip1[0]
    tip2 = tip2[0]
    
    # Get the path from tip1 to root
    path1 = []
    current = tip1
    while current:
        path1.append(current)
        if current == tree.root:
            break
        # Find the parent by checking all clades in the tree
        found = False
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                found = True
                break
        if not found:
            break
    
    # Get the path from tip2 to root
    path2 = []
    current = tip2
    while current:
        path2.append(current)
        if current == tree.root:
            break
        # Find the parent by checking all clades in the tree
        found = False
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                found = True
                break
        if not found:
            break
    
    # Find the MRCA by walking up from both tips until paths intersect
    mrca = None
    for node1 in path1:
        for node2 in path2:
            if node1 == node2:
                mrca = node1
                break
        if mrca:
            break
    
    if not mrca:
        return -1  # No MRCA found (error case)
    
    # Count edges from tip1 to MRCA
    dist1 = 0
    current = tip1
    while current != mrca:
        dist1 += 1
        # Find the parent
        found = False
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                found = True
                break
        if not found or not current:
            return -1  # Shouldn't happen if MRCA is found
    
    # Count edges from tip2 to MRCA
    dist2 = 0
    current = tip2
    while current != mrca:
        dist2 += 1
        # Find the parent
        found = False
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                found = True
                break
        if not found or not current:
            return -1  # Shouldn't happen if MRCA is found
    
    return dist1 + dist2

def find_siblings(tree, tip):
    """
    Find sibling tips (other tips sharing the same parent) for a given tip.
    
    Parameters:
    - tree: Bio.Phylo tree object
    - tip: Bio.Phylo Clade object (a terminal node)
    
    Returns:
    - siblings: List of Bio.Phylo Clade objects (sibling tips)
    """
    if tip == tree.root:
        return []  # Root has no siblings
    
    # Find the parent by checking all clades
    parent = None
    for clade in tree.find_clades():
        if tip in clade.clades:
            parent = clade
            break
    
    if not parent:
        return []  # No parent found
    
    # Find all terminal siblings under the parent
    siblings = [sibling for sibling in parent.clades if sibling.is_terminal() and sibling != tip]
    return siblings

def count_same_tissue_pairs(tree, group_dc, distance_threshold=2):
    """
    Count the percentage of same-tissue pairs in a phylogenetic tree, considering only:
    - Tip-tip pairs (siblings) with the same tissue and cophenetic distance <= distance_threshold.
    
    Parameters:
    - tree: Bio.Phylo tree object (loaded from Newick format)
    - group_dc: Dictionary mapping sample IDs to group/tissue labels
    - distance_threshold: Maximum cophenetic distance for considering pairs (default=2)
    
    Returns:
    - percentage: Percentage of unique same-tissue pairs out of unique actual comparisons made
    - comparisons: Number of unique actual comparisons made (siblings evaluated once)
    - pairs: Number of unique same-tissue pairs
    """
    # Get all tip names
    tip_names = [tip.name for tip in tree.get_terminals()]
    
    unique_pairs = set()  # Use a set to store unique same-tissue pairs
    unique_comparisons = set()  # Use a set to store unique comparisons (siblings evaluated once)
    
    # Iterate over each tip
    for tip1_name in tip_names:
        tissue1 = group_dc.get(tip1_name, None)  # Use .get() to handle missing keys
        if tissue1 is None:
            print(f"Warning: No group label found for tip {tip1_name}, skipping.")
            continue
        
        # Find the tip object
        tip1 = [t for t in tree.get_terminals() if t.name == tip1_name][0]
        
        # Find siblings
        siblings = find_siblings(tree, tip1)
        
        # Check for tip-tip pairs (siblings)
        for sibling in siblings:
            tip2_name = sibling.name
            # Use sorted tuple to ensure uniqueness (e.g., ("A", "B") instead of ("B", "A"))
            comparison = tuple(sorted([tip1_name, tip2_name]))
            unique_comparisons.add(comparison)
            
            tissue2 = group_dc.get(tip2_name, None)
            if tissue2 is None:
                print(f"Warning: No group label found for tip {tip2_name}, skipping.")
                continue
            
            # Calculate cophenetic distance manually
            dist = calculate_cophenetic_distance(tree, tip1_name, tip2_name)
            if dist == -1 or dist > distance_threshold:
                continue
            
            if tissue1 == tissue2:
                # Use a sorted tuple to ensure uniqueness for same-tissue pairs
                pair = tuple(sorted([tip1_name, tip2_name]))
                unique_pairs.add(pair)
    
    pairs = len(unique_pairs)  # Number of unique same-tissue pairs
    comparisons = len(unique_comparisons)  # Number of unique actual comparisons made
    
    # Return percentage of unique same-tissue pairs out of unique actual comparisons made
    if comparisons > 0:
        percentage = (pairs / comparisons) * 100
    else:
        percentage = 0.0
    
    return percentage, comparisons, pairs

def parse_arguments():
    """
    Parse command-line arguments for the permutation test.
    """
    parser = argparse.ArgumentParser(description="Perform permutation test on phylogenetic tree to assess tissue clustering.")
    
    # Required arguments
    parser.add_argument('--samplesheet', required=True, help='Path to the samplesheet CSV file with sample IDs and group/tissue information.')
    parser.add_argument('--tree_file', required=True, help='Path to the phylogenetic tree file in Newick format.')
    
    # Optional arguments
    parser.add_argument('--distance_threshold', type=int, default=2, help='Maximum cophenetic distance for considering pairs (default=2).')
    parser.add_argument('--num_permutations', type=int, default=1000, help='Number of permutations to perform (default=1000).')
    parser.add_argument('--random_seed', type=int, default=42, help='Random seed for reproducibility (default=42).')
    parser.add_argument('--output_prefix', default='permutation_test', help='Prefix for output files (default="permutation_test").')
    
    return parser.parse_args()

# Example usage
if __name__ == "__main__":
    print("\n========== Permutation Test for Phylogenetic Tree Group Clustering ==========\n")
    
    # Parse command-line arguments
    args = parse_arguments()
    print(f"Input parameters:")
    print(f"  Samplesheet: {args.samplesheet}")
    print(f"  Tree file: {args.tree_file}")
    print(f"  Distance threshold: {args.distance_threshold}")
    print(f"  Number of permutations: {args.num_permutations}")
    print(f"  Random seed: {args.random_seed}")
    print(f"  Output prefix: {args.output_prefix}")
    
    # Set random seed
    random.seed(args.random_seed)
    
    # Create output file names
    stats_file = f"{args.output_prefix}_stats.txt"
    histogram_file = f"{args.output_prefix}_histogram.png"
    
    # Load group mapping from samplesheet
    print("Loading group mapping from samplesheet...")
    group_dc = read_tissue_mapping(args.samplesheet)
    
    if not group_dc:
        print("Unable to proceed without group mappings. Exiting.")
        sys.exit(1)
    
    # Load the phylogenetic tree
    print("Loading phylogenetic tree...")
    try:
        tree = Phylo.read(args.tree_file, "newick")
        print(f"Tree loaded with {len(tree.get_terminals())} tips.")
    except FileNotFoundError:
        print(f"Error: Tree file '{args.tree_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading tree file: {e}")
        sys.exit(1)
    
    # Print basic information
    group_counts = Counter(group_dc.values())
    print("\nGroup distribution:")
    for group, count in group_counts.items():
        print(f"  {group}: {count} samples")
    
    # Calculate the percentage of same-tissue pairs with the specified threshold
    print(f"\nCalculating observed percentage of same-group pairs (distance threshold = {args.distance_threshold})...")
    observed_percentage, comparisons, pairs = count_same_tissue_pairs(tree, group_dc, args.distance_threshold)
    print(f"Observed percentage: {observed_percentage:.2f}%")
    print(f"Unique comparisons made: {comparisons}")
    print(f"Unique same-group pairs: {pairs}")
    
    # For the permutation test, randomize tissue labels and repeat
    print(f"\nStarting permutation test with {args.num_permutations} permutations...")
    n_perm = args.num_permutations
    random_percentages = []
    tip_names = list(group_dc.keys())
    group_values = list(set(group_dc.values()))  # Unique group labels
    
    start_time = time.time()
    
    # In a typical scenario, we'll have fewer unique groups than tips
    if len(group_values) < len(tip_names):
        # Using random.choices (sampling with replacement) when we have fewer groups than tips
        for _ in tqdm(range(n_perm), desc="Permutation progress"):
            perm_tissues = dict(zip(tip_names, random.choices(group_values, k=len(tip_names))))
            percentage, _, _ = count_same_tissue_pairs(tree, perm_tissues, args.distance_threshold)
            random_percentages.append(percentage)
    else:
        # Using random.sample (sampling without replacement) in the unusual case where groups ≥ tips
        for _ in tqdm(range(n_perm), desc="Permutation progress"):
            perm_tissues = dict(zip(tip_names, random.sample(group_values, len(tip_names))))
            percentage, _, _ = count_same_tissue_pairs(tree, perm_tissues, args.distance_threshold)
            random_percentages.append(percentage)
    
    elapsed_time = time.time() - start_time
    print(f"Permutation test completed in {elapsed_time:.2f} seconds.")
    
    # Calculate p-value 
    p_value = np.mean(np.array(random_percentages) >= observed_percentage)
    
    # Print results
    print("\nPermutation Test Results:")
    print(f"Observed percentage of same-group pairs: {observed_percentage:.2f}%")
    print(f"P-value: {p_value}")
    print(f"Random percentages - Mean: {np.mean(random_percentages):.2f}%, Median: {np.median(random_percentages):.2f}%, Max: {np.max(random_percentages):.2f}%")
    
    # Save important statistics to a text file
    with open(stats_file, "w") as f:
        f.write("Permutation Test (Shuffle Group Labels) Statistics:\n")
        f.write(f"Group distribution: {Counter(group_dc.values())}\n")
        f.write(f"Number of permutations: {n_perm}\n")
        f.write(f"Observed Percentage of Same-Tissue Pairs: {observed_percentage:.2f}%\n")
        f.write(f"Unique Comparisons Made: {comparisons}\n")
        f.write(f"Unique Same-Tissue Pairs: {pairs}\n")
        f.write(f"P-value: {p_value}\n")
        f.write(f"Random Percentages - Mean: {np.mean(random_percentages):.2f}%, Median: {np.median(random_percentages):.2f}%, Max: {np.max(random_percentages):.2f}%\n")
    
    print(f"\nStatistics saved to '{stats_file}'")
    
    # Visualize the distribution and save the plot with improved binning
    plt.figure(figsize=(8, 5), dpi=300)

    bins = np.arange(0, 100 + 5, 5)  # Bins from 0 to 100 with 5% increments
    plt.hist(random_percentages, bins=bins, edgecolor='black')
    plt.axvline(x=observed_percentage, color='red', linestyle='dashed', linewidth=2)
    plt.title("Permutation Test: Shuffle Group Labels\n", fontsize=14)
    plt.xlabel("\nPercentage of Same-Group Pairs (%)")
    plt.ylabel("Frequency\n")
    plt.xlim(0, 100)
    plt.xticks(np.arange(0, 101, 10))

    # Remove right and top spines
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Save figure
    plt.savefig(histogram_file, dpi=300, bbox_inches='tight') 
    print(f"Histogram saved to '{histogram_file}'")
    
    # Show plot if display is available
    if os.environ.get('DISPLAY'):
        plt.show()
    
    plt.close()  # Close the figure to free memory
    
    print("\n========== Permutation Test Complete ==========\n")