#!/usr/bin/env python3
"""
build_phylo.py - Phylogenetic Tree Building Script with bootstrapping for microsatellite targets
"""

import pickle
import more_itertools
import skbio
import ete3
import numpy as np
import statistics
import argparse
import re
import random

def calcBulk(alleleDict):
    """
    Groups allele sizes into consecutive repeat units for phylogenetic analysis.

    Converts allele sizes (in base pairs) to repeat counts by dividing by motif length,
    collects unique repeats across samples, groups consecutive repeat counts, and stores
    them as base pair sizes in 'allele_groups' for each target ID.
    """
    for target_id in sorted(list(alleleDict.keys())):
        motif_len = alleleDict[target_id]["motif_len"]
        all_alleles = set()
        for sample in sorted(alleleDict[target_id]["sample"].keys()):
            try:
                all_alleles.update([int(allele / motif_len) for allele in alleleDict[target_id]["sample"][sample]["allelotype"]])
            except:
                continue
        alleleDict[target_id]["allele_groups"] = {}
        for idx, group in enumerate(more_itertools.consecutive_groups(sorted(all_alleles))):
            alleleDict[target_id]["allele_groups"][idx] = [allele * motif_len for allele in list(group)]
    return alleleDict

def calcDist(alleleDict, distDict, sample_pair, sample1, sample2, shared_targets, dist_metric, sample_list):
    """
    Calculate genetic distance between two samples based on their microsatellite profiles.
            
    Distance Metrics:
    ----------------
    - minComp: Compares only the minimum number of alleles between samples (assumes allele dropout and choose the minimum possible distance)
    - allComp: Compares all possible pairs of alleles between samples (assumes both alleles are present)
    - EqorNot: Binary distance - counts alleles as either matching (0) or different (1)
    - Abs: Continuous distance - uses the absolute difference between allele sizes
    """
    if sample1 != sample2:
        total_dist = 0
        num_comp = 0
        for target_id in shared_targets:
            for allele_idx, allele_group in alleleDict[target_id]["allele_groups"].items():
                allelotype1 = list(set(alleleDict[target_id]["sample"][sample1]["allelotype"]).intersection(set(allele_group)))
                allelotype2 = list(set(alleleDict[target_id]["sample"][sample2]["allelotype"]).intersection(set(allele_group)))
                if "minComp" in dist_metric:
                    dist_list = []
                    n = min(len(allelotype1), len(allelotype2)) # number of alleles to compare
                    if n == 0:
                        continue
                    for i in range(0, len(allelotype1), n):
                        temp_allelotype1 = allelotype1[i:i + n]
                        for j in range(0, len(allelotype2), n):
                            temp_allelotype2 = allelotype2[j:j + n]
                            if "EqorNot" in dist_metric:
                                dist = int(len(set(temp_allelotype1).symmetric_difference(set(temp_allelotype2))) / 2)
                            elif "Abs" in dist_metric:
                                dist = sum(abs(x - y) for x, y in zip(sorted(temp_allelotype1), sorted(temp_allelotype2)))
                            dist_list.append(dist)
                    total_dist += min(dist_list)
                    num_comp += n
                elif "allComp" in dist_metric:
                    for allele1 in sorted(allelotype1):
                        for allele2 in sorted(allelotype2):
                            if "EqorNot" in dist_metric:
                                if allele1 != allele2:
                                    total_dist += 1
                            elif "Abs" in dist_metric:
                                total_dist += abs(allele1 - allele2)
                            num_comp += 1
        if sample_pair == tuple(sorted(sample_list[:2])):
            print(f"\tDistance for {sample_pair}: targets={len(shared_targets)}, total_dist={total_dist}, num_comp={num_comp}, dist={total_dist/num_comp if num_comp > 0 else 0}")
    else:
        total_dist = 0
        num_comp = len(shared_targets)
    distDict["sampleComp"][sample_pair] = {
        "dist": float(total_dist / num_comp) if num_comp > 0 else 0,
        "num_targets": len(shared_targets),
        "num_alleles": int(num_comp),
        "num_diff": int(total_dist)
    }
    return distDict

def makeDistMatrix(sharedDict, alleleDict, sample_list, dist_metric):
    distDict = {"samples": sample_list, "sampleComp": {}}
    seen = set()
    for sample1 in sorted(sample_list):
        for sample2 in sorted(sample_list):
            sample_pair = tuple(sorted([sample1, sample2]))
            if sample_pair not in seen:
                seen.add(sample_pair)
                distDict = calcDist(alleleDict, distDict, sample_pair, sample1, sample2, sharedDict[sample_pair], dist_metric, sample_list)
    return distDict

def clean_node_names(tree):
    """Remove quotes from node names in an ETE3 tree.
    """
    for node in tree.traverse():
        if node.name and node.name.startswith("'") and node.name.endswith("'"):
            node.name = node.name[1:-1]
    return tree

def drawTree(distDict, sample_list, outgroup, prefix, bootstrap):
    print(f"\tDrawing tree: Building distance matrix for {len(sample_list)} samples...")
    distMatrix = [[distDict["sampleComp"][tuple(sorted([s1, s2]))]["dist"] 
                   for s2 in sorted(sample_list)] for s1 in sorted(sample_list)]
    if not bootstrap:
        targetMatrix = [[distDict["sampleComp"][tuple(sorted([s1, s2]))]["num_targets"] 
                         for s2 in sorted(sample_list)] for s1 in sorted(sample_list)]
        pairwise_numTargets = [distDict["sampleComp"][tuple(sorted([s1, s2]))]["num_targets"] 
                               for s1 in sorted(sample_list) for s2 in sorted(sample_list) if s1 != s2]
        sample_numTargets = [distDict["sampleComp"][(s, s)]["num_targets"] for s in sorted(sample_list)]
        with open(prefix + ".buildPhylo.stats.txt", 'w') as statsOutput:
            statsOutput.write(f"Number of Samples Analyzed:\t{len(sample_list)}\n{','.join(sample_list)}\n")
            avg_targets = float(sum(pairwise_numTargets) / len(pairwise_numTargets)) if pairwise_numTargets else 0
            avg_sample_targets = float(sum(sample_numTargets) / len(sample_numTargets)) if sample_numTargets else 0
            statsOutput.write(f"Avg targets captured per single cell:\t{avg_sample_targets}\t[Min:{min(sample_numTargets)}, Max:{max(sample_numTargets)}]\n")
            statsOutput.write(f"Avg targets shared per pair of cells:\t{avg_targets}\t[Min:{min(pairwise_numTargets)}, Max:{max(pairwise_numTargets)}']\n")
            
            statsOutput.write("\nDistance Matrix (pairwise distances between samples):\n")
            statsOutput.write("Sample," + ",".join(sorted(sample_list)) + "\n")
            for dist_idx, dist_list in enumerate(distMatrix):
                statsOutput.write(f"{sorted(sample_list)[dist_idx]},{','.join(str(round(i, 4)) for i in dist_list)}\n") # change to 4 decimal places from 3
            
            statsOutput.write("\nTarget Matrix (number of shared targets between samples):\n")
            statsOutput.write("Sample," + ",".join(sorted(sample_list)) + "\n")
            for target_idx, target_list in enumerate(targetMatrix):
                statsOutput.write(f"{sorted(sample_list)[target_idx]},{','.join(str(j) for j in target_list)}\n")
        pickle.dump(distDict, open(prefix + ".buildPhylo.distDict.pkl", "wb"))

    print(f"\tDrawing tree: Creating DistanceMatrix...")
    distObj = skbio.DistanceMatrix(distMatrix, sorted(sample_list))
    
    print(f"\tDrawing tree: Running NJ algorithm...")
    try:
        skbio_tree = skbio.tree.nj(distObj)
        # Convert skbio TreeNode to newick string
        newick_str = str(skbio_tree)
        print(f"\tDrawing tree: Converting to ETE tree...")
        ete_tree = ete3.Tree(newick_str)
        # Clean node names
        ete_tree = clean_node_names(ete_tree)
        
    except Exception as e:
        print(f"\tNJ failed: {e}")
        raise
    
    if outgroup != "NA":
        print(f"\tDrawing tree: Setting outgroup {outgroup}...")
        if outgroup == "Midpoint":
            tree_midpoint = ete_tree.get_midpoint_outgroup()
            if tree_midpoint is not None:
                ete_tree.set_outgroup(tree_midpoint)
            else:
                print("Warning: Midpoint outgroup not found. Proceeding without outgroup.")
        else:
            ete_tree.set_outgroup(outgroup)
    print(f"\tDrawing tree: Completed.")
    return ete_tree

def bootstrap_iteration(args):
    i, target_list, sharedDict, alleleDict, sample_list, dist_metric, outgroup, prefix = args
    print(f"\tIteration {i}: Sampling {len(target_list)} targets...")
    sample_size = int(0.8 * len(target_list))  # specify sample 80% of targets
    bootstrap_targets = set(np.random.choice(target_list, sample_size, replace=False))
    unique_targets = len(bootstrap_targets)
    print(f"\tIteration {i}: Sampled {unique_targets} unique targets ({unique_targets/len(target_list)*100:.2f}%)")
    print(f"\tIteration {i}: Building bootstrap_sharedDict for {len(sharedDict)} pairs...")
    bootstrap_sharedDict = {sample_pair: [t for t in sharedDict[sample_pair] if t in bootstrap_targets] 
                            for sample_pair in sharedDict}
    
    distDict_temp = makeDistMatrix(bootstrap_sharedDict, alleleDict, sample_list, dist_metric)
    print(f"\tIteration {i}: Drawing tree...")
    tree_temp = drawTree(distDict_temp, sample_list, outgroup, prefix, True)
    # Note: No need to clean node names here as drawTree already does it
    return tree_temp

def bootstrapTree(nodeDict, treeTemp, bootstrap_samples):
    tempNodes = set()
    for node in treeTemp.search_nodes():
        leaf_list = [leaf.name for leaf in node]
        tempNodes.add(tuple(sorted(set(leaf_list))))
    for node in nodeDict.keys():
        node_intersect = tuple(sorted(set(node).intersection(bootstrap_samples)))
        if node_intersect == node:
            nodeDict[node]["Num_sampled"] += 1
            if node_intersect in tempNodes:
                nodeDict[node]["Num_verified"] += 1
    return nodeDict

def buildPhylo(sample_list_file, prefix, alleleDict_file, dist_metric, outgroup, bootstrap, bootstrap_iterations=100):
    print("Building phylogenetic tree...")
    print("NOTE: Tree topology may vary between runs due to tie-breaking in neighbor-joining algorithm")

    print("Loading allele dictionary...")
    with open(alleleDict_file, 'rb') as f:
        alleleDict = pickle.load(f)
    
    print("Processing bulk alleles...")
    alleleDict = calcBulk(alleleDict)
    
    # Get samples from the sample_list file
    print("Reading sample list...")
    with open(sample_list_file, 'r') as f:
        sample_list = f.read().splitlines()
    print(f"Found {len(sample_list)} samples in sample list")

    print("\tPre-computing shared target_id between each pairwise sample")
    sharedDict = {}
    for sample1 in sample_list:
        for sample2 in sorted(sample_list):
            sample_pair = tuple(sorted([sample1, sample2]))
            if sample_pair not in sharedDict:
                shared_targets = [
                    target_id for target_id in alleleDict
                    if sample1 in alleleDict[target_id]["sample"] and
                    sample2 in alleleDict[target_id]["sample"] and
                    "allelotype" in alleleDict[target_id]["sample"][sample1] and
                    "allelotype" in alleleDict[target_id]["sample"][sample2]
                ]
                sharedDict[sample_pair] = shared_targets

    print("\tCalculating distance matrix and building original NJ tree")
    distDict_original = makeDistMatrix(sharedDict, alleleDict, sample_list, dist_metric)
    
    print("Drawing tree...")
    MStree = drawTree(distDict_original, sample_list, outgroup, prefix, False)

    print("Saving newick tree...")
    with open(prefix + '.buildPhylo.newick-original.txt', 'w') as f:
        f.write(MStree.write(format=0) + "\n")
    print("Saving distance dictionary...")
    with open(prefix + ".buildPhylo.distDict.pkl", "wb") as f:
        pickle.dump(distDict_original, f)

    if bootstrap:
        print("Setting up bootstrap...")
        nodeDict = {tuple(sorted([leaf.name for leaf in node])): {"Num_verified": 0, "Num_sampled": 0}
                    for node in MStree.search_nodes() if not node.is_leaf()}
        target_list = list(alleleDict.keys())
        print(f"\tStarting bootstrap with {len(target_list)} targets, {bootstrap_iterations} iterations...")
        for i in range(bootstrap_iterations):
            print(f"\tRunning iteration {i+1}/{bootstrap_iterations}...")
            tree_temp = bootstrap_iteration((i, target_list, sharedDict, alleleDict, sample_list, dist_metric, outgroup, prefix))
            if tree_temp is not None:
                nodeDict = bootstrapTree(nodeDict, tree_temp, sample_list)

        print(f"\tUpdating node support values...")
        bootstrap_values = []
        for node in MStree.search_nodes():
            if not node.is_leaf():
                leaf_list = tuple(sorted([leaf.name for leaf in node]))
                if nodeDict[leaf_list]["Num_sampled"] > 0:
                    node_support = round(float(nodeDict[leaf_list]["Num_verified"] / nodeDict[leaf_list]["Num_sampled"]), 2)
                    bootstrap_values.append(node_support)  # Collect valid support values
                else:
                    node_support = 99.0  # Placeholder for untested nodes
                node.add_features(support=node_support)

        print(f"\tWriting bootstrap newick file...")
        with open(prefix + '.buildPhylo.newick-bootstrap.txt', 'w') as f:
            f.write(MStree.write(format=0) + "\n")

        # Calculate and write bootstrap statistics
        if bootstrap_values:
            average_bootstrap = statistics.mean(bootstrap_values)
            median_bootstrap = statistics.median(bootstrap_values)
            with open(prefix + "_bootstrap_stats.txt", 'w') as f:
                f.write(f"Average Bootstrap Value: {average_bootstrap:.2f}\n")
                f.write(f"Median Bootstrap Value: {median_bootstrap:.2f}\n")
                f.write(f"Number of nodes with support values: {len(bootstrap_values)}\n")
                f.write(f"Number of nodes without support (value = 99.0): {sum(1 for node in MStree.traverse() if not node.is_leaf() and node.support == 99.0)}\n")
                f.write("\nAll Bootstrap Values (0-1 only):\n")
                for i, value in enumerate(bootstrap_values, 1):
                    f.write(f"{i}. {value:.2f}\n")
        else:
            print("No bootstrap values found to calculate statistics.")

def main():
    parser = argparse.ArgumentParser(description='Build phylogenetic tree from HipSTR allele data')
    parser.add_argument('--sample_list', required=True, help='Input file with sample list')
    parser.add_argument('--alleleDict', required=True, help='Input pickle file with allele dictionary')
    parser.add_argument('--prefix', required=True, help='Output file prefix')
    parser.add_argument('--dist_metric', default='minComp_EqorNot', help='Distance metric (default: minComp_EqorNot)')
    parser.add_argument('--outgroup', default='Midpoint', help='Outgroup for tree rooting (default: Midpoint)')

    parser.add_argument('--bootstrap', action='store_true', help='Run bootstrap analysis')
    parser.add_argument('--bootstrap_iterations', type=int, default=100, help='Number of bootstrap iterations (default: 100)')
    
    args = parser.parse_args()
    
    buildPhylo(
        args.sample_list,
        args.prefix,
        args.alleleDict,
        args.dist_metric,
        args.outgroup,
        args.bootstrap,
        args.bootstrap_iterations
    )

if __name__ == "__main__":
    main() 