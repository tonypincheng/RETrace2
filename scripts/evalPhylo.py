#!/usr/bin/env python3
import argparse
from tqdm import tqdm
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import pandas as pd
import matplotlib
import seaborn as sns
matplotlib.use('Agg')
import multiprocessing
import random
import pickle
manager = multiprocessing.Manager()

def multi_calcTriplets(triplet_group, NJTree, rootDist_pd, sampleDict, tripletDict):
    for triplet in tqdm(triplet_group):
        tree_dist = []
        ref_dist = []
        for sample_pair in ([0,1], [0,2], [1,2]):
            MRCA = NJTree.get_common_ancestor([triplet[sample_pair[0]], triplet[sample_pair[1]]])
            MRCA_dist = NJTree.get_distance(MRCA)
            tree_dist.append(MRCA_dist)
            ref_dist.append(rootDist_pd.loc[sampleDict[triplet[sample_pair[0]]]["clone"]][sampleDict[triplet[sample_pair[1]]]["clone"]])
        #We want to save the count information in encoded string treeDict[triplet]: ref_dist,[1 (correct) or 0 (incorrect)]
        if tree_dist.index(max(tree_dist)) == ref_dist.index(max(ref_dist)): #If using MRCA depth as comparison metric (find most related pair)
            tripletDict[triplet] = str(max(ref_dist)) + ",1"
        else:
            tripletDict[triplet] = str(max(ref_dist)) + ",0"
    return

def calc_tripletAccuracy(tripletDict, prefix):
    errorDict = {}
    for triplet in tqdm(sorted(tripletDict.keys())):
        (ref_dist, correct_bool) = tripletDict[triplet].split(',')
        if int(ref_dist) != 0: #We want to ignore ref_dist of 0, meaning that either pair of triplet could potentially be closest to MRCA
            if ref_dist not in errorDict.keys():
                errorDict[ref_dist] = {}
                errorDict[ref_dist]["Correct"] = 0
                errorDict[ref_dist]["Total"] = 0
            errorDict[ref_dist]["Correct"] += int(correct_bool)
            errorDict[ref_dist]["Total"] += 1
    #Plot propotion of correct triplets per distance metric
    dist_list = sorted(errorDict.keys())
    corr_list = [] #Number of correct triplets
    rate_corr_list = [] #Rate of correct triplets
    num_list = [] #Total number of triplets analyzed
    [total_corr, total_triplets] = [0,0]
    for dist in dist_list:
        corr_list.append(errorDict[dist]["Correct"])
        num_list.append(errorDict[dist]["Total"])
        rate_corr_list.append(errorDict[dist]["Correct"] / errorDict[dist]["Total"])
        total_corr += errorDict[dist]["Correct"]
        total_triplets += errorDict[dist]["Total"]
    dist_list.append("All Triplets")
    corr_list.append(total_corr)
    rate_corr_list.append(total_corr / total_triplets)
    num_list.append(total_triplets)
    corr_df = pd.DataFrame({"Distance": dist_list, "Correct Triplets": rate_corr_list})
    sns_barplot = sns.barplot(x="Distance", y="Correct Triplets", data=corr_df)
    fig = sns_barplot.get_figure()
    fig.savefig(prefix + ".evalPhylo.eps")

    #Write correct triplet rate into text file
    f_output = open(prefix + ".evalPhylo.txt", 'w')
    f_output.write("\n".join([str(dist_list[i]) + "\t" + str(corr_list[i]) + "\t" + str(num_list[i]) + "\t" + str(float(corr_list[i] / num_list[i])) for i in range(len(dist_list))]) + "\n")
    f_output.close()

    return

def load_samplesheet(samplesheet_path):
    """
    Read the samplesheet CSV file and extract sample IDs and group information.
    Expects columns named 'sample_id' and 'group'.
    Returns a dictionary similar to the old sampleDict format.
    """
    sampleDict = {}
    try:
        df = pd.read_csv(samplesheet_path)
        
        # Check if required columns exist
        if 'sample_id' not in df.columns:
            print("Error: 'sample_id' column not found in samplesheet.")
            return {}
        
        # Use 'group' column if available, otherwise use 'unknown'
        for _, row in df.iterrows():
            sample_id = row['sample_id']
            group = row['group'] if 'group' in row and not pd.isna(row['group']) else 'unknown'
            
            sampleDict[sample_id] = {
                'group': group,  # Using 'group' instead of 'clone'
                'clone': group   # Keep 'clone' for backward compatibility
            }
        
        if not sampleDict:
            print("Warning: No sample mappings found in samplesheet.")
        return sampleDict
        
    except FileNotFoundError:
        print(f"Error: File '{samplesheet_path}' not found.")
        return {}
    except Exception as e:
        print(f"Error reading samplesheet: {e}")
        return {}

def evalPhylo(samplesheet, prefix, exVivo_rootDist, tree_file, nproc, distDict_file):
    '''
    This script is from the original RETrace repository, specifically for our known ex vivo cell culture tree.  It will be used to determine the accuracy of any phylogenetic tree we calculate.  To do this, we need to input the following files:
        1) samplesheet = CSV file containing sample information with 'sample_id' and 'group' columns
        2) exVivo_dist = csv file containing MRCA distance from root, as approximated in units of cell divisions
        3) newick_tree = file containing Newick tree output
        4) prefix = output prefix for error calculation statistics comparing calculated Newick tree to given ex vivo tree
    We will then use the above input to calculate an errorDict which contains the following structure:
        errorDict
            cell_div = reference cell division difference between nodes (ex: [2-1-G10_3-1-A2, 2-1-G10_3-1-B1, 2-2-B1_3-2-A6] = abs(max())) 
                "Correct" = number of triplets correct
                "Total" = total number of triplets analyzed
    '''
    sampleDict = load_samplesheet(samplesheet)
    #We want to extract only the group information for sampleDict
    distDict = pickle.load(open(distDict_file, 'rb'))

    #Import exVivo tree distances into pandas dataframe
    rootDist_pd = pd.read_csv(exVivo_rootDist, delimiter=',', index_col=0)

    #Import Newick tree and calculate all triplets of leaves in tree [sample] (with at least two clones per triple).  Determine whether distances between leaves is correct based on rootDist_pd
    with open(tree_file, 'r') as f:
        newick_tree = f.read().replace("\n",'')
    NJTree = Tree(newick_tree)

    #Create list of all triplets and run through to determine whether each triplet is correct
    tripletDict = manager.dict() #This allows for parallel processing of triplet errors for multiple triplets at once
    jobs = []
    triplet_set = set()
    print("Naming all triplets in tree")
    for sample1 in tqdm(distDict["samples"]):
        clone1 = sampleDict[sample1]["clone"]
        for sample2 in sorted(distDict["samples"]):
            clone2 = sampleDict[sample2]["clone"]
            for sample3 in sorted(distDict["samples"]):
                clone3 = sampleDict[sample3]["clone"]
                triplet = tuple(sorted([sample1, sample2, sample3]))
                if len(set([clone1, clone2, clone3])) >= 2 and len(set(triplet)) == 3:
                # if len(set([clone1, clone2, clone3])) == 3 and triplet not in triplet_list: #We want only triplets where each of the three leaves stem from different clones
                    triplet_set.add(triplet)
    triplet_list = list(triplet_set)

    print("Calculating correct triplet rate")
    random.shuffle(triplet_list) #We want to randomize in order to even processing time
    for triplet_group in [triplet_list[i::nproc] for i in range(nproc)]:
        p = multiprocessing.Process(target = multi_calcTriplets, args = (triplet_group, NJTree, rootDist_pd, sampleDict, tripletDict))
        jobs.append(p)
        p.start()
    #Join tripletDict
    for p in jobs:
        p.join()

    #Plot percentage of correct triplets
    print("Plotting/printing correct triplet rate")
    calc_tripletAccuracy(tripletDict, prefix)

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate accuracy of phylogeny given exVivo tree data"
    )
    
    parser.add_argument("--samplesheet", required=True,
        help="CSV file with sample information (sample_id, group columns)")
    parser.add_argument("--tree", required=True,
        help="Newick tree file")
    parser.add_argument("--prefix", required=True,
        help="Output prefix for results")
    parser.add_argument("--distDict", required=True,
        help="Pickle file containing distance dictionary")
    parser.add_argument("--exVivo_rootDist", 
        default="~/software/RETrace/Data/exVivo.rootDist.csv",
        help="CSV file with MRCA distances from root")
    parser.add_argument("--nproc", type=int, default=10,
        help="Number of processors to use")
    
    args = parser.parse_args()
    
    # Call the evalPhylo function with parsed arguments
    evalPhylo(
        samplesheet=args.samplesheet,
        prefix=args.prefix,
        exVivo_rootDist=args.exVivo_rootDist,
        tree_file=args.tree,
        nproc=args.nproc,
        distDict_file=args.distDict
    )

if __name__ == "__main__":
    main()
