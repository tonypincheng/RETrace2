#!/usr/bin/env python3
import argparse
import pickle
import random

def import_targetDict(probe_file):
    '''Import probe file into targetDict'''
    targetDict = {}
    with open(probe_file) as f:
        for line in f:
            #target_info = line.split(',')[0]
            target_info = line.split()[0].split(",")[0] # 2023/01 update by Tony, split "\t" and "," so both new and old IDT output work
            info_list = target_info.split('_')
            #Reformat target_id and save to targetDict
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]
            target_id = chrom + ":" + chromStart + "-" + chromEnd

            #Save relevant information to targetDict
            targetDict[target_id] = {}
            (targetDict[target_id]["chrom"], targetDict[target_id]["chromStart"], targetDict[target_id]["chromEnd"]) = (chrom, chromStart, chromEnd)
            (targetDict[target_id]["num_sub"], targetDict[target_id]["sub_seq"]) = info_list[-1].split('x')
            #Also save the expected up/down_seq around the microsatellite
            MS_frag = line.split()[1]
            targetDict[target_id]["MS_frag"] = MS_frag
            frag_seq = targetDict[target_id]["sub_seq"]*int(targetDict[target_id]["num_sub"])
            frag_split = MS_frag.split(frag_seq)
            #There are a few cases (i.e. chr6:55179847-55179875) in which microsatellite frag_seq is found more than once in reference fragment.  Thus, we need to choose the optimal up/down-seq based on which pair has the greater length
            (up_seq, down_seq) = ('', '')
            for i in range(len(frag_split) - 1):
                if len(frag_split[i]) > len(up_seq) and len(frag_split[i + 1]) > len(down_seq):
                    (up_seq, down_seq) = frag_split[i:i+2]
            (targetDict[target_id]["up_seq"], targetDict[target_id]["down_seq"]) = (up_seq, down_seq)
            # targetDict[target_id]["sample_msCount"] = {} #Place holder for sample msCounts in targetDict
    return targetDict


def import_targetDict_from_bed(bed_file):
    '''Import BED file into targetDict (compatible with RETrace2 BED format)'''
    targetDict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = fields[1]
            end = fields[2]

            try:
                motif_len = int(fields[3])
                num_copies = int(fields[4])
                target_id = fields[5]  # Use full target ID as in RETrace2

                # Create targetDict entry compatible with downsampling logic
                # We only need the keys to match alleleDict - the values structure is flexible
                targetDict[target_id] = {
                    "chrom": chrom,
                    "chromStart": start,
                    "chromEnd": end,
                    "motif_len": motif_len,
                    "num_copies": num_copies
                }

            except (ValueError, IndexError):
                # Skip if we can't parse the fields correctly
                continue

    return targetDict

def downsampleTarget(alleleDict_file, target_input, numTargets, output_alleleDict, output_sample_list, random_seed=None, use_bed=False):
    '''
    This script is modified from the orginal RETrace repo but now supports BED file input for consistency with the rest of the RETrace2 pipeline.
    It will downsample the alleleDict file such that samples have the user-specified numTargets.
    If samples have < numTargets, remove from alleleDict["samples"]
    Also outputs a sample list containing only samples that remain after downsampling.
    '''
    
    # Set random seed for reproducible results (if specified)
    if random_seed is not None:
        random.seed(random_seed)
        print(f"Setting random seed to {random_seed} for reproducible downsampling results")
        print("NOTE: Using fixed seed - results will be identical across runs")
    else:
        print("Using random downsampling - results will vary between runs")
        print("TIP: Use --random_seed for reproducible results, omit for robustness testing")

    alleleDict = pickle.load(open(alleleDict_file, 'rb'))
    if use_bed:
        targetDict = import_targetDict_from_bed(target_input)
    else:
        targetDict = import_targetDict(target_input)

    #We want to run through the alleleDict to obtain initial list of target_ids per sample
    sampleTargets = {} #Keep track of the targets that were captured per sample
    for target_id in sorted(alleleDict.keys()):
        if target_id in targetDict.keys(): #We only want to analyze target_id specified in targetDict
            for sample in sorted(alleleDict[target_id]["sample"].keys()):
                if sample not in sampleTargets.keys():
                    sampleTargets[sample] = {}
                    sampleTargets[sample]["original"] = set()
                sampleTargets[sample]["original"].add(target_id)

    #Downsample target_ids in sampleTargets
    for sample in sampleTargets.keys():
        if len(sampleTargets[sample]["original"]) >= numTargets:
            sampleTargets[sample]["downsample"] = set(random.sample(list(sampleTargets[sample]["original"]), numTargets))
        else:
            sampleTargets[sample]["downsample"] = set() #Create empty set if sample does not satisfy numTargets requirement

    #Run through alleleDict again and remove samples from calls
    for target_id in sorted(alleleDict.keys()):
        original_numSamples = len(alleleDict[target_id]["sample"].keys())
        for sample in sorted(alleleDict[target_id]["sample"].keys()):
            if target_id not in sampleTargets[sample]["downsample"]:
                del alleleDict[target_id]["sample"][sample]

    #Collect remaining samples after downsampling and write sample list
    remaining_samples = set()
    for target_id in sorted(alleleDict.keys()):
        remaining_samples.update(alleleDict[target_id]["sample"].keys())
    
    #Write sample list (one sample per line, compatible with build_phylo.py)
    with open(output_sample_list, 'w') as f:
        for sample in sorted(remaining_samples):
            f.write(f"{sample}\n")
    print(f"Written {len(remaining_samples)} samples to sample list: {output_sample_list}")

    #Save alleleDict into output file
    pickle.dump(alleleDict, open(output_alleleDict, 'wb'))

def main():
    parser = argparse.ArgumentParser(description="Randomly downsample targets from alleleDict and output sample list")
    parser.add_argument('--alleleDict', action="store", dest="alleleDict_file", help="Original alleleDict file")
    parser.add_argument('--numTargets', action="store", dest="numTargets", type=int, help="Number of targets for random downsampling")
    parser.add_argument('--output_alleleDict', action="store", dest="output_alleleDict", help="Output alleleDict file containing only randomly downsampled targets")
    parser.add_argument('--output_sample_list', action="store", dest="output_sample_list", help="Output sample list file containing samples that remain after downsampling (compatible with build_phylo.py)")
    parser.add_argument('--random_seed', action="store", dest="random_seed", type=int, default=None, help="Random seed for reproducible downsampling results. If not specified, uses true randomness for robustness testing")

    # Mutually exclusive group for target input formats
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument('--target_info', action="store", dest="target_info", help="Target_info file location (legacy format)")
    target_group.add_argument('--target_bed', action="store", dest="target_bed", help="Target BED file location (RETrace2 format)")

    args = parser.parse_args()

    # Determine which format to use
    if args.target_bed:
        use_bed = True
        target_input = args.target_bed
    else:
        use_bed = False
        target_input = args.target_info

    downsampleTarget(args.alleleDict_file, target_input, args.numTargets, args.output_alleleDict, args.output_sample_list, args.random_seed, use_bed)

if __name__ == "__main__":
    main()
