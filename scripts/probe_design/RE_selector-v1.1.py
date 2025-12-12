#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import argparse
import os
from tqdm import tqdm
import mmap
from Bio import SeqIO
import re as regex
import matplotlib
matplotlib.use('Agg') #This sets the backend to not display to device (http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined)
import matplotlib.pyplot as plt
import seaborn as sns
import operator

'''
Usage: python script.py --bed input.bed --fasta input.fasta --RE_file RE.txt --outDir outDir

Example:
python RE_selector-v1.1.py \
        --bed ../1_find_mono-nucleotide/mm39_MS_location.mono-nucleotide_10bp-max.txt \
        --fasta /media/Scratch_SSD_Voyager/picheng/GenomeDB/mm39/raw_fasta/mm39.fa \
        --RE_file /media/Home_Raid1_Voyager/picheng/RETrace_REselection/RE_selection_Tony/NEB_RE.MseI.txt \
        --readLen 250 \
        --outDir ./

----------v1.0----------
We want to determine which restriction enzyme(s) will cut the microsatellite such that we obtain two different RE site flanking the MS.
Below are the statistics we want to calculate from the simulated capture:
    1) Number of MS captured per RE added (with the same cut site flanking the sequence)
    2) Number of MS captured by different RE cutting around flanking regions
    3) Plot of fragment length distribution from RE cut (to determine maxRange for the number of fragments captured within 200bp range)
----------v1.1----------
We want to determine which RE will work best paried with MspI to perform the dual sequencing approach of obtaining MS information along with methylation information.
Consequently, we must develop a new script that will perform the following:
    1) Determine MS captured by the secondary RE, output fasta file containing sequence around MS sites to be used to design MS-specific primers for RCA
    2) Plot separate fragment length distributions from the restriction enzyme used
'''
#%%
def get_num_lines(file_path): #For tqdm, we need the total number of lines in file <https://blog.nelsonliu.me/2016/07/29/progress-bars-for-python-file-reading-with-tqdm/>
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(),0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def importFasta(fasta_input):
    fastaDict = {}
    for fasta in SeqIO.parse(fasta_input, "fasta"):
        fastaDict[fasta.id] = {}
        fastaDict[fasta.id]["seq"] = str(fasta.seq).upper()
    return fastaDict

def MSLocator(fastaDict, bed_input):
    MSDict = {}
    with open(bed_input) as f:
        for line in f:
            if not line.startswith('#'):
                line = line.rstrip()
                (bin_name, chrom, chromStart, chromEnd, name) = line.split()
                (num_sub, sub_seq) = name.split('x')
                ms_seq = sub_seq*int(num_sub)
                ms_len = len(ms_seq)
                if fastaDict[chrom]["seq"][int(chromStart):int(chromStart)+ms_len] == ms_seq: #Check that given genomic location is correct
                    targetID = chrom + ':' + chromStart + '_' + str(int(chromStart)+ms_len) + '_' + name
                    if chrom not in MSDict.keys():
                        MSDict[chrom] = {}
                    MSDict[chrom][int(chromStart)] = targetID
                    MSDict[chrom][int(chromStart)+ms_len] = targetID
    return MSDict

def RELocator(fastaDict, RE_info, RE_name):
    '''We want to precompute the locations of cutsites for the given restriction enzyme across the genome'''
    RE_info[RE_name]['pos'] = {}
    for chrom in sorted(fastaDict.keys()):
        RE_info[RE_name]['pos'][chrom] = []
        for match in regex.finditer(RE_info[RE_name]['seq'], fastaDict[chrom]["seq"]):
            RE_info[RE_name]['pos'][chrom].append(match.start() + int(RE_info[RE_name]['cut']))
    return

def genomeCutter(RE1, RE2, MSDict, RE_info):
    '''Merge the RE1, RE2, and MS position lists into a single position list (total_pos_list) along with corresponding characteristics of each position (pos_item_list) for each chrom'''
    posDict = {}
    print("Positioning RE cuts")
    for chrom in tqdm(sorted(MSDict.keys())):
        posDict[chrom] = {}
        RE1_pos = sorted(RE_info[RE1]['pos'][chrom])
        RE2_pos = sorted(RE_info[RE2]['pos'][chrom])
        MS_pos = sorted(MSDict[chrom].keys())

        total_pos_list = [] #This contains all the positions in order
        pos_item_list = [] #This will be a list containing ordered positional items (i.e. [RE1, MS1, MS1, RE2, RE2...])

        RE1_index = 0 #We want to keep a running count for the index position of each list
        RE2_index = 0
        MS_index = 0

        while RE1_index<len(RE1_pos) or RE2_index<len(RE2_pos) or MS_index<len(MS_pos):
            comparison_list = []
            if RE1_index<len(RE1_pos):
                comparison_list.append(RE1_pos[RE1_index])
            else:
                comparison_list.append(float("inf")) #If we have reached the end of the RE1_pos list, then we want to enter a dummy "inf" variable
            if RE2_index<len(RE2_pos):
                comparison_list.append(RE2_pos[RE2_index])
            else:
                comparison_list.append(float("inf"))
            if MS_index<len(MS_pos):
                comparison_list.append(MS_pos[MS_index])
            else:
                comparison_list.append(float("inf"))

            min_comparison = min(comparison_list) #We want to determine which is the next subsequent event in the chrom
            pos_item = []
            if min_comparison == comparison_list[0]:
                pos_item.append("RE_" + RE1)
                RE1_index += 1
            if min_comparison == comparison_list[1]:
                pos_item.append("RE_" + RE2)
                RE2_index += 1
            if min_comparison == comparison_list[2]:
                pos_item.append("MS_" + MSDict[chrom][MS_pos[MS_index]])
                MS_index += 1
            pos_item_list.append('|'.join(pos_item))
            total_pos_list.append(min_comparison)
        posDict[chrom]['pos_item_list'] = pos_item_list
        posDict[chrom]['total_pos_list'] = total_pos_list
    return posDict

def findFrag(RE1, RE2, posDict, readLen, fastaDict, outDir):
    '''This script will determine the fragment sizes across the whole reference genome along with those that contain suitable MS (that can be read by a read with given readLen)'''
    max_fragSize = 2000
    fragDict = {}
    fragDict['ref_' + RE1], fragDict['ref_' + RE2], fragDict['ref_diff'] = {}, {}, {}
    fragDict['MS_' + RE1], fragDict['MS_' + RE2], fragDict['MS_diff'] = {}, {}, {}
    fragList = []

    output_seq = open(outDir + '/fragSeq.' + RE1 + '.fasta', 'a')
    print("Finding Fragments")
    for chrom in tqdm(sorted(posDict.keys())):
        RE_pos, RE_ID, MS_pos, MS_ID = ([] for i in range(4)) #Initialize multiple lists in same line (see <https://stackoverflow.com/questions/2402646/python-initializing-multiple-lists-line>)
        RE_pos.append(0) #Initialize first RE cut in chromosome
        RE_ID.append('Start')
        for i in range(0, len(posDict[chrom]['pos_item_list'])-1):
            if 'RE' in posDict[chrom]['pos_item_list'][i]:
                RE_pos.append(posDict[chrom]['total_pos_list'][i])
                RE_ID.append(posDict[chrom]['pos_item_list'][i].replace('RE_',''))
            elif 'MS' in posDict[chrom]['pos_item_list'][i]:
                MS_pos.append(posDict[chrom]['total_pos_list'][i])
                MS_ID.append(posDict[chrom]['pos_item_list'][i].replace('MS_',''))
                if '|' in posDict[chrom]['pos_item_list'][i]:
                    (MS_pos, MS_ID) == ([] for i in range(2)) #We want to reset if there is a restriction enzyme that cuts at the same position as the start/end of MS
            if len(RE_pos)==2: #We want to determine every fragment as defined as the portion of DNA that is in between two RE cuts (regardless of whether there is MS in between)
                frag_size = RE_pos[1]-RE_pos[0]
                if frag_size < max_fragSize:
                    fragList.append(frag_size)
                    #Determine fragments with different RE flanking the sequence
                    if len(set(RE_ID))==2:
                        try:
                            fragDict['ref_diff'][frag_size] += 1
                        except:
                            fragDict['ref_diff'][frag_size] = 1
                    #Determine fragments with the same RE flanking the sequence
                    elif len(set(RE_ID))==1:
                        try:
                            fragDict['ref_' + RE_ID[0]][frag_size] += 1
                        except:
                            fragDict['ref_' + RE_ID[0]][frag_size] = 1
                    '''Check whether fragment contains a single microsatellite that can be read using given readLen'''
                    if len(MS_pos)==2 and len(set(MS_ID))==1:
                        if (MS_pos[1]-RE_pos[0])<=readLen and (RE_pos[1]-MS_pos[0])<=readLen:
                            #Determine MS fragments with different RE flanking the sequence
                            if len(set(RE_ID))==2: #Different RE cut both sides of flanking region
                                try:
                                    fragDict['MS_diff'][frag_size] += 1
                                except:
                                    fragDict['MS_diff'][frag_size] = 1
                            #Determine fragments with the same RE flanking the sequence.  We also want to output the fragment seq
                            elif len(set(RE_ID))==1:
                                try:
                                    fragDict['MS_' + RE_ID[0]][frag_size] += 1
                                except:
                                    fragDict['MS_' + RE_ID[0]][frag_size] = 1
                                if RE_ID[0] == RE1:
                                    fragSeq = fastaDict[chrom]["seq"][RE_pos[0]:RE_pos[1]]
                                    output_seq.write('>' + MS_ID[0].replace('MS_','') + "\n" + fragSeq + "\n")
#                                print(fragSeq + "\n" + str(frag_size) + ',' + str(fragDict['MS_' + RE_ID[0]][frag_size]) + "\t" + ','.join(str(x) for x in RE_pos) + "\t" + ','.join(str(y) for y in MS_pos) + "\t" + ','.join(RE_ID) + "\t" + ','.join(MS_ID))
                RE_pos, RE_ID, MS_pos, MS_ID = ([] for i in range(4))
                RE_pos.append(posDict[chrom]['total_pos_list'][i]) #Re-initialize RE_pos to start off where the last cut was
                RE_ID.append(posDict[chrom]['pos_item_list'][i].replace('RE_',''))
    output_seq.close()
    return fragDict, fragList

def calcStats(RE1, RE2, fragDict, fragList, outDir):
    '''Return the following stats for RE1: (RE1 is the sticky end to which we ligate on our MS adapters)
        1) Total microsatellite fragments
        2) Total fragments from reference genome
        3) 200bp range with most MS fragments (number MS fragments; enrichment)
        4) 200bp range with highest enrichment of MS fragments (number MS fragments; enrichment)
    Return the following stats for RE2: (RE2 is MspI, which is used to capture methylation data)
        1) Total fragments cut with RE2 on both sides (this is number maximum CpGs we can sequence)
        2) 200bp range with most RE2 fragments
    '''
    output = open(outDir + '/RE_stats.txt', 'a')
    rangeDict = {}
    for key in fragDict.keys():
        print(key + "\t" + str(sum(fragDict[key].values())))
        rangeDict[key] = {}
        if 'MS' in key and 'diff' not in key:
            rangeDict['enrichment_' + key] = {}
    print("Calculating Stats")
    for i in tqdm(range(min(fragList), max(fragList)-200)):
        range_name = "[" + str(i) + "," + str(i+200) + ")"
        #Initialize rangeDict range_name to calculate number of fragments in each size range
        for fragSize in sorted(set(fragList)):
            if fragSize>=i and fragSize<i+200:
                for key in fragDict.keys():
                    if fragSize in fragDict[key].keys():
                        try:
                            rangeDict[key][range_name] += fragDict[key][fragSize]
                        except:
                            rangeDict[key][range_name] = fragDict[key][fragSize]
            elif fragSize>=i+200:
                break
        for RE in [RE1,RE2]:
            try:
                rangeDict['enrichment_MS_' + RE][range_name] = rangeDict['MS_' + RE][range_name]/(rangeDict['MS_' + RE][range_name] + rangeDict['ref_' + RE][range_name])
            except:
                continue

    #Calculate RE1 (microsatellite capture) stats
    frag_RE1 = sum(fragDict['MS_' + RE1].values())
    total_ref_RE1 = sum(fragDict['ref_' + RE1].values())
    total_enrichment_RE1 = frag_RE1/(frag_RE1 + total_ref_RE1)
    maxRange_frag_RE1 = max(rangeDict['MS_' + RE1], key = rangeDict['MS_' + RE1].get)
    maxRange_enrichment_RE1 = max(rangeDict['enrichment_MS_' + RE1], key = rangeDict['enrichment_MS_' + RE1].get)
#    print(maxRange_frag_RE1 + "\t" + str(rangeDict['MS_' + RE1][maxRange_frag_RE1]))
#    print(maxRange_enrichment_RE1 + "\t" + str(rangeDict['enrichment_MS_' + RE1][maxRange_enrichment_RE1]))

    frag_diff = sum(fragDict['MS_diff'].values())
    total_ref_diff = sum(fragDict['ref_diff'].values())

    #Calculate RE2 (methylation information) stats
    frag_RE2 = sum(fragDict['ref_' + RE2].values())
    maxRange_frag_RE2 = max(rangeDict['ref_' + RE2], key = rangeDict['ref_' + RE2].get)

    output.write(RE1 + "\t" + str(frag_RE1) + "\t" + str(total_ref_RE1) + "\t" + str("%.4f" % total_enrichment_RE1) + "\t" +
                 str(frag_diff) + "\t" + str(total_ref_diff) + "\t" +
                 maxRange_frag_RE1 + "(" + str(rangeDict['MS_' + RE1][maxRange_frag_RE1]) + "," + str("%.4f" % rangeDict['enrichment_MS_' + RE1][maxRange_frag_RE1]) + ")" + "\t" +
                 maxRange_enrichment_RE1 + "(" + str(rangeDict['MS_' + RE1][maxRange_enrichment_RE1]) + "," + str("%.4f" % rangeDict['enrichment_MS_' + RE1][maxRange_enrichment_RE1]) + ")" + "\t" +
                 str(frag_RE2) + "\t" + maxRange_frag_RE2 + "(" + str(rangeDict['ref_' + RE2][maxRange_frag_RE2]) + ")\n")
    output.close()

def plotFrag(fragDict, RE, outDir):
    '''Plot fragment sizes for the following:
        1) MS_RE fragments = those that capture the microsatellite region
        2) ref_RE fragments = all fragments that were cut out with the same given RE on both flanking regions
        3) ref_MspI_HpaII fragments = fragments cut with same MspI cut site on flanking regions; these are the fragment sizes containing desired methylation information
    '''
    print("Plotting Fragments")
    #In order to plot the distribution of fragment sizes, we want to make a new 1D list containing fragSize for all targets (only care about and fragDict['MS_diff'])
    fragList_MS_RE = []
    for fragSize in sorted(fragDict['MS_' + RE].keys()):
        fragList_MS_RE.extend([int(fragSize)]*fragDict['MS_' + RE][fragSize])
    fragList_ref_RE = []
    for fragSize in sorted(fragDict['ref_' + RE].keys()):
        fragList_ref_RE.extend([int(fragSize)]*fragDict['ref_' + RE][fragSize])
    fragList_ref_MspI_HpaII = []
    for fragSize in sorted(fragDict['ref_MspI_HpaII'].keys()):
        fragList_ref_MspI_HpaII.extend([int(fragSize)]*fragDict['ref_MspI_HpaII'][fragSize])
    if len(fragList_MS_RE)>100 and len(fragList_ref_RE)>100 and len(fragList_ref_MspI_HpaII)>100:
        #Plot the KDE distribution using seaborn
        fig = plt.figure()

        ax1 = fig.add_subplot(3,1,1)
        ax1.set_title("Microsatellite fragments with " + RE + " sticky ends")
        ax1.set_xlim([0,2000])
        ax1.locator_params(nbins=5, axis='y')
        sns.distplot(fragList_MS_RE, bins=40, kde=False, hist=True, norm_hist=False, ax=ax1)

        ax2 = fig.add_subplot(3,1,2)
        ax2.set_title("All fragments with " + RE + " sticky ends")
        ax2.set_xlim([0,2000])
        ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax2.locator_params(nbins=5, axis='y')
        sns.distplot(fragList_ref_RE, bins=40, kde=False, hist=True, norm_hist=False, ax=ax2)

        ax3 = fig.add_subplot(3,1,3)
        ax3.set_title("All fragments with MspI_HpaII sticky ends")
        ax3.set_xlim([0,2000])
        ax3.set_xlabel("Fragment Sizes")
        ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax3.locator_params(nbins=5, axis='y')
        sns.distplot(fragList_ref_MspI_HpaII, bins=40, kde=False, hist=True, norm_hist=False, ax=ax3)

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tight_layout()
        fig.savefig(outDir + '/fragSize.' + RE + ".png")
        plt.close()

#%%
def main():
    parser = argparse.ArgumentParser(description="Determine fragment sizes cut around MS using different RE")
    parser.add_argument('--bed', action="store", dest="bed_input")
    parser.add_argument('--fasta', action="store", dest="fasta_input")
    parser.add_argument('--RE_file', action="store", dest="RE_file")
    parser.add_argument('--readLen', action="store", dest="readLen", default=250) #Default set by HiSeq RapidRun 2x250bp PE reads
    parser.add_argument('--outDir', action="store", dest="outDir")
    args = parser.parse_args()

    '''Import fasta file'''
    print("Import Fasta")
    fastaDict = importFasta(args.fasta_input)

    '''Locate MS positions in each chromosome'''
    print("Import Microsatellites")
    MSDict = MSLocator(fastaDict, args.bed_input)

    '''Make outDir (will create *.png files containing plots of each RE combo)'''
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    '''Pre-compute the restriction enzyme cut positions'''
    print("Import RE Cut Sites")
    RE_info = {}
    with open(args.RE_file) as f:
        for line in tqdm(f, total=get_num_lines(args.RE_file)):
            (RE_name, RE_seq) = line.split()
            RE_info[RE_name] = {}
            RE_info[RE_name]['cut'] = RE_seq.index('|')
            RE_info[RE_name]['seq'] = RE_seq.replace('|','')
            RELocator(fastaDict, RE_info, RE_name)

    output = open(args.outDir + '/RE_stats.txt', 'w')
    output.write("RE_name\tTotal MS Frag\tTotal Ref Frag\tTotal Enrichment\tTotal MS Frag (Diff RE)\tTotal Ref Frag (Diff RE)\tRange w/ Max MS (Num MS, Enrichment)\tRange w/ Max Enrichment (Num MS, Enrichment)\tTotal MspI Frag\tRange w/ Max MspI Frag (Num Frag)\n")
    output.close()
    '''Determine how each RE pair captures MS loci'''
    RE_analyzed = []
    for RE in sorted(RE_info.keys()):
        if (RE != "MspI_HpaII") and (RE_info[RE]['seq'] not in RE_info["MspI_HpaII"]['seq']) and (RE_info["MspI_HpaII"]['seq'] not in RE_info[RE]['seq']):
            RE_name = ','.join(sorted([RE,"MspI_HpaII"]))
            if RE_name not in RE_analyzed:
                print("----------Analyzing: " + RE_name + "----------")
                RE_analyzed.append(RE_name)
                posDict = genomeCutter(RE, "MspI_HpaII", MSDict, RE_info)
                fragDict, fragList = findFrag(RE, "MspI_HpaII", posDict, int(args.readLen), fastaDict, args.outDir)
                # NOTE: Update 1/3/2020: I commented out the following section of code in order to skip calculation of stats in order to speed up processing
                if sum(fragDict['MS_' + RE].values()) > 100:
                    calcStats(RE, "MspI_HpaII", fragDict, fragList, args.outDir)
                    plotFrag(fragDict, RE, args.outDir)
                else:
                    print("Too few MS fragments:\t" + str(sum(fragDict['MS_' + RE].values())))

#%%
if __name__ == "__main__":
    main()
