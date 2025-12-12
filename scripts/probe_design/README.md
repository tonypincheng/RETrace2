# RETrace v2 Probe Design

**Last updated:** Tony 2022-08-31

## Quick Start: 1nt Probe Design

### Step 1: Find Mononucleotide MS Locations
Run `1_find_mono-nucleotide.sh` to identify 1nt microsatellite (MS) locations in the reference genome using a custom script.

**Output:**
- `a.out`
- `hg19_MS_location.mono-nucleotide.txt`

### Step 2: Select Restriction Enzyme Fragments
Run `2_RE_selector_command.sh` to get restriction enzyme fragments containing MS.

**Output:** FASTA file with MS location + MS-containing fragment sequences for probe design

### Step 3: Probe Design
Create directory and prepare environment:
```bash
mkdir 3_probe_design
cd 3_probe_design
cp ../target_filter.sh ../seqkit .
conda activate retrace  # or environment with seqkit and seqtk installed
```

#### 3.1 Filter and Split Targets
Run `target_filter.sh` to:
- Filter by length (>80 or 120bp fragments)
- Perform alignment check
- Split targets into 1000 target chunks for IDT probe design input

**Output:**
- `filtered fragSeq.MseI.fasta`
- Split FASTA files

**Note:** Alternatively, you can send the whole `filtered fragSeq.MseI.fasta` to IDT and let them process the big file (no need to split).

#### 3.2 IDT Probe Design
Run IDT probe design tool with updated parameters:
- Change default target length to 120bp
- v2 tool has optimized tiling method

**Output:** Combined `.xlsx` file

#### 3.3 Generate Oligo Pool Order
Run `Oligo_orderPool.ipynb` to convert IDT output to oligo pool order format.

**Note:** IDT v2 output format has changed. This script is modified from Chris's `Custom-Array_Order.v2.py`.

---

## RETrace v1 Background (Chris's Original)

### Goal
Determine the number of microsatellite loci that can be captured per restriction enzyme digestion, now including mononucleotide repeats (increased mutation rate - Boyer et al, Human Mol Gen 2002).

### Pipeline

1. **Download trf data for hg19**
   - Source: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.trf.bed.gz

2. **Reformat BED file**
   - Keep only: bin, chrom, chromStart, chromEnd, name
   - Filtering criteria:
     - Subunit length: 1-6
     - Only perfect repeats (col 8 of trf bed must equal 100)
     - Total MS length <100bp (avoids difficult-to-sequence fragments)

3. **Run RE_selector-v1.1.py**
   - Fragment the genome
   - Output sequences resulting from various restriction enzyme digestions

4. **Calculate statistics**
   - Count each microsatellite subunit length in outputted `fragSeq.fasta` for each restriction enzyme

### Important Note: Mononucleotide Detection Issue

The trf data is missing many mononucleotides (especially G's and C's). RepeatMasker data has the same issue. 

**Solution:** Manual search using custom C script `find_mono-nucleotide.c` (adapted from https://gist.github.com/eade421ffda5e0c71055b514a8f86822.git)

**Modifications:**
1. Changed requirement from 5 to **≥15 mononucleotides** in a row (targets longer sequences)
2. Output format matches `hg19_MS_location.txt` format (bin, chrom, chromStart, chromEnd, name)

**Usage:**
```bash
gcc find_mono-nucleotide.c  # Creates "a.out" compiled program
cat /path/to/hg19_reference.fa | ./a.out > hg19_MS_location.mono-nucleotide.txt
cat hg19_MS_location.mono-nucleotide.txt hg19_MS_location.txt > hg19_MS_location.all.txt
```

### References
- https://www.biostars.org/p/267241/
- Old analysis: `/media/NAS3/Chris.NAS3/genome_analyzer-3TB_slot2/Fate_Mapping/pp_design/microsatellite_locations/hg19_all_MS_filtering_new`
- RE selection: `/media/NAS3/Chris.NAS3/genome_analyzer-6TB_slot4/Fate_Mapping/MS+Methyl/RE_selection`
