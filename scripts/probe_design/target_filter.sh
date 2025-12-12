# Update 2/2023
# No need to split frag file, can ask IDT to process the large file for you
for f in ../2_RE_selector_1nt_10bp-max/fragSeq.MseI.fasta
do
        fbname=$(basename "$f" .fasta)
        ./seqkit seq -m 120 $f >$fbname.rm-short.fasta # min 120bp filter show hihger on target %
        sed -i 's/\:/\_/g' $fbname.rm-short.fasta
        bwa mem -t 10 /media/Scratch_SSD_Voyager/picheng/GenomeDB/mm39/bwa-index/mm39.fa $fbname.rm-short.fasta | samtools view -Sbh -q1 - | samtools bam2fq - | seqtk seq -A > $fbname.filtered.fasta  

rm $fbname.rm-short.fasta

done


# Use following IDT design parameter: 
## Input format: FASTA
## Species: mm39
## Run probe QC: Yes
## Panel name: 1ntMS_10-50bp_mm39
## Probe length: 120
## Probe tiling density: 1x
## Maximum end gap: 30bp





# Chris previous script
# Update 9/2022: Tony add some comments 

# Step1
#for f in ../2_RE_selector_output-1nt_greater14/fragSeq.MseI.fasta
#do
#	fbname=$(basename "$f" .fasta)
#	./seqkit seq -m 120 $f >$fbname.rm-short.fasta # Tony update bc IDT retrict only 120bp or submit custom design
#	sed -i 's/\:/\_/g' $fbname.rm-short.fasta
#	bwa mem -t 10 /media/Scratch_SSD_Voyager/picheng/GenomeDB/hg38/bwa-index/hg38.fa $fbname.rm-short.fasta | samtools view -Sbh -q1 - | samtools bam2fq - | seqtk seq -A >$fbname.fasta # Tony: this is to do aligment check?
	
	#We need to split the output fasta file to smaller files (of 1000 targets each) in order to use xGEN probe online design portal (see <https://www.idtdna.com/site/order/ngs>)
#	split -l 1800 $fbname.fasta $fbname.split  # Tony:change split from 2000-->1800, filezie limit changed in IDT
#	for g in ./$fbname.split*
#	do
#		mv $g $g.fasta
#	done
#done


# When running the xGEN target capture probe design online tool, we want to use the following parameters: (https://www.idtdna.com/site/order/designtool/index/XGENDESIGN)
#	1) Input Format: FASTA Sequence
#	2) Target Species: hg38
#	3) Probe Length: 120
#	4) Probe Tiling Density: 1x
#Targets that are discarded are normally due to target being too short for the desired probe size (120bp)
#We do this manually at <https://www.idtdna.com/site/order/ngs> for each split fasta file and download the following two files:
#	1) NGS-Targets.xls = Target information for the probes
#	2) NGS_Design.zip = Design output of the probes (on average ~3 probes/target were designed)
#We then want to combine all the information into one single file describing NGS-Targets and NGS-Design


# Step2
#for f in ./*.xlsx; do fbname=$(basename "$f" .xlsx); ssconvert $f $fbname.NGS-Targets.csv; cat $fbname.NGS-Targets.csv >>NGS-Targets.csv; rm $fbname.NGS-Targets.csv; done



#for h in fragSeq.MseI.split*.NGS-Design.zip
#do
#	hbname=$(basename "$h" .NGS-Design.zip)
#	unzip $h
#	ssconvert Probe-Sequences.xls Probe-Sequences.csv
#	cat Probe-Sequences.csv >>NGS-Design.csv
#	rm Probe-Sequences*
#
#	ssconvert $hbname.NGS-Targets.xls $hbname.NGS-Targets.csv
#	cat $hbname.NGS-Targets.csv >>NGS-Targets.csv
#	rm $hbname.NGS-Targets.csv
#done
