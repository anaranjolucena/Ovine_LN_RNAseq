##########################################
#  Ovine Lymph Node Liver Fluke RNA-seq  #
##########################################
 
# Author: Amalia Naranjo
# Last updated on: 31/03/2020

########################
# Perform MD5 checksum #
########################

# Enter working directory:
cd /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine

# Perform md5sum check:
md5sum -c values.md5sum >> \
/home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/md5check_UCD.txt

# Check that all files passed the check:
grep -c 'OK' md5check_UCD.txt

###########################################
# FastQC quality check of raw FASTQ files #
###########################################


# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering
cd !$

# Run FastQC in one file to check if it works:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering \
--noextract --nogroup -t 20 \
/home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz

# Transfer compressed folder to personal laptop via SCP
# and check the HTML report:
scp alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/N12_S29_L002_R1_001_fastqc.zip .

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine \
-name *.fastq.gz`; \
do echo "fastqc --noextract --nogroup -t 20 \
-o /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering $file" \
>> fastqc.sh; \
done

# Check number of lines in script:
wc -l fastqc.sh

# Split and run all scripts on Rodeo:
split -d -l 19 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Analysis complete" >> succesful_fastqc.txt
done

# Check number of lines in created document
wc -l succesful_fastqc.txt

# Deleted all HTML files:
rm -r *.html

# Collect FastQC stats:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/tmp


for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/tmp; \
done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> summary_pre-filtering.txt; \
done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Check sequence quality:
grep 'Per base sequence quality' summary_pre-filtering.txt >> seq_quality.txt
wc -l seq_quality.txt # 38 lines
grep PASS seq_quality.txt | wc -l # 38 lines
grep WARN seq_quality.txt | wc -l # 0 line

# Check if sequence contain adapters:
grep 'Adapter Content' summary_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt # 38 lines
grep PASS adapter_content.txt | wc -l # 30 lines
grep WARN adapter_content.txt | wc -l # 8 lines

# Check sequence length:
grep 'Sequence length' basic_stats_pre-filtering.txt >> seq_length.txt
wc -l seq_length.txt # 38 lines
less seq_length.txt # 151 for all

# Check per base sequence quality
grep 'Per base sequence quality' summary_pre-filtering.txt >> base_quality.txt
wc -l base_quality.txt # 38 lines
grep PASS base_quality.txt | wc -l # 38 lines
grep WARN base_quality.txt | wc -l # 0 lines


# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line)
# and check HTML reports:
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/pre-filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir /home/workspace/alucena/ovineLN_RNAseq/filt_fastq
cd !$


# I put file with adapters in this directory from my terminal (not from rodeo) (Using only first 13 bases from adapter)
scp -r Illumina_PE_adapters13bp.txt alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/filt_fastq

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 151 \
-pe1 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
-pe2 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz \
-o /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters13bp.txt -5a_mp 100 -5a_del 0 \
-5a_ins 0 -5a_fmi 140 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create a bash script to perform filtering of each pair of FASTQ files (to do all samples):
for file in `find /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/ \
-name *R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_L002_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 151 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters13bp.txt -5a_mp 100 \
-5a_del 0 -5a_ins 0 -5a_fmi 140 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.sh; done

# Split and run all scripts on Rodeo:
split -d -l 10 filtering.sh filtering.sh.
for script in `ls filtering.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all files were processed:
for file in `ls filtering.sh.0*.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files with discarded reads:
# Do bash script:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/ \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done

#Run script:
for script in `ls discarded_compression.sh`; \
do 
chmod 755 $script 
nohup ./$script &
done

# Gather ngsShoRT reports from all samples into one file:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/ \
-name final_PE_report.txt`; \
do echo echo \
"\`dirname $file | perl -pe 's/.*(N\d+_S\d\d)/\$1/'\` \
\`grep 'Read Pair Count:' $file\` \
\`grep 'Removed PE Pair\* Count:' $file\` >> \
ngsshort.txt" >> ngsshort_summary.sh
done

chmod 755 ngsshort_summary.sh
./ngsshort_summary.sh

# Add header to facilitate wrangling in R later on: 
sed -i $'1 i\\\nSample_name X1 X2 X3 Read_pair_count X4 X5 X6 X7 Removed_PE_pair_count Percent_removed' ngsshort.txt


# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line)
# and check HTML reports:
scp \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/ngsshort.txt .


################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering
cd !$


# Run FastQC in one sample to see if it's working well:
#For R1:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R1_001.fastq.gz

#In new rodeo tab, R2:

fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R2_001.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line)
# and check HTML reports:
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/*fastqc.zip .

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/ \
-name *_L002_R*.fastq.gz`; \
do echo "fastqc --noextract --nogroup -t 20 \
-o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering $file" \
>> fastqc_filt.sh; \
done

# Check number of lines in script:
wc -l fastqc_filt.sh

# Split and run all scripts on Rodeo:
split -d -l 19 fastqc_filt.sh fastqc_filt.sh.
for script in `ls fastqc_filt.sh.0*`;
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
for file in `ls fastqc_filt.sh.0*.nohup`; \
do more $file | grep "Analysis complete" >> succesful_fastqc.txt
done

# Check number of lines in created document
wc -l succesful_fastqc.txt

# Deleted all HTML files:
rm -r *.html

#Send to laptop
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/*fastqc.zip .

# Check all output from FastQC:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/tmp; \
done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/tmp/*_fastqc/ \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/tmp/*_fastqc/ \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

#Remove tmp directory
rm -r tmp

#Transfer to laptop
scp \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/*.txt .

##############################################################################
# Alignment of FASTQ files against the Ovis aries reference genome with STAR #
##############################################################################


# Required software is STAR 2.7.3a, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf


# Download Ovis aries reference genome, Oar rambouillet version 1.0
mkdir -p /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/source_file
cd !$
nohup wget ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCF_002742125.1_Oar_rambouillet_v1.0/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.fna.gz &

#unzip file
gunzip GCF_002742125.1_Oar_rambouillet_v1.0_genomic.fna.gz

# Download annotation file for NCBI Ovis aries Annotation Release 103

mkdir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/annotation_file
cd !$
wget ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCF_002742125.1_Oar_rambouillet_v1.0/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gtf.gz
gunzip GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gtf.gz


# Generate genome indexes files using annotations:
mkdir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150
cd !$


nohup STAR --runThreadN 40 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--genomeFastaFiles \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/source_file/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.fna \
--sjdbGTFfile /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/annotation_file/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gtf \
--sjdbOverhang 150 \
--outFileNamePrefix \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 &


# Create and enter alignment working directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment
cd !$

# Mapping reads from one FASTQ file to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R1_001.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R2_001.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/ \
-name *R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/trimmed_(N\d+_S\d\d_)L002_R1_001.fastq.gz/$1/'`; \
foldername=`echo $sample | perl -p -e 's/(N\d+_S\d\d)_$/$1/'`; \
echo "mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/$foldername; \
cd /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn $file $file2 \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./$sample \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx" \
>> alignment.sh; done

# Split and run all scripts on Rodeo:
split -d -l 10 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`;
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup
# 10: It worked
grep -c 'finished successfully' alignment.sh.01.nohup
# 9: It worked

#Copy script for perl, in alignment directory:
cp /home/workspace/acampos/flukeRNAseq/star_report_opener.pl .
chmod 755 star_report_opener.pl

# Merge all STAR log.final.out files into a single file:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment \
-name *Log.final.out`; \
do perl /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/star_report_opener.pl -report $file; done;

#Transfer to personal laptop
scp \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/All_star_log_final_out.txt .


#############################################
# FastQC quality check of aligned BAM files #
#############################################


# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/


#Create directory in quality_check for FASTQC of BAM file
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/1st_filtered/N12_S29Aligned.out.bam

# Transfer compressed folder to personal laptop via SCP (tab in personal laptop terminal)
# and check the HTML report:
scp alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/N12_S29Aligned.out_fastqc.zip .

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/*_S*/ \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 20 \
-o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts on Rodeo
split -d -l 10 fastqc_aligned.sh fastqc_aligned.sh.
chmod 755 fastqc_aligned.sh.*

for script in `ls fastqc_aligned.sh.*`;
do
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/tmp; \
done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/tmp/*_fastqc/ \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/tmp/*_fastqc/ \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Secure copy tmp folder to laptop (share folder) and rename folder as post_alignment.
# I use my own laptop (out of Rodeo)
scp -r alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/tmp .

# Remove temporary folder in Rodeo:
rm -r tmp/


#########################
# Calculate insert size #
#########################

#Software employed 
### Picard version 2.18.27: java -jar /usr/local/src/picard/build/libs/picard.jar

# Create and go to working directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/insert_size
cd !$


# Sort aligned BAM files:
for file in \
`find /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment \
-name *Aligned.out.bam`; \
do outfile=`basename $file | perl -p -e 's/Aligned.out.bam/Sorted.out.bam/'`; \
echo "java -jar /usr/local/src/picard/build/libs/picard.jar \
SortSam I=$file O=$outfile SORT_ORDER=coordinate" >> sort.sh; \
done

# Split and run all scripts on Rodeo
split -d -l 10 sort.sh sort.sh.
chmod 755 sort.sh.*

for script in `ls sort.sh.*`;
do
nohup ./$script > ${script}.nohup &
done


# Collect insert sizes:
for file in `ls *_Sorted.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Sorted.out.bam//'`; \
echo "java -jar /usr/local/src/picard/build/libs/picard.jar \
CollectInsertSizeMetrics \
I=$file \
O=${sample}_insert_size_metrics.txt \
H=${sample}_insert_size_histogram.pdf M=0.5" >> collect_insert_size.sh; \
done

# Split and run all scripts on Rodeo
chmod 755 collect_insert_size.sh

for script in `ls collect_insert_size.sh`;
do
nohup ./$script > ${script}.nohup &
done


# Collect insert size metrics for all samples into one file:
for file in `ls /home/workspace/alucena/ovineLN_RNAseq/insert_size/*_insert_size_metrics.txt`; \
do sample=`basename $file | perl -p -e 's/_insert_size_metrics.txt//'`; \
stats=`sed -n '/MEDIAN_INSERT_SIZE/{n;p;}' $file`; \
printf "${sample}\t${stats}\n" >> All_insert_size.txt; \
done

wc -l All_insert_size.txt # 19 lines

# Add header to summary stats file:
header=`grep 'MEDIAN_INSERT_SIZE' /home/workspace/alucena/ovineLN_RNAseq/insert_size/N12_S29_insert_size_metrics.txt`; \
sed -i $"1 i\Sample_id\t${header}" \
All_insert_size.txt

wc -l All_insert_size.txt # 20 lines


# Transfer stats to laptop:
scp \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/insert_size/All_insert_size.txt .

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################


# Required package is featureCounts, which is part of Subread software (featureCounts v2.0.0),
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf


# Download GFF annotation file for NCBI Ovis aries Annotation Release 103 (GTF file not working with featureCounts)

cd /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/annotation_file
wget ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCF_002742125.1_Oar_rambouillet_v1.0/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gff.gz
gunzip GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gff.gz


# Create working directories:
cd /home/workspace/alucena/ovineLN_RNAseq
mkdir -p Count_summarisation/sense
cd /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense


# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/annotation_file/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gff \
-B -p -C -R BAM -s 2 -T 15 -t gene -g Dbxref -F GFF -o ./counts.txt \
/home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/N12_S29/N12_S29_Aligned.out.bam

# Bash script to run featureCounts in all samples:
for file in `find /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/$sample; \
cd /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/$sample; \
featureCounts -a \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/annotation_file/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.gff \
-B -p -C -R BAM -s 2 -T 20 -t gene -g Dbxref -F GFF \
-o ${sample}_counts.txt $file" >> counts.sh; \
done


# Split and run all scripts on Rodeo
split -d -l 10 counts.sh counts.sh.
chmod 755 counts.sh.*

for script in `ls counts.sh.*`;
do
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Summary of counting results can be found in file' counts.sh.00.nohup # 10 lines
grep -c 'Summary of counting results can be found in file' counts.sh.01.nohup # 9 lines

# Copy all *counts.txt files to temporary folder:
mkdir /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp

for file in `find /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/ \
-name *_counts.txt`; do cp $file \
-t /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp; \
done

# Transfer count files to laptop:
scp -r alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp .

# Remove tmp folder:
rm -r tmp

# Copy all *summary files to temporary folder:
mkdir /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp2

for file in `find /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense \
-name *summary`; do cp $file \
-t /home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp2; \
done

# Transfer count files to laptop:
scp -r alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/Count_summarisation/sense/tmp2 .

# Remove tmp folder:
rm -r tmp2

