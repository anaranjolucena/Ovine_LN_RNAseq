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


# I put file with adapters in this directory from my terminal (not from rodeo)
scp -r Illumina_PE_adapters.txt alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/filt_fastq

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 151 \
-pe1 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
-pe2 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz \
-o /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

################################################################
# FastQC quality check of filtered FASTQ files from 1st sample #
################################################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering
cd !$


# Run FastQC in one file to see if it's working well:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R1_001.fastq.gz

#In new rodeo tab

fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R2_001.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line)
# and check HTML reports:
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-filtering/*fastqc.zip .

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

#Create directory for alignment of one file with 1st filter conditions:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/1st_filtered
cd !$

# Mapping reads from one FASTQ file to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R1_001.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/N12_S29_L002/trimmed_N12_S29_L002_R2_001.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29 \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

#Create directory in quality_check for FASTQC of BAM files
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/1st_filtered/N12_S29Aligned.out.bam

# Transfer compressed folder to personal laptop via SCP (tab in personal laptop terminal)
# and check the HTML report:
scp alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-alignment_BAM/N12_S29Aligned.out_fastqc.zip .

################################################################################
# Adapter-contamination and quality filtering of raw FASTQ file 2nd conditions #
################################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html


# Run ngsShoRT in one pair of reads changing number of bp covered to 140 (furthest matching index)
# First, navigate to filt_fastq folder, where adapter sequences are saved (or give path):
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 151 \
-pe1 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
-pe2 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz \
-o /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/2nd_filtered \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 140 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

########################################################################
# FastQC quality check of option 2 filtering FASTQ files from 1 sample #
########################################################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_filtering
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/2nd_filtered/trimmed_N12_S29_L002_R1_001.fastq.gz

#In new rodeo tab
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/2nd_filtered/trimmed_N12_S29_L002_R2_001.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line) and check HTML reports:
# Navigate locally to quality check directory and create new directory:
mkdir post-2nd_filtering
cd !$
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_filtering/*fastqc.zip .


###########################################################################################
# Alignment of 2nd filtered FASTQ files against the Ovis aries reference genome with STAR #
###########################################################################################


#Create directory for alignment of one file with 2nd filter conditions:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/2nd_filtered
cd !$

# Mapping reads from one FASTQ file to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/2nd_filtered/trimmed_N12_S29_L002_R1_001.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/2nd_filtered/trimmed_N12_S29_L002_R2_001.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29 \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

#Create directory in quality_check for FASTQC of BAM files
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_alignment_BAM
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_alignment_BAM --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/2nd_filtered/N12_S29Aligned.out.bam

# Transfer compressed folder to personal laptop via SCP (tab in personal laptop terminal)
# and check the HTML report:
scp alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-2nd_alignment_BAM/N12_S29Aligned.out_fastqc.zip .

#############################################################################################
# 3rd- 4rth Adapter-contamination and quality filtering of raw FASTQ file SOFTWARE:cutadapt #
#############################################################################################

# Required software is cutadapt (version 2.9). 
# https://cutadapt.readthedocs.io/en/stable/guide.html#

#Create output directory
mkdir /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/3rd_filtered
cd !$

#Run cutadapt for 100 bp minimum read legnth
cutadapt -q 20 -m 100 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o N12_S29_L002_R1_trimmed.fastq.gz \
-p N12_S29_L002_R2_trimmed.fastq.gz /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
/home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz


#Create output directory
mkdir /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/4rth_filtered
cd !$

#Run cutadapt for 120 bp minimum read legnth (-j to include number of cores...)
#Remember to do nohup to obtain the summary in a text file 
cutadapt -q 20 -m 120 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 20 -o N12_S29_L002_R1_trimmed.fastq.gz \
-p N12_S29_L002_R2_trimmed.fastq.gz /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
/home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz


##############################################################################
# FastQC quality check of option 3 and 4 filtering FASTQ files from 1 sample #
##############################################################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory 3rd filtering:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-3rd_filtering
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-3rd_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/3rd_filtered/N12_S29_L002_R1_trimmed.fastq.gz

#In new rodeo tab
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-3rd_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/3rd_filtered/N12_S29_L002_R2_trimmed.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line) and check HTML reports:
# Navigate locally to quality check directory and create new directory:
mkdir post-3rd_filtering
cd !$
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-3rd_filtering/*fastqc.zip .


# Create and enter the quality check output directory 4th filtering:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-4th_filtering
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-4th_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/4rth_filtered/N12_S29_L002_R1_trimmed.fastq.gz

#In new rodeo tab
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-4th_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/4rth_filtered/N12_S29_L002_R2_trimmed.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line) and check HTML reports:
# Navigate locally to quality check directory and create new directory:
mkdir post-4th_filtering
cd !$
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-4th_filtering/*fastqc.zip .


###################################################################################################
# Alignment of 3rd and 4th filtered FASTQ files against the Ovis aries reference genome with STAR #
###################################################################################################


#Create directory for alignment of one file with 3rd filter conditions:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/3rd_filtered
cd !$

# Mapping reads from one FASTQ file (3rd filtered option) to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/3rd_filtered/N12_S29_L002_R1_trimmed.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/3rd_filtered/N12_S29_L002_R2_trimmed.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29 \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &


#Create directory for alignment of one file with 4th filter conditions:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/4th_filtered
cd !$

# Mapping reads from one FASTQ file (4th filtered option) to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/4rth_filtered/N12_S29_L002_R1_trimmed.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/4rth_filtered/N12_S29_L002_R2_trimmed.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29 \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

#############################################################################################
# 5th Adapter-contamination and quality filtering of raw FASTQ file SOFTWARE:ngsShoRT #
#############################################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html


# I put file with adapters in this directory from my terminal (not from rodeo)
scp -r Illumina_PE_adapters13bp.txt alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/filt_fastq

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 151 \
-pe1 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R1_001.fastq.gz \
-pe2 /home/workspace/ccorreia/ovineRNAseq/fastq/Ovine/N12_S29_L002_R2_001.fastq.gz \
-o /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/5th_filtered \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters13bp.txt -5a_mp 100 -5a_del 0 \
-5a_ins 0 -5a_fmi 140 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

########################################################################
# FastQC quality check of option 5 filtering FASTQ files from 1 sample #
########################################################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-5th_filtering
cd !$

# Run FastQC:
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-5th_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/5th_filtered/trimmed_N12_S29_L002_R1_001.fastq.gz

#In new rodeo tab
fastqc -o /home/workspace/alucena/ovineLN_RNAseq/quality_check/post-5th_filtering --noextract --nogroup \
-t 20 /home/workspace/alucena/ovineLN_RNAseq/filt_fastq/5th_filtered/trimmed_N12_S29_L002_R2_001.fastq.gz

# Transfer compressed folders to personal laptop via SCP (in a new tab from your own mac command line) and check HTML reports:
# Navigate locally to quality check directory and create new directory:
mkdir post-5th_filtering
cd !$
scp -r \
alucena@rodeo.ucd.ie:/home/workspace/alucena/ovineLN_RNAseq/quality_check/post-5th_filtering/*fastqc.zip .

###########################################################################################
# Alignment of 5th filtered FASTQ files against the Ovis aries reference genome with STAR #
###########################################################################################


#Create directory for alignment of one file with 5th filter conditions:
mkdir /home/workspace/alucena/ovineLN_RNAseq/STAR-2.7.3a_alignment/5th_filtered
cd !$

# Mapping reads from one FASTQ file (5th filtered option) to the indexed genome,to check it works :
nohup STAR --runMode alignReads --runThreadN 30 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_NCBI/STAR-2.7.3a_index_150 \
--readFilesIn \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/5th_filtered/trimmed_N12_S29_L002_R1_001.fastq.gz \
/home/workspace/alucena/ovineLN_RNAseq/filt_fastq/5th_filtered/trimmed_N12_S29_L002_R2_001.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./N12_S29_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &


