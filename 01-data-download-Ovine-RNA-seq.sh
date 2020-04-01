########################################
# Ovine Lymph Node Liver Fluke RNA-seq #
########################################

# Author: Carolina N. Correia
# Last updated on: 18/03/2020

##############################
# Download raw data from BYU #
##############################

# Create and enter working directory on Rodeo:
mkdir -p /home/workspace/ccorreia/ovineRNAseq/fastq
cd !$

# Download files:
screen -D -R ovine_download
wget -r --user=dnasc_user --ask-password --no-parent --no-check-certificate https://files.rc.byu.edu/fslg_dnasc/data_transfer/Ovine/
# Detach the screen session by pressing Ctrl+A then d
# Terminate the ssh session if desired
exit

# Check if download is complete:
screen -D -R ovine_download # Detach the screen session by pressing Ctrl+A then d.

# After download has been completed:
screen -D -R ovine_download
exit
screen -X -S ovine_download quit

# Modify file permissions to read and execute only:
chmod -R 555 *.fastq.gz

##############
# Next steps #
##############

# Please go to file: 02-QC-Ovine-RNA-seq.sh
