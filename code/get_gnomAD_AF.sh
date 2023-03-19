#!/bin/bash

# this script meant to be used in slurm array 0-23
# the script requires accsess to AWS service (connected user, no need to pay)
# I think we might be able to change it to the vcf files on google cloud, no need a user to ge them

# function which exit the script if the command above did not preformed well
check_exit_code() {
  if [ $? -ne 0 ]; then # $? refferes to the previous exit code, 0 means it ended well
    echo Last command, "$1" failed. Exiting.
    exit 1
  fi
}


# loading the neccesery modules
module load hurcs
check_exit_code "module load"
module load bcftools
check_exit_code "module load"


CHROM_LIST=(chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY)
CHROM=${CHROM_LIST[$SLURM_ARRAY_TASK_ID]} # $SLURM_ARRAY_TASK_ID is the number which slurm gave to the current job in the array

echo "start bcftools for $CHROM"
# bcftools command flags explanation-
# o - the path to the output file
# i - include ( filter only variant with FILTER == PASS
# R - regions file, takes variants that are on the regions described in the given file
# f - the output format 
# last variable is the location of the vcf file. We are using the file which is stored on amazon cloud
bcftools query -o data/${CHROM}_gnomAD_AF_19_03.tsv -R data/peak_bcftools_format.tsv -f "%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AF_popmax\t%FILTER\n" https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${CHROM}.vcf.bgz

check_exit_code "module bcftools"
echo "DONE $CHROM"