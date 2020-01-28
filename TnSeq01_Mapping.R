###
#
# Tn_Seq pipeline 1 : Align to Genome 
#
# Author : Fx_stubbe
#
###


# 0) Library
# -------- -------- -------- -------- #

library(tidyverse)

# 1) Paths and Inputs
# -------- -------- -------- -------- #

#Folder where the read pair for the giving experiment is
main.p <- "/Users/carlo/Desktop/Eme/Data/Reads/"

#genome for the index
bwa_index <- "/Users/carlo/Desktop/Eme/Data/Genome/Melitensis16M.fasta"

#Give a name to the experiement
name_exp <- "meli16M_2yt_R1"


# 2) Make a genomic index
# If index is already made --> SKIP this step
# -------- -------- -------- -------- #

cmd <- paste("/Users/carlo/Documents/TOOLS/bwa/bwa index -a bwtsw", bwa_index , sep = " ")
system(cmd)

# 3) Genome mapping
# Primer based sequencing so no adapter treaming needed
# -------- -------- -------- -------- #

setwd(main.p)
reads <- list.files()

# Align each set of reads onto the reference genome
for(read_f in reads){
  cmd <- paste("/Users/carlo/Documents/TOOLS/bwa/bwa mem -t 4",  bwa_index, read_f, ">" , str_replace(read_f, ".fastq", ".sam") ,sep = " ") 
  system(cmd)
}

# Sort files
sam_id <- list.files(pattern = ".sam") 
bam_list <- "list_bams.txt"
for(sam_f in sam_id){
  cmd <- paste("samtools sort -@ 4 -o sorted_", str_replace(sam_f, ".sam", ".bam")," ", sam_f, sep = "")
  cat(paste(main.p,"sorted_",str_replace(sam_f, ".sam", ".bam"), sep= ""), sep = "\n", file = bam_list, append = T)
  system(cmd)
}

# Merge Bam files
final_bam <- "final_bam.txt"
cmd <- paste("samtools merge -@ 4 -b list_bams.txt ", name_exp ,".bam", sep = "")
cat(paste(main.p, name_exp ,".bam", sep= ""), sep = "\n", file = final_bam, append = T)
system(cmd)

# Get coverage nt by nt
cmd <-  paste("samtools depth -a -f final_bam.txt > depth_",name_exp,".txt",sep = "")
#cmd <- paste("bamtools coverage -in ", name_exp ,".bam ", "-out coverage_" ,name_exp,".txt", sep ="")
system(cmd)
