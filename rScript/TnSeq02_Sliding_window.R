###
#
# Tn_Seq pipeline 2 : The sliding window method
#
# Author : Fx_stubbe
#
###

# 0) Libraries
# -------------------------------- #

library(tidyverse)
library(data.table)

# 1) Paths and Inputs
# -------------------------------- #

# ----- Parameters
rWindow <- 100
rSliding <- 5

# ----- Paths

#Tn_seq position file
Tn_seq <- "/Users/carlo/Desktop/Eme/Data/Reads/depth_meli16M_2yt_R1.txt"
Tn_data <- fread(paste(Tn_seq,sep=""), sep = "\t", header=F, data.table = F, fill = T, 
                 col.names = c("ChX","Pos", "Insertion"))
#annotation
gff <- read_csv("/Users/carlo/Desktop/Eme/Data/Melitensis16M_gene_table_v3.csv")

#File structure
p.main <- "/Users/carlo/Desktop/Eme/"
p.output <- paste(p.main, "Output/", sep = "")
p.chromo <- paste(p.main, "Output/Chromosomes/", sep = "")
p.table <- paste(p.main, "Output/Tables/", sep = "")
p.rFile <- paste(p.main, "Output/Tables/R", rWindow,"/",sep = "")

# 3) Prepare Main folders
# -------------------------------- #

dir.create(p.output)
dir.create(p.chromo)
dir.create(p.table)
dir.create(p.rFile)

# 4) write files by chromosome
# -------------------------------- #

#Get Ids
ChX_id <- Tn_data$ChX %>% unique()
gff.id <- gff$seqname %>%unique()

#Write split Chromosomes
for(chx in 1:length(ChX_id)){
  Tn_data %>% filter(ChX == ChX_id[chx])  %>%
    write_delim(.,paste(p.chromo, ChX_id[chx],".txt",sep=""), delim = "\t")
}

Tn_data <- Tn_data %>% mutate(log = log10(Insertion+1))
write_delim(Tn_data,paste(p.table, "full_pos.txt",sep=""), delim = "\t")
                                      
# 5) Get essential genes
# -------------------------------- #

# --- Define a splitting function

splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

# --- Deal one chromosome at a time

for(chx in 1:length(ChX_id)){
  
  cat(paste("Dealing with Chromosome : ",gff.id[chx] , sep=""),sep="\n")
  
  #Import Chromosome data
  Data <- fread(paste(p.chromo, gff.id[chx],".txt",sep=""),  sep = "\t", header=T, data.table = F, fill = T)

  #Split into sliding window
  Windows_Pos <- splitWithOverlap(Data$Pos, rWindow, rWindow - rSliding)
  Windows <- splitWithOverlap(Data$Insertion, rWindow, rWindow - rSliding)
  
  #Insertion ratio by sliding window
  Windows_ratio <- sapply(seq(length(Windows)), function(i) {log10(sum(Windows[[i]])+1)})
  
  #Average transposon insertion by sliding window
  average_insertion <- rapply(Windows, sum ) %>% mean()
  
  #-- Get essential genes
  
  #Subset the gff
  gff.t <- gff %>% filter(seqname == gff.id[chx])
  
  #Catching variables
  c_New_Locus <- rep(NA,length(Windows))
  c_Old_Locus <- rep(NA,length(Windows))
  c_Name <- rep(NA,length(Windows))
  c_Product <- rep(NA,length(Windows))
  
  #Windows with no transposon insertion
  target_window <- which(Windows_ratio == 0)
  
  #Loop over gff.t
  for(k in 1:nrow(gff.t)){
    
    #cat(paste(k, "/", nrow(gff.t), sep = ""), sep = "\n")
    
    #Get matching
    matching <- sapply(target_window, function(x){
      if(length(intersect( Windows_Pos[[x]] ,c(gff.t$start[k]:gff.t$end[k]))) > 0){
        target_window <- target_window[-which(target_window == x)]
        return(x)
        }else{return(0)}
    })
    matching <- matching[which(matching != 0)]

    #Match info (Seqname, product, name, ...) in catching variables
    c_New_Locus[matching] <- gff.t$GenBank_locus_tag[k]
    c_Old_Locus[matching] <- gff.t$old_locus_tag[k]
    c_Name[matching] <- gff.t$name[k]
    c_Product[matching] <- gff.t$Product[k]
  }
 
  #Make a dataframe
  Data_filtered <- tibble(New_locus = c_New_Locus, Old_locus = c_Old_Locus, Name = c_Name, Product = c_Product) %>%
    drop_na(.)
  Data_filtered <- Data_filtered %>% group_by(New_locus, Old_locus, Name, Product) %>% tally()
  
  #Write_dataframe
  write_delim(Data_filtered,paste(p.rFile,gff.id[chx],"_R", rWindow ,".txt",sep=""), delim = "\t")
  
}

