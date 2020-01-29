###
#Make a readable table from a gff file (file containing genome features )
###

### 1) load packages ###
# -------------------- #
library(tidyverse)
library(data.table)
library(plyr)

# -------------------- #

### 2) path(s) and inputs ###
# ---------------------- #

f.root <- "/Users/stubbf02/Fx_Stubbe/ressources/genomes/Brucella"
f.genome <- paste( f.root , "/fasta/" , sep = "" )
f.gff <- paste( f.root , "/gff/" , sep = "" )
f.table <- paste( f.root , "/gene_tables/" , sep = "" )

gff.f <- list.files(f.gff)
fasta.f <- list.files(f.genome)
table.s <- gsub( ".fasta" , "_gene_table.csv" , fasta.f )

# ---------------------- #

### 3) Functions ###
# ---------------- #

getAttributeField <- function( x , field , attrsep = ";"){
  s = strsplit( x , split = attrsep , fixed = TRUE )
      sapply( s , function( atts ){
        a = strsplit( atts , split = "=", fixed = TRUE)
        m = match( field , sapply( a , "[" , 1 ))
        if ( !is.na( m ) ){
          rv = a[[ m ]][ 2 ]
        }
        else {
          rv = as.character( NA )
        }
        return( rv )
      } )
}

gffRead <- function( gffFile , nrows = -1 ) {
  cat( "Reading " , gffFile , ": " , sep = "" )
  
  gff <- read.table( gffFile , sep = "\t" , as.is = TRUE , quote = "" ,
                   header = FALSE, comment.char = "#" , nrows = nrows ,
                   colClasses = c( "character" , "character" , "character" , "integer",  
                                "integer", "character" , "character" , "character" , "character" ) )
  
  colnames( gff ) = c( "seqname" , "source" , "feature" , "start" , "end" ,
                    "score" , "strand" , "frame" , "attributes" )
  
  cat( "found" , nrow( gff ) , "rows with classes:" ,
      paste( sapply( gff , class ) , collapse = ", ") , "\n")
  stopifnot( !any( is.na( gff$start ) ) , !any( is.na( gff$end ) ) )
  return( gff )
}


# ---------------- #


### 4) Parse Gff ###
# ---------------- #
for(genome_annotation in gff.f){
  
  # read gff file
  setwd( f.gff )
  gff <- gffRead( gffFile = genome_annotation )
  
  # extract attributes
  gff$GenBank_locus_tag <- getAttributeField( gff$attributes , "locus_tag")
  gff$Product <- getAttributeField(gff$attributes, "product")
  gff$gene_biotype <- getAttributeField(gff$attributes, "gene_biotype")
  gff$protein_id <- getAttributeField(gff$attributes, "protein_id")
  gff$name <- getAttributeField(gff$attributes, "Name")
  
  #remove rows and columns with non-relevant info
  gffClean <- gff[ !gff$feature %in% c( "CDS","region", "sequence_feature" , "STS" , "exon" , "transcript" ) ,  ]
  gffClean <- gffClean %>% filter(GenBank_locus_tag != "NA")
  gffClean <- gffClean[ -c( 6,8,9 ) ]
  
  #add product info to main dataframe
  gffCDS <- gff[ gff$feature %in% c("CDS") , ] %>%
    select(start, end, Product, protein_id)
  
  for(i in 1:nrow(gffCDS)){
    x <- which(gffClean$start == gffCDS$start[i] & gffClean$end == gffCDS$end[i])
    gffClean$Product[x] <- gffCDS$Product[i]
    gffClean$protein_id[x] <- gffCDS$protein_id[i]
  }
  
  #add product info to main dataframe
  gffRNA <- gff[ gff$feature %in% c("tRNA", "ncRNA", "tmRNA", "riboswitch", "SRP_RNA", "rRNA") , ] %>%
    select(start, end, Product)
  
  for(i in 1:nrow(gffRNA)){
    x <- which(gffClean$start == gffRNA$start[i] & gffClean$end == gffRNA$end[i])
    gffClean$Product[x] <- gffRNA$Product[i]
  }

  #output file
  write.csv( gffClean, file = paste( f.table , gsub( ".gff" , "_gene_table.csv" , genome_annotation ) , sep = "" ) , row.names = FALSE , quote = FALSE) 
  
  # ---------------- #
}

  