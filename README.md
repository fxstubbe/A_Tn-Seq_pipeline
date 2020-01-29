# A Tn-Seq pipeline

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Requirements

**Software requirements**

- Burrows Wheeler Aligner (BWA)
- Samtools
- R (tidyverse, data.table, patchwork)

**Input Data**

- A reference genome (.fasta) and it's annotation file (.gff). Genomes can be downloaded on [NCBI][ncbi], see [Melitensis 16M][Melitensis_16M] for an exemple.
- Samtools
- R (tidyverse, data.table, patchwork)

[ncbi]: https://www.ncbi.nlm.nih.gov/genome/
[Melitensis_16M]: https://www.ncbi.nlm.nih.gov/genome/?term=melitensis%2016M&utm_source=gquery&utm_medium=search


## Read Alignment

To determine the location and abundance of each transposon insertion, short sequencing reads are mapped to a reference genome. This allows to identify the 'correct' genomic loci from which the read originated. 

* **Make a Genomic Index**

Read alignment by BWA uses an FM-index (Full-text index in Minute space) of the reference genome. Both single and paired read libraries are supported. The bwtsw algorithm should be used for short (<100bp) sequencing reads. 

```
bwa index  [-p prefix] [-a algoType] <genome.fasta>
```

Arguments : 

- `-p` (STR) : Prefix of the output database (same as database file name)
- `-a`(STR) : Algorithm for constructing BWT index. Are available **is** and **bwtsw** (recommemded)
- `<genome.fasta>` : Reference genome 

* **Genome Alignment**

The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary

```
bwa mem [-t nThreads] <bwa_index> <Reads> > <output.sam>
```
bwa mem -t 4",  bwa_index, read_f, ">" , str_replace(read_f, ".fastq", ".sam") ,sep = " ") 

### Convert GFF to human readable table

### Make a genome index

### BAM generation using SamTools

### Get coverage by nt and export in .txt format

### Find essential genes using a sliding window strategy

### Graph generation for essential genes
