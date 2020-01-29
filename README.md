# A Tn-Seq pipeline

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Requirements

**Software requirements**

- Burrows Wheeler Aligner (BWA)
- Samtools
- R (tidyverse, data.table, patchwork)

**Input Data**

- A reference genome (.fasta) and it's annotation file (.gff)

[Melitensis 16M]: https://www.ncbi.nlm.nih.gov/genome/?term=melitensis%2016M&utm_source=gquery&utm_medium=search

- Samtools
- R (tidyverse, data.table, patchwork)

## Align sequencing reads

To determine the location and abundance of each transposon insertion, short sequencing reads are mapped to a reference genome. This allows to identify the 'correct' genomic loci from which the read originated. 


**Choose the Reference Genome **

* **Retrieve Genome and Make index**



```
bwa index -a bwtsw <genome.fasta>
```

* **Retrieve Genome annoatation and convert to readable table**


### Convert GFF to human readable table

### Make a genome index

### BAM generation using SamTools

### Get coverage by nt and export in .txt format

### Find essential genes using a sliding window strategy

### Graph generation for essential genes
