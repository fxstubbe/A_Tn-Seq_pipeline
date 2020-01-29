# A Tn-Seq pipeline

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Requirements

**Software requirements**

- Burrows Wheeler Aligner ([BWA][bwa])
- [Samtools][samtools]
- R ([tidyverse][tidy], [data.table][d.table], [patchwork][patch])

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[tidy]: https://www.tidyverse.org/
[d.table]: https://github.com/Rdatatable/data.table
[patch]: https://github.com/thomasp85/patchwork

**Input Data**

- A reference genome (.fasta) and it's annotation file (.gff). Genomes can be downloaded on [NCBI][ncbi], see [Melitensis 16M][Melitensis_16M] for an exemple.
- Samtools
- R (tidyverse, data.table, patchwork)

[ncbi]: https://www.ncbi.nlm.nih.gov/genome/
[Melitensis_16M]: https://www.ncbi.nlm.nih.gov/genome/?term=melitensis%2016M&utm_source=gquery&utm_medium=search


## Read Alignment

To determine the location and abundance of each transposon insertion, short sequencing reads are mapped to a reference genome. This allows to identify the 'correct' genomic loci from which the read originated. 

### Make a Genomic Index

Read alignment by BWA uses an FM-index (Full-text index in Minute space) of the reference genome. Both single and paired read libraries are supported. The bwtsw algorithm should be used for short (<100bp) sequencing reads. 

```
bwa index  [-p prefix] [-a algoType] <genome.fasta>
```

Arguments : 

- `-p` (STR) : Prefix of the output database (same as database file name)
- `-a`(STR) : Algorithm for constructing BWT index. Are available **is** and **bwtsw** (recommemded)
- `<genome.fasta>` : Reference genome 

### Genome Alignment

The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary

**Single-end**

```
bwa mem [-t nThreads] <db.prefix> <reads.fq> > <output.sam>
```

**Paired-end**

```
bwa mem [-t nThreads] <db.prefix> <reads_F1.fq> <reads_F2.fq> > <output.sam>
```

Arguments : 

- `-t` (INT) : Number of threads
- `<db.prefix>` : Index database outputed at the previous step (**Make a Genomic Index**) 
- `<Reads.fq>` : Sequencing Reads (F1, F2 if paired-end)
- `<output.sam>` : Name of the output alignment file (ex: L001_R2.sam) 

For additionnal informations, see [BWA Manual][bwa_manual].

[bwa_manual]: http://bio-bwa.sourceforge.net/bwa.shtml

### Location and Abundance of transposon insertions 

**Sort SAM files**

Sort alignments by leftmost coordinates 

```
samtools sort [-@ nThreads] [-o Output file] <input.sam>
```

Arguments : 

- `-@` (INT) : Number of threads
- `-o` (STR) : Name of the output (ex: L001_sorted_R2.bam)
- `<input.sam>` : SAM file to sort (Generated at the previous step **Genome Alignment**)


**Merge BAM files**

Merge aligned (**same** genome) and sorted reads into a single alignment file

```
samtools merge [-@ nThreads] [-b list.txt] <output.bam>
```

Arguments : 

- `-@` (INT) : Number of threads
- `-b` (FILE) : List containing path to the sorted BAM files to merge
- `<onput.sam>` : Name of the output BAM file (ex: L001_merged.bam)


**Get coverage**

```
samtools depth [-a] [-f merged_file.bam] [-o <FILE.txt>]
```

Arguments : 
- `-a` : Output all positions (including those with zero depth)
- `-f` (FILE) : Use the BAM files specified in the FILE (a file of filenames, one file per line)
- `-o` (STR) : Name of the output file (ex: L001_coverage.txt)

Exemple of Output table : 

Chromosome      | Coordinates                        | Coverage
---            | ---                            | ---
I          | 1                     |  0
I          | 2                    |  0
I          | 3                    |  117
...          | ...                   |  ...
II          | 1917                    |  98
II          | 1918                    |  98

## The sliding window strategy
