# A Tn-Seq pipeline

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Requirements

**Software requirements**

- [Burrows Wheeler Aligner][bwa]
- [Samtools][samtools]
- [Qualimap][quali]
- R ([tidyverse][tidy], [data.table][d.table], [patchwork][patch])

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[tidy]: https://www.tidyverse.org/
[d.table]: https://github.com/Rdatatable/data.table
[patch]: https://github.com/thomasp85/patchwork
[quali]: http://qualimap.bioinfo.cipf.es/

**Input Data**

- A reference genome (.fasta) and it's annotation file (.gff). Genomes can be downloaded on [NCBI][ncbi], see [Melitensis 16M][Melitensis_16M] for an exemple.
- Short reads sequencing files (.fastqz)

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

The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary.

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
- `<onput.bam>` : Name of the output BAM file (ex: L001_merged.bam)


**Quality assessment**

Get merged alignment statistics

```
qualimap bamqc [-outdir quality] [-bam <file.bam>] [-gff <annotaion.gff>][-c] [-nt nThread] [--java-mem-size=8G]
```

Arguments : 

- `-outdir` (FOLDER) : Number of threads
- `-bam` (FILE) : Input mapping file in BAM format
- `-gff` (FILE) : Feature file with regions of interest in GFF
- `-nt`(INT) : Number of threads
- `--java-mem-size=8G` : Allows 8 Gb of RAM memory 


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


**Convert annotation (.gff) to a readble table**

Extraction of biotypes (protein_coding,pseudogene, tRNA,rRNA,ncRNA,RNase_P_RNA,mRNA,SRP_RNA) along with their product and  genomic informations.

An exemple gff converted file is provided.

**Script** : gff_converter.R


## The insertion density approach

Each gene is assigned an insertion index value equal to the count of reads or unique insertions mapped to this gene, divided by the length of the gene. Many essential genes can tolerate insertion in the 3'region, therefore producing truncated but functionnal product. 

### Looking for an essential peak

To account for trucnated but functionale product ans misannotated start site; 5%, 10%, 15% and 20% of the coding sequence from both the 5′ end and 3′were disraguarded while calculating the insertion index. The internal plot density was plotted on a histogram, which  produces an essential gene peak and a non‐essential gene peak. The point separating the two peaks in the distribution can be used as a cutoff value where genes with lower insertion densities are considered essential whereas genes with higher insertion densities are assigned putative non-essential status. In our approach, the first histogram on the distribution is defined as the essential peak.

Playing around with binsize and sequence truncation harness the poweer to better define the essential from non essential peaks. 

![](https://user-images.githubusercontent.com/43237088/73437231-41477900-434c-11ea-8b25-4152d100bfed.png)

### Peak definition, a question of saturation

**Low Saturation**

![](https://user-images.githubusercontent.com/43237088/73656045-9a3a4880-468f-11ea-8dea-9fdff7e54a05.png)

**High Saturation**

![](https://user-images.githubusercontent.com/43237088/73462728-8b485300-437c-11ea-84b6-3a3086f345c0.png)

**Sequencing note**

Since hihg saturation Tn-seq results in better essential peak definition, one might think that sequencing a very high number of mutants is preferable. This is true to some extend but has its limits. Indeed, a Tn-seq not saturating enough will fail to resolve low fitness genes, a very high saturation will raise the probability to detect rare variant (Transposon in essential genes) therefore lowering the dectection strength. 

Although this approach can be seen as problematic because it compares insertion densities among genes of varying sizes while ignoring the fact that random variance in insertion densities is higher for smaller genes, it has been successfully used in situations where the insertion densities are high, that is, for mutant libraries with a high level of saturation.

https://genome.cshlp.org/content/19/12/2308

https://onlinelibrary.wiley.com/doi/full/10.1111/mmi.12686


## The sliding window strategy

**Working on it!**

### Compute an essentiallity index

To assess for essentiality of a gene, we used a sliding window approach. Instead of counting insertions in genes, we count insertions in overlapping windows of a fixed size.

**Script** : Sliding_window.R

**User defined parameters**  

- `rWindow` : Size of the sliding window
- `rSliding` : Sliding window shift

**Method** 

1) The coverage file and the converted GFF files are loaded into R.  

2) The coverage file is splitted into chromosomes, which are processed individually.

3) The coverage file is split into coordinates windows matching the `rWindow` and `rSliding` parameters. For exemple, with a rWindow of 200 and a rsliding of 5, a genome of length 3 278 307 bp is split into 655,662 windows. 

4) For each window the sum of aligned reads (coverage) is computed. The logarithm in base 10 is computed for the sum.

5) Each annotated gene in the converted GFF file is assigned an essentiallity index. The index corresponds to the number of empty (0 insertion) window overlaping (even by 1 nt) the annotated gene. If there are no empty window overlapping a given gene (essentiallity index = 0), the gene is skipped.

6) For each chromosome, a file containing a list of potentially essential genes (essentiallity index > 0) is created. Using an R100, 296 essential genes were found.


## Make coverage graphs

**Script** : Automated_Cgrpahs.R
**Script** : Locus_Plots.R

The scripts require the coverage file and the list of potentially essential genes

**User defimed parameters**  

- `boundaries` : Number of nucleotides to represent in the upstream and downstream region of a given gene
- `S_factor` : Smoothing factor (size of average positions to plot) 
- `rWindow` (AUTOMATED): Size of window used to generate the list (used to read the good set of files)
- `Gene` (LOCUS PLOT) : Locus (from converted GFF) to plot 

**Output example** 

![](https://user-images.githubusercontent.com/43237088/73352315-2237e100-4291-11ea-9353-d1c435a9c44b.png)

### Locus identification

Increasing genomic sequence truncation increases essential peak resolution, with a binsize equal to 30. Therefore, the essential peak (centered around 0) comprises genes with a maximum insertion index value of 15. Using the 20% end truncation and an insertion index <= 15, we retrieved 493 essential genes which correspond to 14.63% of the genome. 

![](https://user-images.githubusercontent.com/43237088/73656223-f69d6800-468f-11ea-9bce-e97cfb3dd174.png)
