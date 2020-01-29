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

## The sliding window strategy

### Convert annotation (.gff) to a readble table

Extraction of biotypes (protein_coding,pseudogene, tRNA,rRNA,ncRNA,RNase_P_RNA,mRNA,SRP_RNA) along with their product and  genomic informations.

An exemple gff converted file is provided.

**Script** : gff_converter.R


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

6) For each chromosome, a file containing a list of potentially essential genes (essentiallity index > 0) is created

### Make coverage graphs

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


## The insertion density approach

**Working on it!**

Idea : 
1. Compute for each gene : (read Count)/gene length. Add parameter allowing the exclusion of 5,10,15 or 20 % from both end
2. Plot bimodal distribution and fit model to find pit
3. Classify essentiality based on the peaks

Each gene is assigned a value equal to the count of reads or unique insertions mapped to this gene, divided by the length of the gene.

If the insertions/reads are distributed randomly throughout the genome the density of insertions is expected to be approximately equal for all genes (apart from random variance).

However, the empirical distribution of insertion densities among different genes is generally bimodal, separating non-essential genes that are not affected or only weakly affected by selection from essential or advantageous genes that feature low insertion densities due to selective elimination of the corresponding mutants.

The point separating the two peaks in the distribution can be used as a cutoff value where genes with lower insertion densities are considered essential whereas genes with higher insertion densities are assigned putative non-essential status. Although this approach can be seen as problematic because it compares insertion densities among genes of varying sizes while ignoring the fact that random variance in insertion densities is higher for smaller genes, it has been successfully used in situations where the insertion densities are high, that is, for mutant libraries with a high level of saturation [10,11].

existing methods

10.As the number of insertion sites for any gene is dependent upon the gene length, the values were made comparable by dividing the number of insertion sites by the gene length to give an “insertion index” for each gene. The distribution of insertion indices is bimodal, corresponding to the essential (mode at 0) and nonessential models. For the original mutant pool and each passage condition, we fitted gamma distributions for the two modes using the R MASS library (http://www.r-project.org). Log2-likelihood ratios (LR) were calculated between the essential and nonessential models for each condition and we called a gene essential if it had a log2-LR of less than −2, indicating it was at least four times more likely according to the essential model than the nonessential model. Genes were assigned “nonessential,” if they had a log2-LR of greater than 2. 

11. An analysis similar to that performed for S. typhi, which used transposon density, was utilized (Langridge et al., 2009). There were two complicating factors in using transposon density in this study. First, determination of ORFs in these organisms was performed by auto‐annotation and very often had incorrect start sites. Many genes were seen that had a large number of transposon insertions at the 5′ end, but very few insertions in the bulk of the coding sequence. Second, many essential genes can tolerate insertions in the 3′ portion of the coding region, producing a truncated, but partially functional gene product. To compensate for these issues, 5%, 10%, 15% and 20% of the coding sequence from both the 5′ end and 3′ end of each ORF was disregarded in calculating transposon density. Excluding 20% from each end of genes from the analysis provided the best distinction between essential and non‐essential genes. The transposon density for the internal 60% of each ORF was calculated and plotted on a histogram. While the histogram produced an essential gene peak and a non‐essential gene peak, the area between the peaks was more extensive and contained a large number of unresolved genes. Internal transposon densities were subjected to Ward's clustering analysis to categorize genes as essential, non‐essential and unresolved. 
