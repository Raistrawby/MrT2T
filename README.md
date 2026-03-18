# MrT2T


This pipeline is mainly based on this previous work:
- https://github.com/JanaSperschneider/GenomeAssemblyTools/tree/master to create the assembly

- https://github.com/JanaSperschneider/NuclearPhaser : to separate the haplotype
Note that it work on dikariotic species due to the physical separation of the nucleus.  

# Download
Please be sure that you have this tools installed:
/!\ For blast research be sure to be on a computing cluster 

Tools:
- hifiasm: https://github.com/chhylp123/hifiasm , https://hifiasm.readthedocs.io/en/latest/

# I - Assembly

We used hifiasm for genome assemblies and phasing and applied the following filtering steps:
a) hifiasm step
b) Removed low coverage to ensure only high-quality sequences were retained.
c) Filtered BLAST contamination to remove any potential non-target sequences.
d) Mitochondrial sequences were removed to focus on nuclear genome assembly.
e) Created a Hi-C contact matrix to visualize haplotype separation and chromosomal interactions.

## a) hifiasm

```hifiasm -o species.hifi_reads.asm -t 40  ---h1 species_R1.fastq.gz --h2 species_R2.hifi_reads.fq.fastq.gz species.hifi_reads.fq.fastq.gz```
Then map back
```minimap2 -ax asm20 ${genome} ${longreads} --secondary=no -o mapping.sam```
check the coverage
```pileup.sh in=mapping.sam out=contig_coverage2.txt basecov=basecoverage.txt binsize=1000 bincov=bincoverage.txt```
## b) low coverage and cleaning
Based on the plot from hifiasm you can select the "noise" and remove them.
Combined:
species.hic.hap1.p_ctg.gfa
species.hic.hap2.p_ctg.gfa

## c, d, e) based on J. Sperschneider scripts
Then with the combined genomes follow:
https://github.com/JanaSperschneider/GenomeAssemblyTools/tree/master/CollapsedGenomicRegions
https://github.com/JanaSperschneider/GenomeAssemblyTools/blob/master/ContaminantScreening/README.md

At the end you should have a cleaned genome.

```makeblastdb -in mitochondrion.1.1.genomic.fna -dbtype nucl -out MITO ```

```blastn -query $genome -db MITO -dust yes -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > mito.screen.txt```

```cat mito.screen.txt | grep -v "^#" >> MITO_CONTIGS.txt```

```cut -f1 MITO_CONTIGS.txt | sort | uniq -c | sort -nr | head```

```sbatch minimap.sh```
```pileup.sh in=mapping.sam out=contig_coverage.txt```
Use BBMap's pileup.sh tool to calculate read coverage and GC content per contig.

contig >=5X with the awk 

```awk '$3 >= 2' pileup_output.txt > filtered_coverage.txt```


cat genome | seqkit grep -f ListID.txt > newGenome.fa
The remaining contigs were then BLASTed to the NCBI nucleotide database
with blast+ v2.13.0 and contigs without ‘Puccinia’, ‘Medioppia’, ‘Phakopsora’,
‘Melampsora’, ‘Uromyces’, or ‘ribosomal’ in their top five hits were removed

#### hit table of genes
genome="Pt_Clean_Genome.fasta"
genes="Puctr1_GeneCatalog_transcripts_20131203.nt.fasta"

biokanga index --threads=4 -i ${genome} -o biokanga_index -r gene_mapping
biokanga blitz --sensitivity=2 --mismatchscore=1 --threads=4 -o Clean_Genome_GeneMapping.txt --in=${genes} --sfx=biokanga_index

### table of busco hit
run_BUSCO.py -i ${genome} -o buscov3_clean_assembly -l basidiomycota_odb9 -m geno -sp coprinus -c4
## BUSCO converter
```python busco_converter.py```


## Assuming HiCExplorer is installed in ~/programs maybe no need bcse already in h5
hicConvertFormat --matrix ~/programs/HiCExplorer/test/test_data/Li_et_al_2015.h5 \ --inputFormat h5
-o GInteration_example --outputFormat GInteractions


## run nuclearphaser
python NuclearPhaser.py -g Pt_Clean_Genome_GeneMapping.txt -b full_table_buscov3_clean_assembly.tsv \
-c HiC_MAPQ30.clean_assembly.20000.matrix.tsv -f ${genome} -o /path_to_output_dir/Genome_Phasing_Output/

### visualisation 
It is a good idea to visualize the two preliminary haplotypes to see if the dot plot has nice synteny. Here we use Dgenies (http://dgenies.toulouse.inra.fr/) for visualization.

## At this stage phase switch correction
and odo the nuclear phaser as weel as the previous analysies again (bowtie etc)

# II - Phased the assembly
Follow this pipeline: https://github.com/JanaSperschneider/NuclearPhaser

# III - Curation

At the end you can detect phase swithes:

### a) circos plot

https://github.com/JustinChu/JupiterPlot

use reference genomes and each haplotype to perform a circos plot

Then use ragtag to separate and scaffold again
https://github.com/malonge/RagTag
