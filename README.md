# Genome Assembly & Haplotype Phasing Pipeline

> A step-by-step guide for Hi-C-based genome assembly, contamination filtering, and nuclear haplotype separation.

---

## Overview

This pipeline covers:
- Genome assembly from HiFi reads using **hifiasm**
- Quality filtering and contamination removal
- Hi-C contact matrix construction
- Haplotype phasing using **NuclearPhaser**

> ⚠️ NuclearPhaser is designed for **dikaryotic species**, exploiting physical nuclear separation to resolve haplotypes.

> ⚠️ For BLAST searches, ensure you are running on a **computing cluster**.

---

## References & Tools

| Tool | Repository |
|------|-----------|
| hifiasm | https://github.com/chhylp123/hifiasm — [docs](https://hifiasm.readthedocs.io/en/latest/) |
| GenomeAssemblyTools | https://github.com/JanaSperschneider/GenomeAssemblyTools/tree/master |
| NuclearPhaser | https://github.com/JanaSperschneider/NuclearPhaser |
| JupiterPlot (circos) | https://github.com/JustinChu/JupiterPlot |
| RagTag | https://github.com/malonge/RagTag |
| Dgenies (dot plot) | http://dgenies.toulouse.inra.fr/ |

---

## Part I — Assembly

The assembly workflow consists of five steps:
1. hifiasm assembly
2. Low-coverage contig filtering
3. BLAST contamination screening
4. Mitochondrial sequence removal
5. Hi-C contact matrix construction

---

### a) hifiasm

Run hifiasm with Hi-C reads for phased assembly:

```bash
hifiasm -o species.hifi_reads.asm -t 40 \
    --h1 species_R1.fastq.gz \
    --h2 species_R2.fastq.gz \
    species.hifi_reads.fq.gz
```

Map reads back to the assembly with minimap2:

```bash
minimap2 -ax asm20 ${genome} ${longreads} --secondary=no -o mapping.sam
```

Check per-contig coverage using BBMap `pileup.sh`:

```bash
pileup.sh in=mapping.sam out=contig_coverage.txt \
          basecov=basecoverage.txt \
          binsize=1000 bincov=bincoverage.txt
```

---

### b) Low-coverage filtering

Based on the hifiasm coverage plot, identify and remove low-coverage contigs ("noise"). Then combine the two haplotype assemblies:

```bash
cat species.hic.hap1.p_ctg.gfa species.hic.hap2.p_ctg.gfa > combined.gfa
```

Filter contigs by minimum coverage (example threshold: ≥ 2×):

```bash
awk '$3 >= 2' pileup_output.txt > filtered_coverage.txt
```

Extract retained contig sequences:

```bash
seqkit grep -f ListID.txt genome.fa > newGenome.fa
```

---

### c–e) Contamination & mitochondrial filtering

Follow the Sperschneider pipelines for:
- [Collapsed genomic region detection](https://github.com/JanaSperschneider/GenomeAssemblyTools/tree/master/CollapsedGenomicRegions)
- [Contaminant screening](https://github.com/JanaSperschneider/GenomeAssemblyTools/blob/master/ContaminantScreening/README.md)

#### Mitochondrial sequence removal

Build a BLAST database from a mitochondrial reference:

```bash
makeblastdb -in mitochondrion.1.1.genomic.fna -dbtype nucl -out MITO
```

Screen the assembly:

```bash
blastn -query ${genome} -db MITO -dust yes -perc_identity 90.0 \
    -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' \
    > mito.screen.txt
```

Collect mitochondrial contig IDs:

```bash
grep -v '^#' mito.screen.txt >> MITO_CONTIGS.txt
cut -f1 MITO_CONTIGS.txt | sort | uniq -c | sort -nr | head
```

#### NCBI BLAST contamination screening

BLAST remaining contigs against the NCBI nucleotide database (blast+ v2.13.0). Remove contigs whose top 5 hits do not include at least one of: `Puccinia`, `Medioppia`, `Phakopsora`, `Melampsora`, `Uromyces`, or `ribosomal`.

#### Gene mapping (biokanga)

```bash
genome="Pt_Clean_Genome.fasta"
genes="Puctr1_GeneCatalog_transcripts_20131203.nt.fasta"

biokanga index --threads=4 -i ${genome} -o biokanga_index

biokanga blitz --sensitivity=2 --mismatchscore=1 --threads=4 \
    -o CleanGenomeGeneMapping.txt \
    --in=${genes} --sfx=biokanga_index
```

#### BUSCO assessment

```bash
run_BUSCO.py -i ${genome} -o buscov3_clean_assembly \
             -l basidiomycota_odb9 -m geno -sp coprinus -c 4

python busco_converter.py
```

#### Hi-C contact matrix

Convert the Hi-C matrix format if needed (HiCExplorer):

```bash
hicConvertFormat \
    --matrix HiC_MAPQ30.clean_assembly.20000.matrix.h5 \
    --inputFormat h5 \
    -o GInteraction_example \
    --outputFormat GInteractions
```

---

## Part II — Haplotype Phasing

Follow the full NuclearPhaser pipeline: https://github.com/JanaSperschneider/NuclearPhaser

Run NuclearPhaser using the gene mapping, BUSCO, and Hi-C data:

```bash
python NuclearPhaser.py \
    -g Pt_Clean_Genome_GeneMapping.txt \
    -b full_table_buscov3_clean_assembly.tsv \
    -c HiC_MAPQ30.clean_assembly.20000.matrix.tsv \
    -f ${genome} \
    -o /path_to_output_dir/Genome_Phasing_Output/
```

### Visualisation & phase-switch correction

Visualise the two preliminary haplotypes using a dot plot to assess synteny. We recommend **Dgenies** (http://dgenies.toulouse.inra.fr/).

> If phase switches are detected, correct them and **repeat the full NuclearPhaser run** along with all upstream steps (read mapping, etc.).

---

## Part III — Curation

After phasing, inspect the assemblies for residual phase switches and scaffold them against reference genomes.

### a) Circos plot (JupiterPlot)

Use a reference genome and each haplotype to generate a circos plot and visually identify large-scale structural inconsistencies:

```
https://github.com/JustinChu/JupiterPlot
```

### b) Scaffolding with RagTag

Use RagTag to re-scaffold each haplotype against a reference:

```
https://github.com/malonge/RagTag
```
