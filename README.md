# RESA: single cell Recurrently Expressed SNV Analysis
Copyright:

## Summary
RESA detects somatic mutation with high precision directly from scRNA-seq data. The key assumption of this method is that cancer cells evolve in a clonal manner and thus expressed somatic mutations have cross-cell recurrence, whereas the noise and artefacts are likely distributed randomly with small probability of recurrence across the cell population. RESA is composed of three main steps: initial variant calling, filtering and labeling, and modeling and refinement.


## RESA workflow
![image](https://user-images.githubusercontent.com/8051136/136513663-8e0f5a8f-29d2-44d2-a7a4-5bed334c3124.png)


## Installation
### Dependencies
Perl (version >=5)

R >= 3.4

R Packages:

Python >=3.7

Python Packages: \
pandas>=1.3.0, rpy2 >=2.9.4, imblearn >= 0.8.0, argparse, sys, os, sklearn, numpy

CTAT-mutation

Minimap2

Strelka

bcftools

annovar

### Install RESA

## Running RESA
### 1. Variant calling

  Please refer to variant calling pipeline instructions for CTAT-mutation and Strelka. Note this step includes alignment and variant calling, and is memory intensive. Below is example code for a single cell sample.
  
      ctat_mutations --left sample1_1.fastq --right sample1_2.fastq --out_dir outdir/sample1 --threads=4
  
      minimap2 -ax splice -t 4 -G 50k -k 21 -w 11 -sr -A 2 -B 8 -O 12,32 -E 2,1 -r 200 -p.5 -N20 -f 1000,5000 -n 2 -m 20 -s40 -g2000 -2K50m -secondary=no GRCh38.p12.genome.fa sample1_1.fastq sample1_2.fastq | samtools view -hb - | samtools sort -o alignment_minimap2/sample1.bam && samtools index alignment_minimap2/sample1.bam
  
      python strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam alignment_minimap2/sample1.bam --referenceFasta GRCh38.p12.genome.fa --rna --runDir sample1_out_path
  
      sample1_out_path/runWorkflow.py -m local -j 4
  

### 2. RESA-filter

      perl filter_and_label_v2.pl  INPUT_folder  RESULTs_folder  file_sample_list  sample_name  minDP  minRecur  maxRecur  PON_filter(optional)

### 3. RESA-refine

      RESA_refine.py -N ./Negative_file -P ./Positive_file -U ./Undefined_file -O dir_out/ -S=True (default: True)

