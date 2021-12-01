# RESA: single cell Recurrently Expressed SNV Analysis
Copyright: MIT license.

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

bcftools

annovar

Dependencies for step1: CTAT-mutation; Minimap2; Strelka

### Install RESA

  Before installing RESA, please registe on https://www.openbioinformatics.org/annovar/annovar_download_form.php to download ANNOVAR. 
  To download and install RESA, use the commands.
  
      docker pull tianyunz/resa:latest
      
      docker run -it tianyunz/resa:latest /bin/bash
      
  The following steps are installing the ANNOVAR and downloading the default database used in the RESA package. The size of database gnomad30_genome is large,  and we suggest setting up the docker image size greater than 90G. 
      
      docker cp annovar.latest.tar.gz containerID:/bin/annovar.latest.tar.gz
      
      docker exec -it containerID /bin/bash
      
      tar -zxvf /bin/annovar.latest.tar.gz
      
      cd /bin/annovar/
      
      ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/
      
	    ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/
      
	    ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/

  hg38_RNAedit.txt cannot be downloaded directly from the ANNOVAR, and here shows one way to convert it from vcf file to the ANNOVAR readable database.
  
	  wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/GRCh38.RNAediting.vcf.gz
    
	  sed -i '1,5d' GRCh38.RNAediting.vcf
    
	  awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}' GRCh38.RNAediting.vcf > temp.txt 
    
	  sed '1i\#Chr\tStart\tEnd\tRef\tAlt' temp.txt > hg38_RNAedit.txt
  
## Running RESA
### 1. Variant calling

  This step applies standard reference alignment and variant calling pipelines developed by others. For each cell, two aligners and variant callers need to be applied. This step is memory intensive and time consuming. To allow more flexibility to make best use of the resource available on the user's side, this step is not incoporated in the RESA package. Please refer to variant calling pipeline instructions for CTAT-mutation and Strelka. 
  
  Here are example codes to run for a single cell sample:
  
      ctat_mutations --left sample1_1.fastq --right sample1_2.fastq --out_dir outdir/sample1 --threads=4
  
      minimap2 -ax splice -t 4 -G 50k -k 21 -w 11 -sr -A 2 -B 8 -O 12,32 -E 2,1 -r 200 -p.5 -N20 -f 1000,5000 -n 2 -m 20 -s40 -g2000 -2K50m -secondary=no GRCh38.p12.genome.fa sample1_1.fastq sample1_2.fastq | samtools view -hb - | samtools sort -o alignment_minimap2/sample1.bam && samtools index alignment_minimap2/sample1.bam
  
      python strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam alignment_minimap2/sample1.bam --referenceFasta GRCh38.p12.genome.fa --rna --runDir sample1_out_path
  
      sample1_out_path/runWorkflow.py -m local -j 4
  

### 2. RESA-filter

      perl filter_and_label.pl  INPUT_folder  RESULTs_folder  file_sample_list  sample_name  minDP  minRecur  maxRecur  PON_filter(0 or 1) PON_file(if PON_filter is 1)
      
  Here are the example codes using RESA-filter to filter example dataset.
      
      mkdir /bin/RESA/examples/result/
      
      perl /bin/RESA/src/filter_and_label_v2.pl \
	    /bin/RESA/examples/input/K562_primers/ \
	    /bin/RESA/examples/result/ \
	    /bin/RESA/examples/input/K562_primers/list_samples_K562_primers.txt K562_primers 3 3 42

### 3. RESA-refine

      RESA_refine.py -N ./Negative_file -P ./Positive_file -U ./Undefined_file -O dir_out/ -S=True (default: True)
  
  Here are the example codes using RESA-refine to predict mutations of the example dataset. And predicted mutations have saved in prdicted_value.csv

      python3.8 RESA_refine.py \
	    -N /bin/RESA/examples/result/K562_primers/folder_confident/Double_neg_qualrecur_fail.txt \
	    -P /bin/RESA/examples/result/K562_primers/folder_confident/Double_pos_qualrecur_pass.txt \
	    -U /bin/RESA/examples/result/K562_primers/folder_unsure/Undefined_B.txt \
	    -O /bin/RESA/examples/result/K562_primers/

