#!/usr/bin/perl -w
###############################
# Given a list of SNV vcf files
# Extract the mutation type,
# and 3bp sequence context
# Written by Ning Shen
# 2020-07-10
###############################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("MutationalPatterns", quietly = TRUE))
  BiocManager::install("MutationalPatterns")
if (!requireNamespace("BSgenome", quietly = TRUE))
  BiocManager::install("BSgenome")
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")

library(MutationalPatterns)
library(BSgenome)
library(VariantAnnotation)
library("BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)
library(stringr)


anno = function(vcf_file,sample_name,dir_out){
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
  vcf_file <- c(vcf_file)
  input <- readVcf(vcf_file)
  input_df <- data.frame(row.names(input))
  vcfs = read_vcfs_as_granges(vcf_file, sample_name, ref_genome) 
	type_context = type_context(vcfs[[1]],ref_genome)
	out_matrix = data.frame("Variant"=names(vcfs[[1]]), "Mut_type"=type_context$types, "Seq_context"=type_context$context, "Mut_spec"= paste0(type_context$types, ".", type_context$context))
	row.names(out_matrix)<-names(vcfs[[1]])
	out_matrix<-out_matrix[input_df$row.names.input.,]
	file_out = paste0(dir_out,sample_name[1])
	write.table(out_matrix, file_out, sep="\t", quote=F, row.names=FALSE, col.names=TRUE)
}


