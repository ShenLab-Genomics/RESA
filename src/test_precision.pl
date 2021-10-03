#!/usr/bin/perl -w

#####################################
#Given 2 vcf.gz, test precision
#12-28-2020
#####################################

$file_scRNA=$ARGV[0];
$file_wes = $ARGV[1];
$temp_dir = "temp_folder_vcf_overlap";
`bcftools isec -p $temp_dir -c none $file_scRNA $file_wes`;
$count_total = `zcat $file_scRNA |grep -v '^#' |wc -l |cut -d " " -f 1`; chomp $count_total;
$f_pos = $temp_dir."/0002.vcf";
$count_pos = `grep -v '^#' $f_pos |wc -l |cut -d " " -f 1`; chomp $count_pos;
print "Total mutation number is $count_total, of which $count_pos are positive\n";
print "Precision is ",sprintf("%.3f",$count_pos/$count_total),"\n";

