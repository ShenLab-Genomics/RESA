#!/usr/bin/perl -w
################################
# Read all SNV info, print 1 vcf per group of samples
# Only print uniq var per group
# Written by Ning Shen
# 2020-10-08 (unfinished)
###############################

#my $dir_var = "/home/ns362/park_lab/Ning/GEO/Wang_cPCs_MM_FEBS_2020_GSE137545/results/RNA-seq/";
my $dir_var = "/home/ns362/park_lab/Ning/GEO/Wang_cPCs_MM_FEBS_2020_GSE137545/results/RNA-seq_v2/";
my $dir_res = $ARGV[0]; #"../results_analysis/variants/vcf_filter_v6/";
my $dir_somatic = $ARGV[1]; # "../results_analysis/variants/vcf_filter_v6/temp_folder_filter_dp5_recur5/";
my $dir_out = $ARGV[2]; #"../results_analysis/variants/vcf_filter_v6/temp_folder_filter_dp5_recur5_vcf/";
my $mut_type = $ARGV[3]; #somatic or germline

my $list = $dir_var."list_samples.txt";
my $dp = 5;
my %saw_gr=();
my %map_sample=();
open(L,"$list")||die;
while(<L>){
	chomp;
	if($_=~/(SRR\d+)_(\S+)/){
		$sample = $1."_".$2; $group = $2;
		$map_sample{$sample}=$group;
	}else{
		#die "Pattern matching problem for identifying sample id and patient id\n";}
		next;
	}
	$dir_isec = $dir_res."filter_".$sample."_dp$dp";
	unless(-e $dir_isec){next;}
	$vcf_input = $dir_isec."/0000.vcf";
	$file_out_v1 = $dir_somatic.$sample.".txt";
	unless(-e $file_out_v1){
		print "Cannot find file $file_out_v1\n";
		next;
	}
	$count_mut = `wc -l $file_out_v1 | cut -d " " -f 1`; chomp $count_mut;
	if($count_mut == 0){
		print "No mutation found for file $file_out_v1\n";
		next;
	}
	$file_vcf_out = $dir_out."SNV_merged_".$group.".vcf";
	unless(-e $file_vcf_out){
		$gr=$group;
		$saw_gr{$gr}=1;
		%$gr=();
#		%saw_var=(); %dup_var=();
		`grep '^#' $vcf_input > $file_vcf_out`;
	}
#	%saw_var=(); $var="";
	open(F,"$file_out_v1")||die;
	while(<F>){
		my (@arr)=split /\s+/,$_;
		if($mut_type eq "germline"){
			$var = $arr[2]."\t".$arr[3]."\t".$arr[5]."\t".$arr[6]."\n"; ##For germline and rnaedit variants
		}elsif($mut_type eq "somatic"){
			$var = $arr[0]."\t".$arr[1]."\t".$arr[3]."\t".$arr[4]."\n";
		}else{
			die "Specify mutation type: germline or somatic!\n";
		}
#		if(!defined $dup_var{$var}){
#			$saw_var{$var}=1; $dup_var{$var}=1;
#		}
		if(!defined $$gr{$var}){
			$$gr{$var}=1;
		}
		next;
	}
	close(F);
	next;
}

foreach $group(keys %saw_gr{$group}){
	my %var_print=();
	$file_vcf_out = $dir_out."SNV_merged_".$group.".vcf";
	open(O,">>$file_vcf_out")||die;
	foreach $sample(keys {values %map_sample eq $group}){
		$dir_isec = $dir_res."filter_".$sample."_dp$dp";
		unless(-e $dir_isec){next;}
		$vcf_input = $dir_isec."/0000.vcf";
		
		open(V,"$vcf_input")||die;
		while(<V>){
			if($_=~/^#/){
#				print O $_;
				next;
			}else{
				my (@arr)=split /\t/,$_;
				my $var = $arr[0]."\t".$arr[1]."\t".$arr[3]."\t".$arr[4]."\n";
#			if(defined $saw_var{$var}){
				if(defined $$saw_gr_var{$var}){
					if(!defined $var_print{$var}){
						print O $_;
						$var_print{$var}=1;
					}
				}
			}
			next;
		}
		close(V);
	}
#	undef %saw_var;
	next;
}
		




