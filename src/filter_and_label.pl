#!/usr/bin/perl -w

#########################################################################
# This code takes inital gatk called vcf.gz file as input, we recommend 
# using ctat-mutation pipeline
# The initial variant set goes through a series of filtering, including 
# focusing on exonic region, remove germline & RNA-editing sites, overlap
# with 2nd aligner+variant caller result, quality and recurrence based
# filter.
# Throughout the filter, the variants are labeled with three categories:
# Confident (Pos, Neg), Unsure, Remove
#
# Written by Ning Shen
# 2021-03-19
#########################################################################

my $dir_var = ""; #Folder where inital gatk variant call and 2nd variant caller vcf files are stored
my $dir_temp = ""; #Folder for temporary files
my $dir_out = ""; #Folder for output

my $file_list = ""; #list of sample names

my $ref_genome = "hg38";
my $dp_cut = 3;

my $dir_inter = $dir_temp."folder_int_files";
my $dir_conf = $dir_temp."folder_confident";
my $dir_uns = $dir_temp."folder_unsure";
my $dir_rm = $dir_temp."folder_remove";
unless(-e $dir_inter){`mkdir $dir_inter`;}
unless(-e $dir_conf){`mkdir $dir_conf`;}
unless(-e $dir_uns){`mkdir $dir_uns`;}
unless(-e $dir_rm){`mkdir $dir_rm`;}
my $sample_name; 
my $file_gatk; my $file_strelka;
my $file_annovar_input; my $file_annovar_gene_header; my $file_annovar_gene; my $file_annovar_exonic; my $file_annovar_gnomad_fil_header; my $file_annovar_gnomad_fil; my $file_annovar_dbsnp_fil_header; my $file_annovar_dbsnp_fil; my $file_annovar_rnaedit_header; my $file_annovar_rnaedit_fil; 
open(L1,"$file_list")||die;
while(<L1>){
	chomp;
	$sample_name=$_;
#	$file_ctat_gz = $dir_var."mutation_calls_ctat/".$sample."/work_dir/HaplotypeCaller.raw_variants.vcf_snpeff.vcf_adj.vcf_dbsnp.vcf_RNAedit.vcf_PASSreads.vcf.gz";
#	$file_ctat_gz = $dir_var."mutation_calls_ctat/".$sample_name."/work_dir/HaplotypeCaller.raw_variants.vcf.gz";
	$file_gatk = $dir_var."variant_call_gatk/".$sample_name.".vcf.gz";
	unless(-e $file_gatk){
                print "Cannot find variant vcf file $file_gatk\n";next;
        }
	$file_strelka = $dir_var."variant_call_strelka/".$sample_name.".vcf.gz";
	$file_var_strelka = $dir_inter."/Strelka_SNVs_".$sample_name;
	`zcat $file_strelka | grep -v '^#' |grep '^chr' |grep -v '^chrM' |cut -f 1,2,3,4,5 > $file_var_strelka`;
	my %snv_strelka = ();
	open(FVS,"$file_var_strelka")||die;
	while(<FVS>){
		my @arr = split /\s+/,$_;
		unless($arr[3] eq "A" || $arr[3] eq "T" || $arr[3] eq "C" || $arr[3] eq "G"){next;}
		unless($arr[4] eq "A" || $arr[4] eq "T" || $arr[4] eq "C" || $arr[4] eq "G"){next;}
		my $var_info = $arr[0]."_".$arr[1]."_".$arr[3]."_".$arr[4];
		$snv_strelka{$var_info}=1;
		next;
	}

	$file_annovar_input = $dir_inter."/Input_annovar_format_".$sample_name;
	`zcat $file_gatk | /home/ns362/park_lab/Ning/software/annovar/convert2annovar.pl -format vcf4 - -includeinfo | grep '^chr' |grep -v '^chrM' > $file_annovar_input`; ##convert variants to annovar format, exclude mutations in chrM or "GL" chrs
	##Map variant to gene and filter exonic variant
	$file_annovar_gene_header = $dir_inter."/Input_annovar_ensGeneAnnot_".$sample_name;
	$file_annovar_gene = $file_annovar_gene_header.".variant_function";
	`/home/ns362/park_lab/Ning/software/annovar/annotate_variation.pl -out $file_genloc_ann_header -exonicsplicing -build hg38 $file_scan /home/ns362/park_lab/Ning/software/annovar/humandb/ -dbtype ensGene`; ##Annotate variant gene location
	my $temp_loc = `grep '^exonic' $file_annovar_gene`;
        my (@arr_temp_loc) = split /\n/,$temp_loc;
	my %snv_exonic=(); 
        foreach my $line(@arr_temp_loc){
        	my @arr = split /\t/,$line;
                unless($arr[5] eq "A" || $arr[5] eq "T" || $arr[5] eq "C" || $arr[5] eq "G"){next;}
                unless($arr[6] eq "A" || $arr[6] eq "T" || $arr[6] eq "C" || $arr[6] eq "G"){next;}
        	my $var_info = $arr[2]."_".$arr[3]."_".$arr[5]."_".$arr[6];
		$snv_exonic{$var_info}=1;
	}	
	$file_annovar_exonic = $dir_inter."/Filter_annovar_exonic_SNV_".$sample_name;
	open(FE,">$file_annovar_exonic")||die;
	open(FI,"$file_annovar_input")||die;
	while(<FI>){
		my (@arr)=split /\t/,$_;
		$var_info = $arr[0]."_".$arr[1]."_".$arr[3]."_".$arr[4];
		if(defined $snv_exonic{$var_info}){
			print FE $_;
		}
		next;
	}
	close(FI);
	close(FE);
	##Filter against common SNP and RNA editing sites
	$file_annovar_gnomad_fil_header = $dir_inter."/Filter_annovar_gnomad_fil_".$sample_name;
	$file_annovar_gnomad_fil = $dir_inter."/Filter_annovar_gnomad_fil_".$sample_name.".hg38_gnomad30_genome_filtered";
	`/home/ns362/park_lab/Ning/software/annovar/annotate_variation.pl -filter -dbtype gnomad30_genome -buildver hg38 -out $file_annovar_gnomad_fil_header $file_annovar_exonic /home/ns362/park_lab/Ning/software/annovar/humandb/`;
	$file_annovar_dbsnp_fil_header = $dir_inter."/Filter_annovar_dpsnp_fil_".$sample_name;
	$file_annovar_dbsnp_fil = $dir_inter."/Filter_annovar_dpsnp_fil_".$sample_name.".hg38_dbsnp_filtered";
	`/home/ns362/park_lab/Ning/software/annovar/annotate_variation.pl -filter -dbtype avsnp147 -buildver hg38 -out $file_annovar_dbsnp_fil_header $file_annovar_gnomad_fil /home/ns362/park_lab/Ning/software/annovar/humandb/`;
	$file_annovar_rnaedit_header = $dir_inter."/Filter_annovar_rnaedit_fil_".$sample_name;
	$file_annovar_rnaedit_fil = $dir_inter."/Filter_annovar_rnaedit_fil_".$sample_name.".hg38_generic_filtered";
	`/home/ns362/park_lab/Ning/software/annovar/annotate_variation.pl -filter -dbtype generic -genericdbfile hg38_RNAedit.txt -buildver hg38 -out $file_annovar_dbsnp_header $file_annovar_dbsnp_fil /home/ns362/park_lab/Ning/software/annovar/humandb/`;

	##Test for overlap with strelka (2nd) variant caller
	$file_double_call_pos = $dir_inter."/Double_call_pos_".$sample_name;
	$file_double_call_neg = $dir_inter."/Double_call_neg_".$sample_name;
	$file_double_pos_qualrecur_pass = $dir_inter."/Double_pos_qualrecur_pass_".$sample_name;
	$file_double_pos_qualrecur_fail = $dir_inter."/Double_pos_qualrecur_fail_".$sample_name;
	$file_double_neg_qualrecur_pass = $dir_inter."/Double_neg_qualrecur_pass_".$sample_name;
	$file_double_neg_qualrecur_fail = $dir_inter."/Double_neg_qualrecur_fail_".$sample_name;
	open(FDP,">$file_double_call_pos")||die;
	open(FDN,">$file_double_call_neg")||die;
	open(FDPP,">$file_double_pos_qualrecur_pass")||die;
	open(FDPF,">$file_double_pos_qualrecur_fail")||die;
	open(FDNP,">$file_double_neg_qualrecur_pass")||die;
	open(FDNF,">$file_double_neg_qualrecur_fail")||die;
	open(FRF,"$file_annovar_rnaedit_fil")||die;
	while(<FRF>){
		my @arr = split /\t/,$_;
		my $var_qual = $arr[10]; my $site_info = $arr[12];
		my $var_info = $arr[0]."_".$arr[1]."_".$arr[3]."_".$arr[4];
		if(defined $snv_strelka{$var_inf}){
			print FDP $_;
			if($var_qual < 30){
				print FDPF $_; next;
			}
			if($site_info=~/;DP=(\d+)/){
				if($1 < $dp_cut){
					print FDPF $_; next;
				}
			}
			if($site_info=~/;FS=(\d\S*);MLEAC/){
				if($1 > 30){
					print FDPF $_; next;
				}
			}
			if($site_info =~/BaseQRankSum=(\S+\d);ClippingRankSum=(\S*\d);DP/){ ##Not all variants have this info
				if($1 < -2.33 || $1 > 2.33 || $2 < -2.33 || $2 > 2.33){ #p-value < 0.01
					print FDPF $_; next;
				}
			}
			if($info =~/MQRankSum=(\S+);QD/){
				if($1 < -2.33 || $1 > 2.33){ 
					print FDPF $_; next;
				}
			}
			if($info =~/ReadPosRankSum=(\S+);SOR/){
				if($1 < -2.33 || $1 > 2.33){
					print FDPF $_; next;
				}
			}
			#print FDPP $_;need recurrence to pass
			push (@double_pos_test_recur, $_);

		}else{
			print FDN $_;
			push (@double_neg_test_recur, $_);
		}
		next;
	}
	close(FRF);
	close(FDN);
	close(FDP);



