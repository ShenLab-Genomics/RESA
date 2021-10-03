#!/usr/bin/perl -w
#From CCLE reported mutation, test max detectable mutation after qc filtering

my $cell =$ARGV[0];
my $condition = $ARGV[1];
my $dir_WES = $ARGV[2]; #dir_data."WES_CCLE/";
my $dir_RNA = $ARGV[3]; #"../temp/folder_int_files/";
my $dir_data = $ARGV[4]; ## input data #"/home/ns362/park_lab/Ning/GEO/Horning_Prostate_LNCaP_cell_line_scRNAseq_CancerResearch_2018_GSE99795/results/";
my $dir_out = $ARGV[5]; 
unless(-e $dir_out){`mkdir $dir_out`;}
my $list_samples = $ARGV[6]; # $dir_data."RNA-seq/list_samples.txt";

print "Reading WES data...\n";
	my $file_wes = $dir_WES.$cell."_mutations_hg38lifted.bed";
	my %saw_var_wes = ();
	open(W,"$file_wes")||die;
	while(<W>){
		my (@arr)=split /\s+/,$_;
		unless($arr[3]=~/SNP/){next;}
		$var = $arr[0]."\t".$arr[1]."\t".$arr[3];
		$saw_var_wes{$var}=1;
	}
	close(W);

#	foreach my $exp(@arr_exp){
#                my $condition = $cell."_".$exp;
print "Reading scRNA-seq...\n";
                my @arr_samples = split /\n/, `grep $condition $list_samples`;
                foreach my $sample(@arr_samples){
			my $f_rna = $dir_RNA."folder_int_files/Filter_annovar_rnaedit_fil_".$sample.".hg38_generic_filtered";
			my $count=0;
			open(FR,"$f_rna")||next;
			while(<FR>){
#				chomp;
				my (@tmps)=split /\t/,$_;
				my $test_var = $tmps[0]."\t".$tmps[1]."\t"."SNP_".$tmps[3]."_".$tmps[4];
				if(defined $saw_var_wes{$test_var}){
					$count++;
					$final_var{$test_var}=1;
#					print FF $_;
				}
				next;
			}
			close(FR);
#print $sample,"\t",$count,"\n";die;
			my $f_rna_vcf = $dir_data."variant_call_ctat/$sample".".vcf.gz";
			my $f_rna_fil = $dir_out.$sample."_SNV_scRNA_WES_confirmed.vcf";
			open(FF,">$f_rna_fil")||die;
			open(FV,"zcat $f_rna_vcf |")||die;
			while(<FV>){
				if($_=~/^#/){
					print FF $_;
				}else{
#					chomp;
					my (@tmps)=split /\t/,$_;
					my $test_var = $tmps[0]."\t".$tmps[1]."\t"."SNP_".$tmps[3]."_".$tmps[4];
					if(defined $final_var{$test_var}){
						print FF $_;
					}
				}
				next;
			}
			close(FV);
			close(FF);
			print "For sample $sample, total SNVs confirmed by WES is: ",$count,"\n";
		}


