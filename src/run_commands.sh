sbatch --mem=2G --partition=park --account=park_contrib -n 1 --time=1-00:00 --wrap "perl filter_and_label.pl ../input_data/Horning_LNCaP_GSE99795/ ../results/Horning_LNCaP_GSE99795/ ~/park_lab/Ning/GEO/Horning_Prostate_LNCaP_cell_line_scRNAseq_CancerResearch_2018_GSE99795/results/RNA-seq/list_samples_0houruntreatedcell.txt LNCaP 3 20" 
sbatch --mem=2G --partition=park --account=park_contrib -n 1 --time=6-00:00 --wrap "perl filter_and_label.pl ../input_data/Rambow_melanoma/ ../results/Rambow_melanoma/ ../input_data/Rambow_melanoma/list_samples.txt Rambow_melanoma 3 3 539"
sbatch --mem=2G --partition=park --account=park_contrib -n 1 --time=1-00:00 --wrap "perl filter_and_label.pl ../input_data/JURKAT_SMART/ ../results/JURKAT_SMART/ ~/park_lab/Ning/GEO/Rodriguez_SET2_cell_line_scRNA-seq_2019_GSE105451/results/RNA-seq/list_samples_JURKAT_SMART.txt JURKAT_SMART 3 3 15"



