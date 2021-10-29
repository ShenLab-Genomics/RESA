

1.Create a container based on the image
	# pull the image
	docker pull xuguorong2016/resa:latest

	# create the container of the image 
	docker run -it xuguorong2016/resa:latest /bin/bash

	## After entering the container, root@containerID will show on the left


2.Install the Annovar and download the customize library 
	# register on https://www.openbioinformatics.org/annovar/annovar_download_form.php and 
	# download the install package. 

	# copy annovar.latest.var 
	docker cp annovar.latest.tar.gz containerID:/bin/annovar.latest.tar.gz

	#enter the container
	docker exec -it containerID /bin/bash

	#Install the annovar
	tar -zxvf /bin/annovar.latest.tar.gz

	# download datasets
	cd /bin/annovar/

	## hg38_gnomad30_genome.txt is a large file, 
	## and make sure the available space is greater than 80G
	./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/

	./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 ensGene/
	./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
	
	## hg38_RNAedit.txt cannot be downloaded directly from the annovar webpage, 
	## and here shows one way to convert it from GRCh38.RNAediting.vcf.gz 
	apt-get install wget
	wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/GRCh38.RNAediting.vcf.gz
	sed -i '1,5d' GRCh38.RNAediting.vcf
	awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}' GRCh38.RNAediting.vcf > temp.txt 
	sed '1i\#Chr\tStart\tEnd\tRef\tAlt' temp.txt > hg38_RNAedit.txt

	
3.Install RESA 
	# download RESA from GitHub and copy to the dock container
	docker cp RESA containerID:/workshop/RESA/


4.Analysis the example dataset K562_primers with RESA
	# RESA-filter
	mkdir /workspace/RESA/examples/result/
	cd /workspace/RESA/src
	perl filter_and_label_v2_1.pl \
	/workspace/RESA/examples/input/K562_primers/ \
	/workspace/RESA/examples/result/ \
	/workspace/RESA/examples/input/K562_primers/sample_list.txt K562_primers 3 3 42

	# RESA-refine
	python3.8 RESA_refine.py \
	-N /workspace/RESA/examples/result/K562_primers/folder_confident/Double_neg_qualrecur_fail.txt \
	-P /workspace/RESA/examples/result/K562_primers/folder_confident/Double_pos_qualrecur_pass.txt \
	-U /workspace/RESA/examples/result/K562_primers/folder_unsure/Undefined_B.txt \
	-O /workspace/RESA/examples/result/K562_primers/
	
	## Predicted mutations have saved in prdicted_value.csv

#Tips
# restart container
docker container start containerID

# If the R packages MutationalPatterns cannot install automatically,
# we suggest install it on R before run the RESA_refine.py 
# and make sure memory is greater than 5GB.




