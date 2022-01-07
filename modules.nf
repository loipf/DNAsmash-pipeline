
// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"


process URMAP_CREATE_INDEX { 

	input:
		val ensembl_release

	output:
		path "Homo_sapiens.GRCh38.dna.genome_smash.ufi", emit: urmap_index

	shell:
	'''
	
	### download necessary genome (without MT and Y)
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
	#wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
	
	cat Homo_sapiens.GRCh38.dna.chromosome.* > Homo_sapiens.GRCh38.dna.genome_smash.fa.gz
	rm Homo_sapiens.GRCh38.dna.chromosome.*
	gunzip Homo_sapiens.GRCh38.dna.genome_smash.fa.gz
	
	### create index
	/usr/src/urmap -make_ufi Homo_sapiens.GRCh38.dna.genome_smash.fa -veryfast -output Homo_sapiens.GRCh38.dna.genome_smash.ufi

	'''
}



process MAPPING_URMAP { 

	input:
		tuple val(sample_id), val(absolute_path), val(file_names) 
		val num_threads
		//path urmap_index

	output:
		//path "*.bam", emit: reads_mapped
		

	shell:
	'''
	
	echo !{sample_id}
	echo !{absolute_path}
	echo !{file_names}
	
	prefix_path=!{absolute_path}/
	prefix_path="test_"
	echo $prefix_path
	
	#read_files=( "${!{file_names}[@]/#/${prefix_path}}" )
	#echo $read_files
	
	test=( !{file_names})
	echo $test
	
	#reads_array=$(echo !{file_names} | sed "s/, / /g; s/[//g"   )
	reads_array=($(echo !{file_names} | sed 's/[][]//g; s/,//g'   ))
	echo $reads_array
	echo ${reads_array[@]}
	
	ARR_PRE=("${reads_array[@]/#/!{absolute_path}/}")
	echo ${ARR_PRE[@]}
	
	### theoretically sorting not needed but maybe speeds up things
	#reads_sorted=$(echo !{file_names} | xargs -n1 | sort | xargs) 
	#echo $reads_sorted
	
	#read_files=( "${!{file_names}[@]/#/${prefix_path}}" )
	#echo $read_files
	
	### combine multiple seq files in the same sample directory with same direction together
	#reads_sorted_1=$(find $reads_sorted -name "*_1.fq.gz" -o -name "*_1.fastq.gz")
	#reads_sorted_2=$(find $reads_sorted -name "*_2.fq.gz" -o -name "*_2.fastq.gz")
	#cat $reads_sorted_1 > raw_reads_connected_1.fastq.gz
	#cat $reads_sorted_2 > raw_reads_connected_2.fastq.gz
	
	samtools samtools
	
	# /usr/src/urmap -veryfast -threads 10 -ufi Homo_sapiens.GRCh38.dna.genome_smash.ufi -map2 /home/stefanloipfinger/Documents/impact/metastasis_rnaseq/reads_raw/B8/B8_1.fq.gz -reverse /home/stefanloipfinger/Documents/impact/metastasis_rnaseq/reads_raw/B8/B8_2.fq.gz -samout sample2_mapped.sam

	## make BAM - try with pipe later
	#samtools view -h -b -@ 10 sample2_mapped.sam | samtools sort -@ 10 > sample2_mapped.bam
	#samtools index -b -@ !{num_threads} sample_mapped.bam
	
	'''
}









process CREATE_T2G_LIST { 
	input:
		path raw_transcripts

	output:
		path "transcript_to_gene_list.csv", emit: t2g_list

	"""
	gunzip -c $raw_transcripts > raw_transcripts.fa

	perl -lne 'print \$1 while /^>(.+gene_symbol:\\S+)/gi' raw_transcripts.fa | sed 's/\\s/,/g' > transcript_to_gene_list.csv
	sed  -i '1i transcript_id,transcript_biotype,chromosome,gene_id,gene_biotype,transcript_biotype,gene_symbol' transcript_to_gene_list.csv  ### add header
	sed -i 's/chromosome://g; s/gene://g; s/gene_biotype://g; s/transcript_biotype://g; s/gene_symbol://g;' transcript_to_gene_list.csv ### remove pre-strings

	"""
}


// removes duplicate transcripts before merging into gene names
process RM_DUPLICATE_TRANSCRIPTS { 
	publishDir "$params.data_dir/kallisto_index", pattern:"Homo*", mode: "copy"

	input:
		path raw_transcripts 

	output:
		path "Homo_sapiens.GRCh38.cdna_ncrna_oneline_unique.txt", emit: trans_oneline_unique
		path "kallisto_removal_info.txt", emit: removal_info

	shell:
	'''
	gunzip -c !{raw_transcripts} | cut -f1 -d" " | awk '/^>/ {printf("\\n%s\\n",$0);next; } { printf("%s",$0);} END {printf("\\n");}' | tail -n +2 | sed 'N;s/\\n/ /' > Homo_sapiens.GRCh38.cdna_ncrna_oneline.txt
	
	remove_duplicate_transcripts.R Homo_sapiens.GRCh38.cdna_ncrna_oneline.txt
	'''
}



process PREPROCESS_READS { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", pattern:"*cutadapt_output.txt", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_3_seq_file
		path adapter_5_seq_file

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''
	
	### combine multiple seq files in the same sample directory with same direction together
	### theoretically sorting not needed but to be ordered
	reads_sorted=$(echo !{reads} | xargs -n1 | sort | xargs)
	reads_sorted_1=$(find $reads_sorted -name "*_1.fq.gz" -o -name "*_1.fastq.gz")
	reads_sorted_2=$(find $reads_sorted -name "*_2.fq.gz" -o -name "*_2.fastq.gz")
	cat $reads_sorted_1 > raw_reads_connected_1.fastq.gz
	cat $reads_sorted_2 > raw_reads_connected_2.fastq.gz

	if [[ (!{adapter_3_seq_file} == "NO_FILE") || (!{adapter_5_seq_file} == "NO_FILE2") ]]
	then
		### run with adapter sequences
		
		### adapter 3' input as string or file
		if [[ !{adapter_3_seq_file} == *".fasta"* ]]
		then
	  		ADAPTER_3=file:!{adapter_3_seq_file}
		else
			ADAPTER_3="!{adapter_3_seq_file}"
		fi
		### adapter 5' input as string or file
		if [[ !{adapter_5_seq_file} == *".fasta"* ]]
		then
	  		ADAPTER_5=file:!{adapter_5_seq_file}
		else
			ADAPTER_5="!{adapter_5_seq_file}"
		fi
		
		cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -a $ADAPTER_3 -A $ADAPTER_5 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz raw_reads_connected_1.fastq.gz raw_reads_connected_2.fastq.gz > !{sample_id}_cutadapt_output.txt
		
	else
		### run without adapter sequence
		cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz raw_reads_connected_1.fastq.gz raw_reads_connected_2.fastq.gz > !{sample_id}_cutadapt_output.txt
	fi
	
	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	fastqc -t !{num_threads} --noextract !{reads}
	'''
}



process FASTQC_READS_PREPRO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	fastqc -t !{num_threads} --noextract !{reads_prepro}
	'''
}


process QUANT_KALLISTO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_quant", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	
	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path kallisto_index

	output:
		path "${sample_id}_abundance.h5", emit: kallisto_abundance
		path "${sample_id}_run_info.json", emit: kallisto_json
		path "${sample_id}_kallisto_output.txt", emit: kallisto_output

	shell:
	'''

	kallisto quant -i !{kallisto_index} -o . -t !{num_threads} !{reads_prepro} &> kallisto_output.txt

	for f in *; do mv -- "$f" "!{sample_id}_${f%}"; done   # rename all with sample_id prefix

	'''
}



// removes duplicate transcripts before merging into gene names
process CREATE_KALLISTO_QC_TABLE { 
	input:
		path kallisto_json

	output:
		path "kallisto_aligned_reads_qc.csv", emit: kallisto_qc_table

	shell:
	'''
	create_kallisto_qc_table.R
	'''
}


process CREATE_GENE_MATRIX { 
	publishDir "$params.data_dir", mode: "copy", saveAs: { filename -> filename.endsWith(".rds") ? "reads_quant/$filename" : filename }

	input:
		path kallisto_qc_table
		path removal_info
		path trans_oneline_unique
		path t2g_list
		path kallisto_abundance
		
	output:
		path "kallisto_gene_counts.csv", emit: gene_matrix
		path "kallisto_gene_counts_norm_sf_vst.csv", emit: gene_matrix_vst
		path "kallisto_aligned_reads_qc.csv", emit: kallisto_qc_table
		path "all_kallisto_abundance_obj.rds"
		path "kallisto_removal_info.txt"
		path "transcript_to_gene_list.csv"

	shell:
	'''
	create_kallisto_gene_matrix.R !{kallisto_qc_table} !{removal_info} !{trans_oneline_unique} !{t2g_list} 
	'''
}




process CREATE_GENE_COUNT_PLOTS { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path kallisto_qc_table
		path gene_matrix
		path gene_matrix_vst
		
	output:
		path "gene_count_analysis_plots.html"
	
	shell:
	'''
	#!/usr/bin/env Rscript

	curr_dir = getwd()
	rmarkdown_file_path = system("which create_gene_count_plots.R", intern = TRUE)
	
	params_list = list(curr_dir = curr_dir, reads_qc_table_file='!{kallisto_qc_table}', gene_matrix_file='!{gene_matrix}',gene_matrix_file_vst='!{gene_matrix_vst}')

	rmarkdown::render(rmarkdown_file_path, output_file=file.path(curr_dir,'gene_count_analysis_plots.html'), params= params_list )

	'''
}





process MULTIQC_RAW { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_raw .
	'''
}


process MULTIQC_PREPRO { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_prepro .
	'''
}


process MULTIQC_QUANT { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_quant .
	'''
}


















