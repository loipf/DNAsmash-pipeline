
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
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
	wget http://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
	
	cat Homo_sapiens.GRCh38.dna.chromosome.* > Homo_sapiens.GRCh38.dna.genome_smash.fa.gz
	rm Homo_sapiens.GRCh38.dna.chromosome.*
	gunzip Homo_sapiens.GRCh38.dna.genome_smash.fa.gz
	
	### create index
	/usr/src/urmap -make_ufi Homo_sapiens.GRCh38.dna.genome_smash.fa -veryfast -output Homo_sapiens.GRCh38.dna.genome_smash.ufi

	'''
}



process MAPPING_URMAP { 
	tag "$sample_folder"

	input:
		path sample_folder
		val num_threads
		path urmap_index

	output:
		tuple path("*.bam"), path("*.bam.bai"), emit: reads_mapped

	shell:
	'''
	
	### theoretically sorting not needed but maybe speeds up things
	reads_sorted=$(ls -d !{sample_folder}/* | xargs -n1 | sort | xargs)
	
	### combine multiple seq files in the same sample directory with same direction together
	reads_sorted_1=$(find $reads_sorted -name "*_1.fq.gz" -o -name "*_1.fastq.gz")
	reads_sorted_2=$(find $reads_sorted -name "*_2.fq.gz" -o -name "*_2.fastq.gz")
	
	# TODO: make cat copy operation optional, reads < urmap index size lead to error
	cat $reads_sorted_1 | seqkit seq --min-len 10 - -j !{num_threads} -o raw_reads_connected_1.fastq.gz
	cat $reads_sorted_2 | seqkit seq --min-len 10 - -j !{num_threads} -o raw_reads_connected_2.fastq.gz

	#/usr/src/urmap -veryfast -threads !{num_threads} -ufi !{urmap_index} -map2 raw_reads_connected_1.fastq.gz -reverse raw_reads_connected_2.fastq.gz -samout - \
	#| samtools view -b -@ !{num_threads} - \
	#| samtools sort -@ !{num_threads} > !{sample_folder}.bam
	
	/usr/src/urmap -veryfast -threads !{num_threads} -ufi !{urmap_index} -map2 raw_reads_connected_1.fastq.gz -reverse raw_reads_connected_2.fastq.gz -samout !{sample_folder}.sam
	
	samtools view -b -@ !{num_threads} !{sample_folder}.sam | samtools sort -@ !{num_threads} > !{sample_folder}.bam
	samtools index -b -@ !{num_threads} !{sample_folder}.bam
	

	'''
}



process RUN_SMASH { 
	publishDir "$params.data_dir", mode: "copy"

	input:
		path urmap_bam_files
		path bam_files
		val num_threads
		val min_snp_read_depth
		
	output:
		path "sample_swap_results.rds"
		path "matching_samples.txt"
		path "matching_samples_dendrogram_clustering.csv"
		path "sample_swap_corr_heatmap_all.pdf"
		path "sample_swap_corr_heatmap_single_removed.pdf"

	shell:
	'''
	run_sample_swap_identification.R !{num_threads} !{min_snp_read_depth}
	'''
}















