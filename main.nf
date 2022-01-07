
/* 
 * DNAsmash PIPELINE 
 * for paired-end reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	URMAP_CREATE_INDEX;
	MAPPING_URMAP;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = -1

params.project_dir	= "$projectDir"
params.sample_list	= "$params.project_dir/sample_list.txt"
params.data_dir	= "$params.project_dir/data"



/*
 * other parameters
 */

params.num_threads	= 5
params.ensembl_release	= 101




log.info """\
DNAsmash PIPELINE
===================================================
sample_list		: $params.sample_list
data_dir		: $params.data_dir
ensembl_version	: $params.ensembl_release
===================================================

"""



/* 
 * main pipeline logic
 */
workflow {

	channel_samples = Channel.fromPath(params.sample_list)
			.splitText() { if(it?.trim()) { file(it.trim())} }		
			.map { it -> tuple( it.getName(), it.listFiles() )  }
			.ifEmpty { error "cannot find any entries in matching: ${params.sample_list}" }
			.take( params.dev_samples )  // only consider a few files for debugging
			.branch {
				fq: it[1].any{ it =~ /.*\.(fastq\.gz|fq\.gz)$/ }
				bam: it[1].any{ it =~ /.*\.(bam|bam\.bai)$/ }
				other: true
			}
			
	// TODO can be optimized to not copy everything into, but filter for following -> string-path struggles
	//.map { it -> tuple( it.getName(), it, it.list().findAll{ it =~ /.*\.(fastq\.gz|fq\.gz|bam|bam\.bai)$/ } )   }

		//channel_samples.view()
		//println "fq"
		channel_samples.fq.view()
		//println "bam"
		channel_samples.bam.view()
		//println "other"
		channel_samples.other.view()
	
	//channel_samples_other = channel_samples.other.toList()
			
	//channel_samples.fq.view()	


	//URMAP_CREATE_INDEX(params.ensembl_release) 
	//MAPPING_URMAP(channel_samples.fq, params.num_threads, URMAP_CREATE_INDEX.out.urmap_index)
	MAPPING_URMAP(channel_samples.fq, params.num_threads)	
	
	
	//PREPROCESS_READS(channel_reads, params.num_threads, params.adapter_3_seq_file, params.adapter_5_seq_file)
	//channel_reads_prepro = PREPROCESS_READS.out.reads_prepro.map{ it -> tuple(it[0], tuple(it[1], it[2])) }
	
	//QUANT_KALLISTO(channel_reads_prepro, params.num_threads, CREATE_KALLISTO_INDEX.out.kallisto_index)
	//CREATE_KALLISTO_QC_TABLE(QUANT_KALLISTO.out.kallisto_json.collect())
	//CREATE_GENE_MATRIX(CREATE_KALLISTO_QC_TABLE.out.kallisto_qc_table, RM_DUPLICATE_TRANSCRIPTS.out.removal_info, RM_DUPLICATE_TRANSCRIPTS.out.trans_oneline_unique, CREATE_T2G_LIST.out.t2g_list, QUANT_KALLISTO.out.kallisto_abundance.collect() )

	//CREATE_GENE_COUNT_PLOTS(CREATE_GENE_MATRIX.out.kallisto_qc_table, CREATE_GENE_MATRIX.out.gene_matrix, CREATE_GENE_MATRIX.out.gene_matrix_vst)

	//MULTIQC_RAW(FASTQC_READS_RAW.out.reports.collect() )

}



workflow.onComplete { 
	println ( workflow.success ? "\ndone!\n" : "oops .. something went wrong" )
	//println ( workflow.success ? "\ndone! following files not found\n" : "oops .. something went wrong" )
	//println("$channel_samples_other")
 } 















