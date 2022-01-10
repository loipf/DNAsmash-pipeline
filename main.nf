
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
	RUN_SMASH;
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

	// // a lot of trouble trying to filter path values and list input flows - complicated handling in modules
	//channel_samples = Channel.fromPath(params.sample_list)
	//		.splitText() { if(it?.trim()) { file(it.trim()) } }		
	//		.map { it -> tuple( it.getName(), it.listFiles() )  }
	//		.branch {
	//			fq: it[1].any{ it =~ /.*\.(fastq\.gz|fq\.gz)$/ }
	//			bam: it[1].any{ it =~ /.*\.(bam|bam\.bai)$/ }
	//			other: true
	//		}
	
	// TODO can be optimized to not copy everything into, but filter for following -> string-path struggles
	//.map { it -> tuple( it.getName(), it, it.list().findAll{ it =~ /.*\.(fastq\.gz|fq\.gz|bam|bam\.bai)$/ } )   }



	channel_samples = Channel.fromPath(params.sample_list)
			.splitText() { if(it?.trim()) { file(it.trim()) } }		
			.ifEmpty { error "cannot find any entries in matching: ${params.sample_list}" }
			.take( params.dev_samples )  // only consider a few files for debugging
			.branch {
				fq: it.list().any{ it =~ /.*\.(fastq\.gz|fq\.gz)$/ }
				bam: it.list().any{ it =~ /.*\.(bam|bam\.bai)$/ }
				other: true
			}


	URMAP_CREATE_INDEX(params.ensembl_release) 
	MAPPING_URMAP(channel_samples.fq, params.num_threads, URMAP_CREATE_INDEX.out.urmap_index)
	RUN_SMASH(MAPPING_URMAP.out.reads_mapped.collect(), channel_samples.bam.collect(), params.num_threads)
	
}



workflow.onComplete { 
	println ( workflow.success ? "\ndone!\n" : "oops .. something went wrong" )
} 















