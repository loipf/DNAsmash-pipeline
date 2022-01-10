# DNAsmash pipeline

a fast paired-end DNA mapping pipeline using URMAP from `.fastq` or `.bam` files for sample swap detection of matched files using SMaSH SNPs in maftools


---
### set up pipeline


before running, you have to set up the attached Docker image:
```sh
docker build -t dnasmash-pipeline https://raw.githubusercontent.com/loipf/DNAsmash-pipeline/master/docker/Dockerfile
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker dnasmash-pipeline` argument.


---
### run pipeline

URMAP index must fit into memory, so at least 32Gb RAM are necessary.

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAsmash-pipeline -r main --project_dir /path/to/folder --sample_list /path/to/list --num_threads 10 -with-docker dnasmash-pipeline
```
for this execution to work properly, you have to be in the current project directory.

the `--sample_list` option file must contain the folder path of each sample where the raw reads (`.fq.gz`|`.fastq.gz`) or already mapped reads (`.bam`+`.bam.bai`) are located (without any empty lines):
```sh
/path/to/reads_raw/sample1
/path/to/reads_mapped/sample3
/path/to/reads_raw/sample2
/path/to/reads_mapped/sample4
```


optional extendable with:
```sh
-resume
-with-report report_DNAsmash-pipeline
-with-timeline timeline_DNAsmash-pipeline
-w work_dir
```

by default, all output will be saved into the `data` folder of the current directory.
best to run with a new clear folder structure as not all new results do overwrite old ones.





