# DNAsmash pipeline

a fast DNA mapping pipeline from `.fastq` files for sample swap detection of matched files using URMAP and SMaSH


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
nextflow run loipf/DNAsmash-pipeline -r main --project_dir /path/to/folder --reads_dir /path/to/samples --num_threads 10 -with-docker dnasmash-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_DNAsmash-pipeline
-with-timeline timeline_DNAsmash-pipeline
-w work_dir
```

by default, all output will be saved into the `data` folder of the current directory.
best to run with a new clear folder structure as not all new results do overwrite old ones.





