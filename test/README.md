# Test

The test dataset contains two samples from SRA. One sample should fail, and one should complete with a "bad" QC status.
The samples are from the v1 Swift panel, so the coverage isn't expected to be 99.7 %, and that's the reason for the
"bad" status. 

The pipeline should be able to run with just docker and nextflow installed, but to use docker, the 
`--use-docker` option is required. The script `test.sh` will run the pipeline test on two samples, 
and write the output and work folder to a subdirectory `run` under this directory (`test`).
The test uses a special `--test` option that overrides the data input in `main.nf`, to load data from 
SRA.

To run the test on a system with fewer than 96 cores / 200 GB RAM, you can use `lower-resources.config`. You
also may want to specify `--use_docker` to use Docker containers instead of singularity.

    bash test.sh -c lower-resources.config --use_docker

The default is to use the default resources, because the test may be used to validate the pipeline on a cluster.


It compares the output report to a [reference output](fasit/report_v8.tsv).

You can re-run the test multiple times, and it will cache the files and the completed jobs.
