# Test

The test dataset contains two samples from SRA. One sample should fail, and one should complete with a "bad" QC status.
The samples are from the v1 Swift panel, so the coverage isn't expected to be 99.7 %, and that's the reason for the
"bad" status. 

The pipeline should be able to run with just docker and nextflow installed, but to use docker, the 
`--use_docker` option is required. The script `test.sh` will run the pipeline test on two samples, 
and write the output and work folder to a subdirectory `run` under this directory (`test`).
The test script uses a dedicated nextflow workflow to download the data into `test_data/`. It
then runs the normal main.nf with a sample list that refers to those downloaded files. [Previously,
main.nf had a special "--test" option, but this has been removed due to caching issues.]

To run the test on a system with fewer than 64 cores / 200 GB RAM, you can use `lower-resources.config`. You
also may want to specify `--use_docker` to use Docker containers instead of singularity.

    bash test.sh -c lower-resources.config --use_docker

The default is to use the standard resource allocations, because the test may be used to validate the pipeline on a cluster.


The test script compares the output report to a [reference output](fasit/report.tsv).

You can re-run the test multiple times, and Nextflow will cache the files and the completed jobs.


## Updating test reference after changes in output

If there are expected changes in the report content, update the file in `fasit/` by replacing it with the
new report file in `run/results/report_vX.tsv` (after running the test, which would have failed).
