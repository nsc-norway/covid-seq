#!/bin/sh

# Script usage:
# bash test.sh [{ADDITIONAL_OPTIONS_TO_NEXTFLOW}]


if [ ! -d test_data ]
then
    echo "Downloading data..."
    if ! nextflow run download-testdata.nf
    then
        echo "ERROR: Downloading test data failed. Make sure to use NCBI API key."
        exit 1
    fi
else
    echo "Test data is already downloaded."
fi

mkdir run || echo "Directory 'run/' already exists, will resume execution."

if nextflow run ../main.nf --lab test --samplelist sampleList_test.csv --outpath run -resume -w run/work "$@"
then
    echo " -- Comparing report file with refrence --"
    if [ -f run/results/report_v*.tsv ] && diff run/results/report_v*.tsv fasit/report.tsv
    then
        echo " -----"
        echo "TEST PASSED"
    else
        echo "TEST FAILED -- Output doesn't match reference output"
        exit 1
    fi
else
    echo "TEST FAILED -- Test aborted because of an error in the pipeline --"
    exit 1
fi
