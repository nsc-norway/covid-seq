#!/bin/sh

mkdir run || echo "Warning: 'run/' dir exists. Will resume execution."

if nextflow run ../main.nf --test --lab test --samplelist sampleList_test.csv --use_docker --outpath run -resume -w run/work
then
    echo " -- Comparing report file with refrence --"
    if diff run/results/report_v8.tsv fasit/report_v8.tsv
    then
        echo " -----"
        echo "TEST PASSED"
    else
        echo "TEST FAILED -- Output doesn't match reference output"
        exit 1
    fi
else
    echo "TEST FAILED -- Test aborted because of Nextflow error --"
    exit 1
fi