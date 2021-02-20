process PICARD_WGSMETRICS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(bam), path(bai)
    path genome
    
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{CollectWgsMetrics.txt,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.trim.CollectWgsMetrics.txt"), emit: PICARD_WGSMETRICS_out

    path "*.{log,sh}"

    script:
    """
    picard -Xmx${task.memory.giga}g -XX:ParallelGCThreads=$task.cpus \
        CollectWgsMetrics \
        COVERAGE_CAP=1000000 \
        INPUT=$bam \
        OUTPUT=${sampleName}.trim.CollectWgsMetrics.txt \
        REFERENCE_SEQUENCE=$genome \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tmp

    cp .command.sh ${sampleName}.CollectWGSmetrics.sh
    cp .command.log ${sampleName}.CollectWGSmetrics.log
    """
}
