process NSCTRIM {
    tag "$sampleName"
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path(primer_file)

    output:
    tuple val(sampleName), path ("${sampleName}_nscTrim_R1.fq.gz"), path ("${sampleName}_nscTrim_R2.fq.gz"), emit: NSCTRIM_out
    path "*.log", emit: NSCTRIM_log
    path "*.{txt,sh}"

    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'*.{txt,log,sh}'

    script:
    """
    NSCtrim \
        --mismatches-per-primer 1 \
        --swapped-primer-pairs \
        <( awk '{ print \$1, \$2, \$3 }' ${primer_file} ) \
        ${read1} ${read2} \
        ${sampleName}_nscTrim_R1.fq.gz ${sampleName}_nscTrim_R2.fq.gz \
        > ${sampleName}.nscTrimStats.txt

    cp .command.sh ${sampleName}.nsc_trim.sh
    cp .command.log ${sampleName}.nsc_trim.log
    """
}