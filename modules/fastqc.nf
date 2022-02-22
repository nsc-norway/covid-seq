process FASTQC {
    tag "$sampleName"
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    val(source)
    
    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'*.{log,sh}'

    output:
    tuple val(sampleName), path ("*zip"), emit: FASTQC_out
    path "*.{log,sh}"
    
    script:
    """
    fastqc -t $task.cpus ${read1} ${read2}
    cp .command.log ${sampleName}.${source}.fastqc.log
    cp .command.sh ${sampleName}.${source}.fastqc.sh
    """
}
