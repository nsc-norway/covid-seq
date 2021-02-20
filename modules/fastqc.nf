process FASTQC {
    tag "$sampleName"
    
    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    val(source)
    
    publishDir "${params.outdir}/1_fastq", mode: 'link', pattern:'*fq.gz'
    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'*.{html,zip,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}*zip"), emit: FASTQC_out
    path "*.html"
    path "*.{log,sh}"
    
    script:
    """
    fastqc -t $task.cpus ${read1} ${read2}
    cp .command.log ${sampleName}.${source}.fastqc.log
    cp .command.sh ${sampleName}.${source}.fastqc.sh
    """
}
