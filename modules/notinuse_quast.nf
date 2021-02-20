process QUAST {
    tag "$sampleName"

    input:
    tuple val(sampleName), path(consensus)
    val(caller)
    path genome

    output:
//    tuple val(sampleName), path ("${sampleName}.trim.mkD.mpileup")
//    path "*.{log,sh}"

//    publishDir "${params.outdir}/4_consensus/quast", mode: 'link', pattern:'*.csv'
//    publishDir "${params.outdir}/4_consensus/quast/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    quast.py \
        --output-dir ${sampleName}.${caller}.quast \
        -r $genome \
        --threads $task.cpus \
        $consensus

    cp .command.sh ${sampleName}.${caller}.quast.sh
    cp .command.log ${sampleName}.${caller}.quast.log
    """
}

