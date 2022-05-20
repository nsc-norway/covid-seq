
process PANGOLIN { 
    tag "$sampleName"

    label 'tiny'

    input:
    tuple val(sampleName), path(consensus)
    val(caller)

    output:
    tuple val(sampleName), val(caller), path ("${sampleName}_${caller}_pangolin.csv"), emit: PANGOLIN_out
    path "*.{log,sh}"

    publishDir "${params.outdir}/5_lineage/pangolin", mode: 'link', pattern:'*.csv'
    publishDir "${params.outdir}/5_lineage/pangolin/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    pangolin $consensus -t $task.cpus
    mv lineage_report.csv ${sampleName}_${caller}_pangolin.csv 

    cp .command.sh ${sampleName}.${caller}.pangolin.sh
    cp .command.log ${sampleName}.${caller}.pangolin.log
    """
}
