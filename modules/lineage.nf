
process PANGOLIN { 
    tag "$sampleName"

    label 'small'

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

process NEXTCLADE { 
    tag "$sampleName"

    label 'small'

    input:    
    tuple val(sampleName), path(consensus)
    val(caller)

    output:
    tuple val(sampleName), val(caller), path("${sampleName}_${caller}_nextclade.csv"), emit: NEXTCLADE_out
    path "*.{log,sh}"

    publishDir "${params.outdir}/5_lineage/nextclade", mode: 'link', pattern:'*.csv'
    publishDir "${params.outdir}/5_lineage/nextclade/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    nextclade.js --version
    
    nextclade.js -j $task.cpus \
        -i $consensus \
        -c ${sampleName}_${caller}_nextclade.csv 
    
    cp .command.sh ${sampleName}.${caller}.nextclade.sh
    cp .command.log ${sampleName}.${caller}.nextclade.log
    """
}