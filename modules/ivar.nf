process IVAR_VARIANTS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(mpileup)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_ivar.vcf.gz"), path ("${sampleName}_ivar.vcf.gz.tbi"), emit: IVAR_VARIANTS_out
    path "*.txt", emit: IVAR_BCFTOOLS_STATS_out
    path "*.{vcf,ivar_variants.log,sh}"

    publishDir "${params.outdir}/3_variants/ivar", mode: 'link', pattern:'*.{vcf}'
    publishDir "${params.outdir}/3_variants/ivar/log", mode: 'link', pattern:'*.{txt,ivar_variants.log,sh}'

    script:
    """
    cat $mpileup | ivar variants -q 20 -m 10 -r $genome -p ${sampleName}_ivar -

    ivar_variants_to_vcf.py ${sampleName}_ivar.tsv ${sampleName}_ivar.vcf > ${sampleName}_ivar.variant.counts.log
    bgzip -c ${sampleName}_ivar.vcf > ${sampleName}_ivar.vcf.gz
    tabix -p vcf -f ${sampleName}_ivar.vcf.gz

    bcftools stats ${sampleName}_ivar.vcf.gz > ${sampleName}_ivar.bcftools_stats.txt

    cp .command.sh ${sampleName}.ivar_variants.sh
    cp .command.log ${sampleName}.ivar_variants.log
    """
}


process IVAR_CONSENSUS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(mpileup)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_ivar.consensus.masked.fa"), path ("${sampleName}_ivar.consensus.masked.qual.txt"), emit: IVAR_CONSENSUS_out
    tuple val(sampleName), path ("${sampleName}_ivar.consensus.masked.fa"), emit: FOR_LINEAGE_out
    path "${sampleName}_ivar.consensus.masked.fa", emit: IVAR_FASTA_ONLY_out
    path "*.{log,sh}"

    publishDir "${params.outdir}/4_consensus/ivar", mode: 'link', pattern:'*.{consensus.masked.fa,consensus.masked.qual.txt}'
    publishDir "${params.outdir}/4_consensus/ivar/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    cat $mpileup | ivar consensus -q 20 -m 10 -n N -p ${sampleName}_ivar.consensus.masked

    header=\$(head -n1 ${sampleName}_ivar.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_ivar_masked/g" ${sampleName}_ivar.consensus.masked.fa

    cp .command.sh ${sampleName}.ivar_consensus.sh
    cp .command.log ${sampleName}.ivar_consensus.log
    """
}
