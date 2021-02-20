process VARSCAN2_VARIANTS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(mpileup)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_varscan2.vcf.gz"), path ("${sampleName}_varscan2.vcf.gz.tbi"), emit: VARSCAN2_VARIANTS_out
    path "*.{vcf,txt,log,sh}"

    publishDir "${params.outdir}/3_variants/varscan2", mode: 'link', pattern:'*.{vcf}'
    publishDir "${params.outdir}/3_variants/varscan2/log", mode: 'link', pattern:'*.{txt,log,sh}'

    script:
    """
    echo "${sampleName}_varscan2" > sample_name.list
    varscan mpileup2cns \
        ${sampleName}.trim.mpileup \
        --min-coverage 10 \
        --min-reads2 5 \
        --min-avg-qual 20 \
        --min-var-freq 0.25 \
        --p-value 0.99 \
        --output-vcf 1 \
        --vcf-sample-list sample_name.list \
        --variants \
        --strand-filter 0 \
        > ${sampleName}_varscan2.vcf \
        2> ${sampleName}.varscan2_variants.log

    bgzip -c ${sampleName}_varscan2.vcf > ${sampleName}_varscan2.vcf.gz
    tabix -p vcf -f ${sampleName}_varscan2.vcf.gz

    bcftools stats ${sampleName}_varscan2.vcf.gz > ${sampleName}_varscan2.bcftools_stats.txt

    cp .command.sh ${sampleName}.varscan2_variants.sh
    """
}


process VARSCAN2_CONSENSUS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path (bam), path (bai), path (vcf), path (tbi)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_varscan2.consensus.masked.fa"), emit: VARSCAN2_CONSENSUS_out
    tuple val(sampleName), path ("${sampleName}_varscan2.consensus.masked.fa"), emit: FOR_LINEAGE_out
    path "*.{log,sh}"

    publishDir "${params.outdir}/4_consensus/varscan2", mode: 'link', pattern:'*.masked.fa'
    publishDir "${params.outdir}/4_consensus/varscan2/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    bedtools genomecov \
        -bga \
        -ibam $bam \
        -g $genome \
        | awk '\$4 < 10' | bedtools merge > ${sampleName}.pre_mask.bed
    
    parse_mask_bed.py $vcf ${sampleName}.pre_mask.bed ${sampleName}.mask.bed
    
    bedtools maskfasta \
        -fi $genome \
        -bed ${sampleName}.mask.bed \
        -fo ${sampleName}_varscan2.ref.masked.fa

    cat ${sampleName}_varscan2.ref.masked.fa | bcftools consensus $vcf > ${sampleName}_varscan2.consensus.masked.fa

    header=\$(head -n 1 ${sampleName}_varscan2.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_varscan2_masked/g" ${sampleName}_varscan2.consensus.masked.fa

    cp .command.sh ${sampleName}.varscan2_consensus.sh
    cp .command.log ${sampleName}.varscan2_consensus.log
    """
}
