process HAPLOTYPECALLER_VARIANTS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(bam)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_hapcaller.vcf.gz"), path ("${sampleName}_hapcaller.vcf.gz.tbi"), emit: hapcaller_VARIANTS_out
    path "*.{vcf,txt,log,sh}"

    publishDir "${params.outdir}/3_variants/haplotypecaller", mode: 'link', pattern:'*.vcf'
    publishDir "${params.outdir}/3_variants/haplotypecaller/log", mode: 'link', pattern:'*.{txt,log,sh}'

    script:
    """
    gatk --java-options "-Xmx${task.memory.gigaBytes}g" HaplotypeCaller \
        -R $genome \
        -I $bam \
        -O ${sampleName}_hapcaller.vcf.gz \
        --create-output-variant-index \
        --ploidy 1 \
        &> ${sampleName}.hapcaller_variants.log
    gunzip -c ${sampleName}_hapcaller.vcf.gz > ${sampleName}_hapcaller.vcf

    bcftools stats ${sampleName}_hapcaller.vcf.gz > ${sampleName}_hapcaller.bcftools_stats.txt

    cp .command.sh ${sampleName}.hapcaller_variants.sh
    """
}



process HAPLOTYPECALLER_CONSENSUS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path (bam), path (bai), path (vcf), path (tbi)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_varscan2.consensus.masked.fa"), emit: VARSCAN2_CONSENSUS_out
    tuple val(sampleName), path ("${sampleName}_varscan2.consensus.masked_Nremoved.fa")
    path "*.{log,sh}"

    publishDir "${params.outdir}/4_consensus/varscan2", mode: 'link', pattern:'*.fa'
    publishDir "${params.outdir}/4_consensus/varscan2/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    bedtools genomecov \
        -bga \
        -ibam $bam \
        -g $genome \
        | awk '\$4 < 10' | bedtools merge > ${sampleName}.pre_mask.bed


    bcftools view -i '2*FORMAT/AD[0] >= FORMAT/DP[0]' -O z -o ${sampleName}_variantsForConsensus.vcf.gz $vcf
    tabix -p vcf -f ${sampleName}_variantsForConsensus.vcf.gz

    parse_mask_bed.py ${sampleName}_variantsForConsensus.vcf.gz ${sampleName}.pre_mask.bed ${sampleName}.mask.bed
    
    bedtools maskfasta \
        -fi $genome \
        -bed ${sampleName}.mask.bed \
        -fo ${sampleName}_varscan2.ref.masked.fa

    cat ${sampleName}_varscan2.ref.masked.fa \
        | bcftools consensus ${sampleName}_variantsForConsensus.vcf.gz > ${sampleName}_varscan2.consensus.masked.fa

    header=\$(head -n 1 ${sampleName}_varscan2.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_varscan2_masked/g" ${sampleName}_varscan2.consensus.masked.fa

    sed -r '2s/^N{1,}//g' ${sampleName}_varscan2.consensus.masked.fa \
            | sed -r '\$ s/N{1,}\$//g' \
            > ${sampleName}_varscan2.consensus.masked_Nremoved.fa

    cp .command.sh ${sampleName}.varscan2_consensus.sh
    cp .command.log ${sampleName}.varscan2_consensus.log
    """
}
