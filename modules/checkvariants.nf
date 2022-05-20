process CHECK_VARIANTS {

    label 'large'

    input:
    path "bam/*"
    path "ivar/*"
    path variantList
    path sampleListFile

    output:
    path "detailed_variant_table.csv", emit: DETAILED_TABLE_out
    path "variant_table.csv"
    path "all.checkvariants.sh"
    path "all.checkvariants.log"
    path "plots/*.pdf"

    publishDir "${params.outdir}/6_VuO_MIK", mode: 'link', pattern:'variant_table.csv'
    publishDir "${params.outdir}/6_VuO_MIK", mode: 'link', pattern:'plots/*.pdf'
    publishDir "${params.outdir}/6_VuO_MIK/log", mode: 'link', pattern:'*.{sh,log}'
    publishDir "${params.outdir}/6_VuO_MIK/log", mode: 'link', pattern:'variants.csv'
    
    script:
    """
    check_variants.py $variantList $sampleListFile bam/ ivar/ 
    plotting.py
    
    cp .command.sh all.checkvariants.sh
    cp .command.log all.checkvariants.log
    """
}
