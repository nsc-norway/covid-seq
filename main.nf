nextflow.enable.dsl=2

pipeline_version = "v4"
nf_mod_path = <NEXTFLOW_DIR> + pipeline_version + "/modules"

// **********************************************************************************

ref_file = <NEXTFLOW_DIR> + pipeline_version + "/util/NC_045512.2.fasta"
primer_bed = <NEXTFLOW_DIR> + pipeline_version + "/util/swift_primers.bed"
primer_master_file = <NEXTFLOW_DIR> + pipeline_version + "/util/sarscov2_v2_masterfile.txt"
vars_under_obs_file = <NEXTFLOW_DIR> + pipeline_version + "/util/variants.csv"

params.ref_id = "NC_045512.2"
params.trim_tool  = "primerclip"
params.align_tool = "bowtie2"
params.outdir = params.outpath + "/results/"

// **********************************************************************************

if (params.outpath.matches('(.*)MIK(.*)')) {
    params.lab = 'MIK'
}
if (params.outpath.matches('(.*)FHI(.*)')) {
    params.lab = 'FHI'
}

File pipeline_tool_file = new File("$params.outpath/pipeline_info.txt")
pipeline_tool_file.write '\n' +
                         'Pipeline:\t' + pipeline_version + '\n' +
                         '\n' +
                         'RunFolder:\t' + params.outpath + '\n' +
                         'SampleSheet:\t' + params.samplelist + '\n' +
                         'Lab:\t\t' + params.lab + '\n' + 
                         'Align:\t\t' + params.align_tool + '\n' +
                         'Trim:\t\t' + params.trim_tool + '\n' +
                         '\n'

include { FASTQC } from "$nf_mod_path/fastqc.nf"
include { FASTP } from "$nf_mod_path/fastp.nf"
include { FASTQC as FASTQC_CLEAN } from "$nf_mod_path/fastqc.nf"

include { MULTIQC } from "$nf_mod_path/multiqc.nf"

include { BOWTIE2_INDEX; BOWTIE2_ALIGN } from "$nf_mod_path/bowtie2.nf"
include { BWA_INDEX; BWA_ALIGN } from "$nf_mod_path/bwa.nf"
include { TANOTI } from "$nf_mod_path/tanoti.nf"

include { IVAR_TRIM } from "$nf_mod_path/ivar.nf"
include { PRIMERCLIP } from "$nf_mod_path/primerclip.nf"

include { PICARD_WGSMETRICS } from "$nf_mod_path/picard.nf"

include { SAMTOOLS_MPILEUP } from "$nf_mod_path/samtools.nf"

include { IVAR_VARIANTS; IVAR_CONSENSUS } from "$nf_mod_path/ivar.nf"
include { VARSCAN2_VARIANTS; VARSCAN2_CONSENSUS } from "$nf_mod_path/varscan2.nf"
include { BCFTOOLS_VARIANTS; BCFTOOLS_CONSENSUS } from "$nf_mod_path/bcftools.nf"

include { PANGOLIN as PANGOLIN_IVAR } from "$nf_mod_path/lineage.nf"
include { PANGOLIN as PANGOLIN_VARSCAN2 } from "$nf_mod_path/lineage.nf"
include { PANGOLIN as PANGOLIN_BCFTOOLS } from "$nf_mod_path/lineage.nf"

include { NEXTCLADE as NEXTCLADE_IVAR } from "$nf_mod_path/lineage.nf"
include { NEXTCLADE as NEXTCLADE_VARSCAN2 } from "$nf_mod_path/lineage.nf"
include { NEXTCLADE as NEXTCLADE_BCFTOOLS } from "$nf_mod_path/lineage.nf"

include { CHECK_VARIANTS } from "$nf_mod_path/checkvariants.nf"

workflow {
    main:
    reads = channel
                .fromPath(params.samplelist)
                .splitCsv(header:true, sep:",")
                .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

    FASTQC(reads, 'raw')
    FASTP(reads)
    FASTQC_CLEAN(FASTP.out.FASTP_out, 'clean')

    if ( params.align_tool == "bowtie2") {
        BOWTIE2_INDEX(ref_file)
        BOWTIE2_ALIGN(FASTP.out.FASTP_out, ref_file, BOWTIE2_INDEX.out.BOWTIE2_INDEX_out)
        ALIGNED = BOWTIE2_ALIGN.out.BOWTIE2_ALIGN_out
    } else if ( params.align_tool == "bwa") {
        BWA_INDEX(ref_file)
        BWA_ALIGN(FASTP.out.FASTP_out, ref_file, BWA_INDEX.out.BWA_INDEX_out)
        ALIGNED = BWA_ALIGN.out.BWA_ALIGN_out
    }    
    
//    if ( params.lab == "FHI") {
//        TANOTI(FASTP.out.FASTP_out, ref_file)
//    }
 
    if ( params.trim_tool == "ivar") {
        IVAR_TRIM(ALIGNED, primer_bed)
        TRIMMED = IVAR_TRIM.out.IVAR_TRIM_out
    } else if ( params.trim_tool == "primerclip") {
        PRIMERCLIP(ALIGNED, primer_master_file)
        TRIMMED = PRIMERCLIP.out.PRIMERCLIP_out
    }	

    PICARD_WGSMETRICS(TRIMMED, ref_file)
    SAMTOOLS_MPILEUP(TRIMMED, ref_file)

    IVAR_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    IVAR_CONSENSUS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    PANGOLIN_IVAR(IVAR_CONSENSUS.out.FOR_LINEAGE_out, 'ivar')
    NEXTCLADE_IVAR(IVAR_CONSENSUS.out.FOR_LINEAGE_out, 'ivar')

    VARSCAN2_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)   
    VARSCAN2_CONSENSUS(TRIMMED.join(VARSCAN2_VARIANTS.out.VARSCAN2_VARIANTS_out), ref_file)
    PANGOLIN_VARSCAN2(VARSCAN2_CONSENSUS.out.FOR_LINEAGE_out, 'varscan2')
    NEXTCLADE_VARSCAN2(VARSCAN2_CONSENSUS.out.FOR_LINEAGE_out, 'varscan2')

    BCFTOOLS_VARIANTS(TRIMMED, ref_file)
    BCFTOOLS_CONSENSUS(TRIMMED.join(BCFTOOLS_VARIANTS.out.BCFTOOLS_VARIANTS_out), ref_file)
    PANGOLIN_BCFTOOLS(BCFTOOLS_CONSENSUS.out.FOR_LINEAGE_out, 'bcftools')
    NEXTCLADE_BCFTOOLS(BCFTOOLS_CONSENSUS.out.FOR_LINEAGE_out, 'bcftools')

    CHECK_VARIANTS(
        TRIMMED.collect { it[1..2] },
        IVAR_VARIANTS.out.IVAR_VARIANTS_out.collect { it[1..2] },
        VARSCAN2_VARIANTS.out.VARSCAN2_VARIANTS_out.collect { it[1..2] },
        BCFTOOLS_VARIANTS.out.BCFTOOLS_VARIANTS_out.collect { it[1..2] },
        vars_under_obs_file,
        Channel.fromPath(params.samplelist)
    )

    // MultiQC
    FILES_FOR_MULTIQC = FASTQC_CLEAN.out.FASTQC_out.collect { it[1] }.mix(
        FASTP.out.FASTP_out_forMULTIQC.collect { it[1] }.mix(
            FASTQC.out.FASTQC_out.collect { it[1] }
        )
    ).collect()
    MULTIQC(FILES_FOR_MULTIQC)

}

