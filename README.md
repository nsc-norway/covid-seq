![NorSeq logo](http://genomics.no/oslo/uploads/images/NorSeqLogo_Acronym%20Colormix.png)

# Covid-seq at [OUS, NSC node](https://www.sequencing.uio.no/) of [NorSeq](https://www.norseq.org/)

## SARS-CoV-2 whole genome sequencing based on multiplexed amplicon method using short-read Illumina sequencers

## Library prep

[Swift Normalase® Amplicon SARS-CoV-2 Panels (SNAP) with Additional Genome Coverage](https://swiftbiosci.com/swif-normalase-amplicon-sars-cov-2-panels/).  

150 bp paired-end sequencing using [Illumina](https://www.illumina.com) [NovaSeq 6000](https://www.illumina.com/systems/sequencing-platforms/novaseq.html) with [NextSeq](https://www.illumina.com/systems/sequencing-platforms/nextseq.html), [MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html) and [HiSeq](https://www.illumina.com/systems/sequencing-platforms/hiseq-2500.html) as backup in the said order.


## Bioinformatics analysis

_This is a snapshot of the production code_

Execution:
```bash
nextflow run main.nf --outpath <Output_folder> --samplelist <SampleList.csv>  --align_tool "bowtie2" -resume
```

References used for the analysis can be found in the folder _util_. See below for details on the sample list.

The pipeline uses singularity or docker containers. To use docker, specify: `--use_docker`. The containers are here: https://github.com/nsc-norway/covid-seq-containers/ . See the file README_script_covid for how to generate the singularity images.


#### In brief:

Primers used in the library prep are trimmed from raw reads using [NSCtrim](https://github.com/nsc-norway/NSCtrim)  
Low quality reads and adapter sequences are trimmed/removed using [fastp](https://github.com/OpenGene/fastp)  
Clean reads are aligned to the genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
Variants and consensus seqeunces are identified using [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html) and [iVar](https://github.com/andersen-lab/ivar)  
Secondary analysis is based on [Pangolin](https://cov-lineages.org/) and [Nextclade](https://clades.nextstrain.org/)  


_Nextflow + Singularity (through Docker) + SLURM executed in linux cluster with 1000+ cores and 5 TB+ Memory_

_Singularity_ image build files will be uploaded soon.


## Sample list

The sample list is a CSV file with one row per sample. Most of the pipeline only uses three columns -- `sample`, `fastq_1` and `fastq_2`.
The report generator and the check_variants_v8.py script use additional information such as the well position.

See [here](doc/sampleList.example.csv) for an example file.

| Column name | Description |
|-------------|-------------|
| `sample`    | Sample name, must be unique. See below for pos/neg controls. |
| `Well`      | Well position, for reports & variants file. Must be in format A1, B1, ... H12. |
| `OrigCtValue`| Ct value for report, or 'NA' if no Ct value is available. |
| `ProjectName`| Project name, for report and QC plots only. |
| `SeqRunId`  | Run ID, for report. |
| `SequencerType` | Sequencer type, for report. |
| `fastq_1`   | Path to read 1 fastq file, relative to where nextflow is started. |
| `fastq_2`   | Path to read 2 fastq file, relative to where nextflow is started. |
| `MIKInputCols` | Optional additional columns to include in the report. Key=Value pairs separated by semicolons. |

If the sample name includes 'neg' / 'NEG' it is treated as a negative control by the report generator and QC plots. Same with 'pos' / 'POS', and positive control.

It is required to have more than one sample passing the alignment QC (BOWTIE2_ALIGN) step. If you only have a single sample, you can disable the
process NSC4FHI_NOISE_NEXTCLADE in `main.nf` by commenting lines 142 to 145 (add //). If you disable that process, the pipeline will run successfully 
with a single sample.

    //NSC4FHI_NOISE_NEXTCLADE(
    //    NOISE_EXTRACTOR.out.NOISE_SUMMARY_FILES_out.collect(),
    //    NEXTCLADE_FOR_FHI.out.NEXTCLADE_FOR_FHI_out
    //    )

