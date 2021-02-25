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

References used for the analysis can be found in the folder _util_.
  
#### In brief:

Primers used in the library prep are trimmed from raw reads using [pTrimmer](https://github.com/DMU-lilab/pTrimmer)  
Low quality reads and adapter sequences are trimmed/removed using [fastp](https://github.com/OpenGene/fastp)  
Clean reads are aligned to the genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
Variants and consensus seqeunces are identified using [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html) and [iVar](https://github.com/andersen-lab/ivar)  
Secondary analysis is based on [Pangolin](https://cov-lineages.org/) and [Nextclade](https://clades.nextstrain.org/)  


_Nextflow + Singularity (through Docker) + SLURM executed in linux cluster with 1000+ cores and 5 TB+ Memory_

_Singularity_ image build files will be uploaded soon.
