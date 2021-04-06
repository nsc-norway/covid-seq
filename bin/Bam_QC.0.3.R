## Bam-QC
## Nacho garcia 
## FHI 2021

## version 0.2 
# Progress bar added 
# Parallelization with doSNOW
# Scores for the lack of reads/ noise identification 
# Single-core function discontinued 


#Libraries
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("Rsamtools")) BiocManager::install("Rsamtools")
if(!require("ShortRead")) BiocManager::install("ShortRead")
if(!require("Rsubread")) BiocManager::install("Rsubread")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("ggpubr")) install.packages("ggpubr")
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("GenomicAlignments")) BiocManager::install("GenomicAlignments")
if(!require("reshape2")) install.packages("reshape2")
if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")
if(!require("foreach")) install.packages("foreach")
if(!require("parallel")) install.packages("parallel")
if(!require("progress")) install.packages("progress")
if(!require("doSNOW")) install.packages("doSNOW")

if(!require("writexl")) install.packages("writexl")

# Bam-QC Parallel -----------------------------------------------------------

Bam_QC_parallel <- function(input.folder= "/media/nacho/Data/NSC/V5/20210208-S3-FHI7/results/2_bam/",
                 primers.bed="/home/nacho/Documents/Corona/NSC/primers.bed",
                 mutation.list="/home/nacho/Documents/Corona/mutationsforQC.csv",
                 number.to.test=259, #Hotspots from previous experiments
                 basal.noise="/media/nacho/Data/NSC/NoiseCtrl_V5.tsv",
                 cutoff=0.25,
                 indexing=FALSE,
                 cores=10,
                 plot.limit=800){
  
  bam_nsc <- input.folder
  positions.to.test<-read.csv(mutation.list)
  
  positions.to.test<-positions.to.test[c(1:number.to.test), ]  
  nucleotides<-positions.to.test$refNtPosition
  
  try(noise<-read.csv(basal.noise, sep ="\t" ))
  
  colnames(noise)<-c("base","Sample","NoiseNSC","ReadsNSC")
  
  if(!exists("noise")){
    print("To be implemented!")
    
  }
  
  #noise.values<-noise$NoiseNSC[-which(is.na(noise$NoiseNSC))]
  #noise.values<-noise.values[-which(noise.values==0)]
  #cutoff<-median(noise.values)+sigma*sd(noise.values)
  
  
  bam_nsc.files<-list.files(bam_nsc, full.names = TRUE,  pattern = "\\.bam$")
  
  if(indexing)indexBam(bam_nsc.files)
  
  plots<-list()
  
  cores.n<-detectCores()
  if(cores>cores.n) cores<- cores.n -2
  if(cores>length(bam_nsc.files)) cores<-length(bam_nsc.files)
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  gc()
  
  ###Progress bar
  samp.id<-gsub("\\..*","",gsub(".*/","",bam_nsc.files))
  pb <- progress_bar$new(format = "Sample: :samp.pb [:bar] :elapsed | eta: :eta",
                         total = length(bam_nsc.files),width = 80)
  opts <- list(progress = function(n){pb$tick(tokens = list(samp.pb = samp.id[n]))} )
  
  out.par<-foreach(i=1:length(bam_nsc.files), .verbose=FALSE, .packages = c("Rsamtools","GenomicRanges",
                                                                       "GenomicAlignments", "ggplot2",
                                                                       "reshape2"),.options.snow = opts) %dopar%{
  #for(i in 1:length(bam_nsc.files)){
    
    bam<-BamFile(bam_nsc.files[i])
    noise.bam<-as.data.frame(matrix(NA, nrow = length(nucleotides), ncol = 3))
    colnames(noise.bam)<-c("Position", "Sample", "Noise")
    noise.bam$Position<-nucleotides    
    noise.bam$Sample<-gsub("\\..*","",gsub(".*/","",bam_nsc.files[i])) 
    counter<-0
    
    for(j in 1:nrow(noise.bam) ){

      base<-noise.bam$Position[j]  
      df<-as.data.frame(matrix(0, ncol = 5, nrow = 1))
      colnames(df)<-c("A","T","C","G","D")
      
      dummybase <- GRanges(seqinfo(bam)@seqnames, IRanges(base , base))
      dummybase<-stackStringsFromBam(bam, param=dummybase)
      dummybase<-unlist(base::strsplit(as.character(unlist(dummybase)),""))  
      dummybase<-table(dummybase)
      
      try(df$A <-dummybase[which(names(dummybase)=="A")],silent = TRUE)
      try(df$T <-dummybase[which(names(dummybase)=="T")],silent = TRUE)
      try(df$C <- dummybase[which(names(dummybase)=="C")],silent = TRUE)
      try(df$G <-dummybase[which(names(dummybase)=="G")],silent = TRUE)
      try(df$D <-dummybase[which(names(dummybase)=="-")],silent = TRUE)
      if(length(dummybase)>0){
        current.noise<- as.numeric((sum(dummybase)-dummybase[which(dummybase==max(dummybase))])/dummybase[which(dummybase==max(dummybase))])
      }else{
        current.noise<-20
      }
            noise.bam$Noise[j]<-current.noise
      
      if(current.noise>=cutoff & current.noise!=20){
        if(!exists("counter")){
          counter<-1
        }else{
          counter<-counter+1
        }
        
        df.agg<-reshape2::melt(df, id.vars=NULL)
        
        plots[[counter]]<-
          ggplot(df.agg)+
          geom_bar(aes(variable,value),stat = "identity", fill="red")+
          ylab("Reads")+
          xlab("Nt")+
          theme_minimal()+
          ggtitle(paste("ID:",gsub("\\..*","",gsub(".*/","",bam_nsc.files[i]))," Base:", base,sep = "" ))
        
      }
      
    }
    
    # if(!exists("noise.out")){
    #   noise.out<-noise.bam
    # }else{
    #   noise.out<-rbind(noise.out, noise.bam)
    # }
    
    output<-list()
    output[[1]]<-plots
    output[[2]]<-noise.bam
    return(output)
    
  }
  stopCluster(cluster.cores)
  
  plots.parallel<-list()
  counter.plots<-0
  try(rm(out.df.par))
  for (i in 1:length(out.par)) {
    
    for (j in 1:length(out.par[[i]][[1]])) {
      if(length(out.par[[i]][[1]])>0){
      if(!is.null(out.par[[i]][[1]][[j]])){
      counter.plots<-counter.plots+1
      plots.parallel[[counter.plots]]<-out.par[[i]][[1]][[j]]}}
    }
    
    if(!exists("out.df.par")){
      out.df.par<-out.par[[i]][[2]]
    }else{
      out.df.par<-rbind(out.df.par, out.par[[i]][[2]])
    }
  }
  
  noise.out<-out.df.par
  
  summary<-as.data.frame(matrix(NA, nrow=length(bam_nsc.files), ncol = 11))
  colnames(summary)<-c("Sample", "Comment","Positions", "NoiseTotal","NoiseMean","NoiseSD","NoisePrimerSites","NoisePrimerSitesMean","NoisePrimerSitesSD",
                       "NoReads","RatioNoisePrimerSite")
  amplicons<-read.csv(primers.bed, sep="\t", header = FALSE)
  colnames(amplicons)<-c("ID", "Start", "End", "Name", "Pool", "Strand")
  samples<-gsub("\\..*","",gsub(".*/","",bam_nsc.files)) 
  summary$Sample<-samples
  #summary$KLD<-NA
  positions.to.test$Primer<-"NO"
  for(l in 1:nrow(positions.to.test)){
    
    if(length(which(amplicons$Start< positions.to.test$refNtPosition[l] & amplicons$End> positions.to.test$refNtPosition[l]))>0)positions.to.test$Primer[l]<-"YES"
    
  }
  
  control.primer<-table(positions.to.test$Primer)
  p.primer.ctrl<-as.numeric(control.primer[names(control.primer)=="YES"]/sum(control.primer))
  qc.plots<-list()
  noise$NoiseNSC[which(noise$NoiseNSC==0)]<-exp(-10)
  for(k in 1:length(samples)){
    dummy.to.test<-noise.out[noise.out$Sample==samples[k],]
    
    dummy.to.test$Primer<-"NO"
    for(l in 1:nrow(dummy.to.test)){
      
      if(length(which(amplicons$Start< dummy.to.test$Position[l] & amplicons$End> dummy.to.test$Position[l]))>0)dummy.to.test$Primer[l]<-"YES"
      
    }
    
    #Code to find out no reads 
    no.reads<-dummy.to.test$Position[dummy.to.test$Noise==20]
    dummy.to.test$Noise[dummy.to.test$Noise==20]<-NA
    
    summary$NoiseTotal[summary$Sample==samples[k]]<-sum(dummy.to.test$Noise, na.rm = TRUE)
    summary$NoiseMean[summary$Sample==samples[k]]<-mean(dummy.to.test$Noise, na.rm=TRUE)
    summary$NoiseSD[summary$Sample==samples[k]]<-sd(dummy.to.test$Noise, na.rm = TRUE)
    
    summary$NoisePrimerSites[summary$Sample==samples[k]]<-sum(dummy.to.test$Noise[dummy.to.test$Primer=="YES"], na.rm = TRUE)
    summary$NoisePrimerSitesMean[summary$Sample==samples[k]]<-mean(dummy.to.test$Noise[dummy.to.test$Primer=="YES"], na.rm=TRUE)
    summary$NoisePrimerSitesSD[summary$Sample==samples[k]]<-sd(dummy.to.test$Noise[dummy.to.test$Primer=="YES"],na.rm=TRUE)
    
    summary$RatioNoisePrimerSite[summary$Sample==samples[k]]<- mean(dummy.to.test$Noise[dummy.to.test$Primer=="YES"], na.rm=TRUE)/
                                                               mean(dummy.to.test$Noise[dummy.to.test$Primer=="NO"], na.rm=TRUE)
      
    
    dummy.to.test$Noise[dummy.to.test$Noise==0]<-exp(-10)
    
    problem.n<-nrow(dummy.to.test[dummy.to.test$Noise>=cutoff,])
    
    if(problem.n==0 & length(no.reads)==0){
      summary$Comment[k]<-"QC PASSED"
    }else{
      
      if(length(which(dummy.to.test$Noise>=cutoff))>0){
        problem.dummy<-dummy.to.test[dummy.to.test$Noise>=cutoff,]
        problem.dummy$Primer<-"NO"
        
        for(l in 1:nrow(problem.dummy)){
          
          if(length(which(amplicons$Start< problem.dummy$Position[l] & amplicons$End> problem.dummy$Position[l]))>0)problem.dummy$Primer[l]<-"YES"
          
        }
        
        
        try(rm(b.test), silent = TRUE)
        try(b.test<-binom.test(as.numeric(table(problem.dummy$Primer)[names(table(problem.dummy$Primer))=="YES"]),
                           sum(table(problem.dummy$Primer)), p=p.primer.ctrl, alternative = "greater"), silent = TRUE)
        if(!exists("b.test")){
          b.test<-list()
          b.test[[1]]<-1
          names(b.test)<-"p.value"
        }
          
        if(b.test$p.value<0.05){
          summary$Comment[k]<-paste("Noise enriched at primer sites p-value=", signif(b.test$p.value, digits = 3), sep = "")  
        }else{
          summary$Comment[k]<- paste("Noise detected [CO:",cutoff, "]",sep = "")
        }
        summary$Positions[k]<-paste(problem.dummy$Position, collapse = ";")}
      
      if(length(no.reads)>0){
        summary$NoReads[k]<-length(no.reads)
        if(is.na(summary$Comment[k])){
          summary$Comment[k]<- paste("No reads found at", length(no.reads), ":",paste(no.reads, collapse = ";"))
        }else{summary$Comment[k]<- paste(summary$Comment[k], "No reads found at", length(no.reads), ":",paste(no.reads, collapse = ";"))} 
        
      }
    }
    
    if(length(which(is.na(dummy.to.test$Noise)))>0) dummy.to.test<-dummy.to.test[-which(is.na(dummy.to.test$Noise)),]
    if(length(which(is.na(noise$NoiseNSC)))>0) noise<-noise[-which(is.na(noise$NoiseNSC)),]
    
    # noise.KL<-noise$NoiseNSC[which(noise$base %in% positions.to.test$refNtPosition)]
    # if(length(which(noise.KL==exp(-10)))>0) noise.KL[which(noise.KL==exp(-10))]<-0
    # 
    # prob.KL<-dummy.to.test$Noise
    # if(length(which(prob.KL==exp(-10)))>0) prob.KL[which(prob.KL==exp(-10))]<-0
    # 
    # hist.ctr<-hist(noise.KL, seq(0,1, length.out=500), plot = FALSE)
    # hist.prblm<-hist(prob.KL, seq(0,1, length.out=500), plot = FALSE)
    # 
    # summary$KLD[k]<-as.numeric(KL(rbind(hist.ctr$counts,hist.prblm$counts),  est.prob = "empirical"))
    
    qc.plots[[k]]<- ggplot()+
      geom_density(data=noise, aes(log(NoiseNSC)), fill="grey")+
      geom_density(data =dummy.to.test, aes(log(Noise)),alpha=0.2, fill="green" )+
      geom_rug(data =dummy.to.test, aes(log(Noise), col=Primer),alpha=0.4 )+
      scale_color_manual(values=c("red","blue"))+
      theme_minimal()+
      geom_vline(xintercept = -10, linetype="dotted",color = "red")+
      xlim(-11,0)+
      xlab("ln(Noise)")+
      ggtitle(paste("ID:",gsub("\\..*","",gsub(".*/","",bam_nsc.files[k])),sep = "" ))
    
  }
  
  #Plotting and saving
  #Dir create
  if(!dir.exists(paste(input.folder,"/BamQC_Output",sep = ""))) dir.create(paste(input.folder,"/BamQC_Output",sep = ""))
  date<-gsub("-","",Sys.Date())
  
  #Cleaning
  summary$Status<-"FAIL"
  summary$Status[summary$Comment=="QC PASSED"]<-"PASS"
  summary$Comment[grep("d N",summary$Comment)]<- gsub("d N", "d. N", summary$Comment[grep("d N",summary$Comment)])
  names(summary)[names(summary)=="Positions"] <- "PositionsNoise"
  summary$PositionsNoReads<-NA
  summary$PositionsNoReads[grep("No reads", summary$Comment)]<-gsub(".*: ","", summary$Comment[grep("No reads", summary$Comment)])
  summary$Comment[grep("No reads", summary$Comment)]<-gsub(" : .*","", summary$Comment[grep("No reads", summary$Comment)])
  summary$PositionsNoise<-gsub("NA;","",summary$PositionsNoise)
  summary$PositionsNoise<-gsub(";NA$","",summary$PositionsNoise)
  
  
  write.csv(summary, paste(bam_nsc,"/BamQC_Output/","ResultsNoise_",date,".csv", sep=""), row.names = FALSE)
  write_xlsx(summary, paste(bam_nsc,"/BamQC_Output/","ResultsNoise_",date,".xlsx", sep=""))
  
  if(length(qc.plots)<=40 ){
    if(length(qc.plots)>0){
    ggarrange(plotlist =  qc.plots[1:length(qc.plots)], ncol = 4, nrow = 10)
    ggsave(paste(bam_nsc,"/BamQC_Output/","ResultsNoisePlots_",date,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
  }  
  }else{
    plotting<-TRUE
    start<-1
    end<-40
    counter<-0
    
    while(plotting){
      if(end==length(qc.plots)) plotting<-FALSE
      ggarrange(plotlist =  qc.plots[start:end], ncol = 4, nrow = 10)
      start<-end+1
      end<-end+40
      if(end>=length(qc.plots)) end<-length(qc.plots)
      counter<-counter+1
      
      ggsave(paste(bam_nsc,"/BamQC_Output/","ResultsNoisePlots",date,"_",counter,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
      
    }
  }
  plots<-plots.parallel
  if(plot.limit>=length(plots)) plot.limit<-length(plots)
  if(plot.limit!=0) plots<-plots[c(1:plot.limit)]
  if(length(plots)>0){
    if(length(plots)<=40){
      ggarrange(plotlist =  plots[1:length(plots)], ncol = 4, nrow = 10)
      ggsave(paste(bam_nsc,"/BamQC_Output/","ResultsNoiseBarPlots",date,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
      
    }else{
      plotting<-TRUE
      start<-1
      end<-40
      counter<-0
      
      while(plotting){
        if(end==length(plots)) plotting<-FALSE
        ggarrange(plotlist =  plots[start:end], ncol = 4, nrow = 10)
        start<-end+1
        end<-end+40
        if(end>=length(plots)) end<-length(plots)
        counter<-counter+1
        
        ggsave(paste(bam_nsc,"/BamQC_Output/","ResultsNoiseBarPlots",date,"_",counter,".pdf", sep=""), width = 420, height = 600, units = "mm") #A4
        
      }
    }
    
  }
}
