library(ggplot2)
library(ggrepel)
library(data.table)
library(plotly)
library(knitr)
library(crosstalk)
library(parallel)
library(ggpubr)
library(ggmosaic)
library(rtracklayer)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]
  
#initial setup
settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)
workdir <- getwd()
cutSites <- rtracklayer::import.bed(con = settings$cutsites)
pipelineOutputDir <- settings$`output-dir`

# remove generation samples - only keep those at 24C at F2
#removedSamples <- c(grep('gen_.*?_F[1345]', sampleSheet$sample_name), 
 #                   grep('gen_16C', sampleSheet$sample_name))
#sampleSheet <- sampleSheet[-removedSamples,]

#remove RNA samples where only a portion of the amplicon was sequenced, which gives wrong values for cut sites at boundaries
removedSamples <- c(grep('middle', sampleSheet$sample_name), 
                    grep('3end', sampleSheet$sample_name))

# remove generation samples - only keep those at F2 (if there is both 24C and 16C samples, keep 
# the one at 24C)
removedSamples <- c(removedSamples, 
                    grep('gen_.*?_F[1345]', sampleSheet$sample_name), # remove any non-F2 samples
                    grep('gen_16C', sampleSheet$sample_name), #remove F2 samples at 16C
                    #for the gen_plates experiments, keeep only those at F2 
                    grep('^gen_plates.*plate[35]', sampleSheet$sample_name)
                    )
# remove selected samples from sample sheet
sampleSheet <- sampleSheet[-unique(removedSamples),]

# match samples to the actual sgRNA guides used for that sample
sampleGuides <- lapply(sampleSheet$sample_name, function(s) {
  target <- sampleSheet[sampleSheet$sample_name == s,]$target_name
  sgRNAs <- unlist(strsplit(x = sampleSheet[sampleSheet$sample_name == s,]$sgRNA_ids, 
                            split = ':'))
  if(sgRNAs[1] == 'none') {
    sgRNAs <- setdiff(unique(unlist(strsplit(x = sampleSheet[sampleSheet$target_name == target,]$sgRNA_ids, split = ':'))), 'none')
  }
  return(sgRNAs)
})
names(sampleGuides) <- as.character(sampleSheet$sample_name)


indelStats <- as.data.table(do.call(rbind, lapply(1:nrow(sampleSheet), function(i) {
  sampleName <- sampleSheet[i, 'sample_name']
  target <- sampleSheet[i, 'target_name']
  
  f <- file.path(pipelineOutputDir, 'indels', sampleName, paste0(sampleName, '.sgRNA_efficiency.tsv'))
  if(file.exists(f)) {
    dt <- data.table::fread(f)
    return(dt)
  } else {
    stop("Can't open sgRNA_efficiency.tsv file for sample ",sampleName,
            " at ",f,"\n")
  }
})))


indelStats$sampleMatchesGuide <- as.factor(apply(indelStats, 1, function(x) {
  s <- as.character(x[['sample']])
  g <- as.character(x[['sgRNA']])
  return(g %in% sampleGuides[[s]])
}))

#find the control/untreated samples
#indelStats$treatment <- ifelse(grepl('N2', indelStats$sample) == TRUE, 'untreated', 'treated')
indelStats$treatment <- ifelse(sampleSheet[match(indelStats$sample, sampleSheet$sample_name),]$sgRNA_ids == 'none', 
                               'untreated', 'treated')

#find which platform was used for each sample
indelStats$tech <- sampleSheet[match(indelStats$sample, sampleSheet$sample_name),]$tech

pdf("indel_efficiencies_at_cut_sites.pdf")
# Compare indel efficiencies at cut sites between treated and untreated samples 
dt <- indelStats[sampleMatchesGuide == TRUE]
ggboxplot(dt, x = "treatment", y = "scores",
          color = "treatment", palette = "jco",
          add = "jitter") + 
  geom_text(stat="count", aes(label=paste('n =',..count..)), y = -1) +
  stat_compare_means() +
  labs(title = paste('Indel efficiencies at cut sites'), 
       y = 'Indel Efficiency (%)')  + 
  facet_grid(~ tech)
dev.off()
