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


# Define some necessary functions
getIndels <- function(pipeline_output_dir, samples) {
  indels <- sapply(simplify = FALSE, samples, function(s) {
    f <- file.path(pipeline_output_dir, 
                   'indels',
                   s, 
                   paste0(s, ".indels.tsv"))
    if(file.exists(f)) {
      dt <- data.table::fread(f)
      dt$sample <- s
      return(dt)
    } else {
      stop("Can't open indels.tsv file for sample",s,
           "at",f,"\n")
    }})
  return(indels)
}

#find cut sites overlapping with the indels
# indels: a data.table object with minimal columns: start, end, 
# cutSites: a GRanges object of cut site coordinates 
# return: data.frame (nrow = nrow(indels), columns are sgRNA ids, 
#         values are 1 if indel overlaps cutsite, otherwise 0. 
overlapCutSites <- function(indels, cutSites, extend = 5) {
  cutSites <- flank(cutSites, width = extend, both = TRUE)
  query.start <- GenomicRanges::makeGRangesFromDataFrame(indels, start.field = 'start', end.field = 'start')
  query.end <- GenomicRanges::makeGRangesFromDataFrame(indels, start.field = 'end', end.field = 'end')
  
  startOverlaps <- as.data.table(findOverlaps(query.start, flank(cutSites, width = extend, both = TRUE), type = 'any', ignore.strand = TRUE))
  endOverlaps <- as.data.table(findOverlaps(query.end, cutSites, type = 'any', ignore.strand = TRUE))
  
  overlaps <- merge(startOverlaps, endOverlaps, by = c('queryHits', 'subjectHits'), all = TRUE)
  
  M <- matrix(data = rep(0, nrow(indels) * length(cutSites)), 
              nrow = nrow(indels), ncol = length(cutSites))
  colnames(M) <- cutSites$name
  
  M[as.matrix(overlaps)] <- 1
  
  return(M)
}


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

# remove pacbio samples 
removedSamples <- c(removedSamples, 
                    grep('pacbio', sampleSheet$sample_name))

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
ggplot(dt, aes(x = treatment, y = scores)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_text(stat="count", aes(label=paste('n =',..count..)), y = -1) +
  geom_jitter(aes(color = treatment), width = 0.2, height = 0, show.legend = FALSE) + 
  labs(title = paste('Indel efficiencies at cut sites'), 
       y = 'Indel Efficiency (%)')  + 
  theme_classic(base_size = 14) + scale_color_brewer(palette = 'Set1')
dev.off()

# prepare indel data 

indels <- do.call(rbind, getIndels(settings$`output-dir`, sampleSheet$sample_name))
indels$freq <- indels$ReadSupport/indels$coverage
indels$indelID <- paste(indels$seqnames, indels$start, indels$end, indels$indelType, sep =  ':')
indels$indelLength <- indels$end - indels$start + 1
indels$treatment <- ifelse(sampleSheet[match(indels$sample, sampleSheet$sample_name),]$sgRNA_ids == 'none', 
                           'untreated', 'treated')

# plot proportion of indels per sample
pdf("proportion_of_indels_per_sample.pdf", height = 20)
ggplot(data = indels[,length(seqnames), by = c('indelType', 'sample')], 
            aes(x = sample, y = V1)) + 
  geom_bar(aes(fill = indelType), stat = 'identity', position = 'fill') + 
  scale_fill_brewer(palette = 'Set1') + 
  labs(y = 'Indel Ratio') + 
  geom_text(aes(x = sample, y = 0.2, label = sample), size = 2, color = 'white') +
  theme(axis.text.y =  element_blank(), legend.position = 'bottom') + coord_flip() 
dev.off()


# make some plots for deletions

deletions <- indels[indelType == 'D']

#find overlaps with cut sites (only considering guides used in the corresponding sample)
deletions <- cbind(deletions, overlapCutSites(deletions, cutSites))
deletions <- do.call(rbind, lapply(unique(deletions$sample), function(sampleName) {
  dt <- deletions[sample == sampleName]
  sgRNAs <- sampleGuides[[sampleName]]
  dt$atCutSite <- apply(subset(dt, select = sgRNAs), 1, function(x) sum(x > 0) > 0)
  return(dt)
}))


deletions.melt <- reshape2::melt(deletions[atCutSite == TRUE], 
                              id.vars = which(grepl(pattern = 'sg[0-9]+', colnames(deletions)) == FALSE), 
                              measure.vars = which(grepl(pattern = 'sg[0-9]+', colnames(deletions)) == TRUE))
sample_guides <- unlist(lapply(names(sampleGuides), function(s) {
  paste(s, sampleGuides[[s]], sep = ":")
}))

deletions.melt$sampleMatchesGuide <- paste(deletions.melt$sample, deletions.melt$variable, sep = ":") %in% sample_guides

#summarize indel counts by sample and treatment condition
#note: indels from untreated samples are compared to all cut sites, while treated samples are compared for a subset of those cut sites
#      so, it might not be fair to compare these numbers at sample level. It should be more fair to compare at sgRNA level. 
#summarize indel counts at cut sites by guide RNA, treatment condition using different frequency thresholds 
pdf("deletion_counts_persample_perguide_plots.pdf")
deletion_counts_plots <- sapply(simplify = FALSE, USE.NAMES = TRUE, X = c(1e-04, 3e-05, 2e-05, 1e-05), 
                                function(freqThreshold) {
                                  
  dt <- deletions.melt[freq >= freqThreshold]
  
  dt1 <- dt[sampleMatchesGuide == TRUE,sum(value), by = c('sample', 'treatment')]
  
  #in case some samples were left out due to freq threshold, fill them up with zero values
  left_out <- unique(deletions[sample %in% setdiff(deletions$sample, dt1$sample),c('sample', 'treatment')])
  if(nrow(left_out) > 0) {
    left_out$V1 <- 0
    dt1 <- rbind(dt1, left_out)
  }
  
  p1 <- ggboxplot(dt1, x = "treatment", y = "V1",
                  color = "treatment", palette = "jco",
                  add = "jitter") + 
    geom_text(stat="count", aes(label=paste('n =',..count..)), y = -max(dt1$V1)/50) +
    stat_compare_means() +
    ggtitle('# deletions\n by sample', 
            subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%')) + 
    theme(axis.title.y = element_blank())
  
  
  dt2 <- dt[sampleMatchesGuide == TRUE, sum(value), by = c('variable', 'treatment')]
  
  p2 <- ggboxplot(dt2, x = "treatment", y = "V1",
                  color = "treatment", palette = "jco",
                  add = "jitter") + 
    geom_text(stat="count", aes(label=paste('n =',..count..)), y = -max(dt2$V1)/50) +
    stat_compare_means() +
    ggtitle('# deletions\n by sgRNA', 
            subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%')) + 
    theme(axis.title.y = element_blank())
  
  dt3 <- dt[sampleMatchesGuide == TRUE, sum(value), by = c('sample', 'variable', 'treatment')]
  p3 <- ggboxplot(dt3, x = "treatment", y = "V1",
                  color = "treatment", palette = "jco",
                  add = "jitter") + 
    geom_text(stat="count", aes(label=paste('n =',..count..)), y = -max(dt3$V1)/50) +
    stat_compare_means() +
    ggtitle('# deletions\n by sgRNA & sample', 
            subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%')) + 
    theme(axis.title.y = element_blank())
  
  p <- cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), ncol = 3)
  print(p)
  return(p)
})
dev.off()






