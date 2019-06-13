library(ggplot2)
library(ggrepel)
library(data.table)
library(ggpubr)
library(rtracklayer)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]
  
#initial setup
settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)

sampleSheet$guide_count <- lengths(strsplit(sampleSheet$sgRNA_ids, ":"))
sampleSheet$guide_count <- ifelse(sampleSheet$sgRNA_ids == 'none', 0, sampleSheet$guide_count)

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

# function to import insertedSequences.tsv files and do some modifications on the resulting table
getInsertions <- function(pipeline_output_dir, samples) {
  insertions <- do.call(rbind, sapply(simplify = FALSE, samples, function(s) {
    f <- file.path(pipeline_output_dir, 
                   'indels',
                   s, 
                   paste0(s, ".insertedSequences.tsv"))
    if(file.exists(f)) {
      dt <- data.table::fread(f)
      dt$sample <- s
      return(dt)
    } else {
      stop("Can't open .insertedSequences.tsv file for sample",s,
           "at",f,"\n")
    }}))
  insertions$end <- insertions$start
  
  #collapse insertions
  insertions <- insertions[,length(name),
                           by = c('seqname', 'sample', 'start', 'end',
                                  'insertedSequence',
                                  'insertionWidth')]
  colnames(insertions)[c(1,7)] <- c('seqnames', 'ReadSupport')
  
  return(insertions)
}

#find cut sites overlapping with the indels
# indels: a data.table object with minimal columns: start, end, 
# cutSites: a GRanges object of cut site coordinates 
# return: data.frame (nrow = nrow(indels), columns are sgRNA ids, 
#         values are 1 if indel overlaps cutsite, otherwise 0. 
overlapCutSites <- function(indels, cutSites, extend = 5) {
  cutSites_ext <- flank(cutSites, width = extend, both = TRUE)
  #check if indel overlaps with the cut site
  query <- GenomicRanges::makeGRangesFromDataFrame(indels)
  overlaps <- as.data.table(findOverlaps(query, cutSites_ext, type = 'any', ignore.strand = TRUE))
  
  M <- matrix(data = rep(0, nrow(indels) * length(cutSites)), 
              nrow = nrow(indels), ncol = length(cutSites))
  colnames(M) <- cutSites$name
  
  M[as.matrix(overlaps)] <- 1
  
  return(M)
}

# remove some problematic samples 
removedSamples <- c( 
  # remove any non-F2 samples
  grep('gen_.*?_F[1345]', sampleSheet$sample_name), 
  #remove F2 samples at 16C
  grep('gen_16C', sampleSheet$sample_name), 
  #for the gen_plates experiments, keep only those at F2 
  grep('^gen_plates.*plate[35]', sampleSheet$sample_name),
  
  ##remove phen*_24 
  grep('phen.*_24$', sampleSheet$sample_name),
  
  ##remove ald_snb-1_ups_CDS2
  grep('ups_CDS2', sampleSheet$sample_name),
  
  ##remove all RNA samples
  grep('pacbio', sampleSheet$sample_name),
  grep('middle', sampleSheet$sample_name), 
  grep('_3end', sampleSheet$sample_name),
  
  ##related to RNA experiment
  grep('lin41_DNA', sampleSheet$sample_name),
  
  ##remove PCR replicates lin-41 pool3
  grep('PCR2', sampleSheet$sample_name),
  grep('PCR3', sampleSheet$sample_name),
  
  ##remove samples where alignment seems to be problematic (seen in N2), 
  ##and mut4&mut8 where it looks like problems with alignment or reference amplicon wrong.
  grep('lin-41_CDS', sampleSheet$sample_name),
  grep('egl-30', sampleSheet$sample_name),
  grep('let-2_CDS', sampleSheet$sample_name),
  grep('mut4', sampleSheet$sample_name),
  grep('mut8', sampleSheet$sample_name),
  
  ##removed samples which were hand-picked for phenotypes and his-72 timecourse experiment 
  grep('picked', sampleSheet$sample_name),
  grep('timecourse', sampleSheet$sample_name),
  
  ##remove snb-1 "ctrl" strain (is a mutants strain and not N2)
  grep('ctrl_strain', sampleSheet$sample_name),
  
  # remove samples with very high off-target alignments 
  which(sampleSheet$target_name == 'par-2_CDS'), 
  which(sampleSheet$target_name == 'rol-6_ups')
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
ggplot(indelStats[sampleMatchesGuide == TRUE], aes(x = treatment, y = scores)) + 
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

# find out if indels overlap cut sites 
#find overlaps with cut sites (only considering guides used in the corresponding sample)
indels <- cbind(indels, overlapCutSites(indels, cutSites, extend = 5))
indels <- do.call(rbind, lapply(unique(indels$sample), function(sampleName) {
  dt <- indels[sample == sampleName]
  sgRNAs <- sampleGuides[[sampleName]]
  dt$atCutSite <- apply(subset(dt, select = sgRNAs), 1, function(x) sum(x > 0) > 0)
  #if an indel overlaps at least two cut sites, then it is considered a double-cut event
  dt$doubleCutEvent <- apply(subset(dt, select = sgRNAs), 1, function(x) sum(x > 0) > 1)
  return(dt)
}))


# make some plots for deletions

deletions <- indels[indelType == 'D']

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
                                  
  dt <- deletions.melt[ReadSupport > 5 & freq >= freqThreshold]
  
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

# check length distribution of deletions at cut sites for treated samples
ggplot2::ggplot(deletions[atCutSite == TRUE & treatment == 'treated'][ReadSupport > 5 & freq > 10^-5], 
                aes(x = doubleCutEvent, y = log10(indelLength))) + 
  geom_violin() + 
  geom_jitter(aes(color = doubleCutEvent), width = 0.2, alpha = 0.25) + 
  labs(title = paste('Length distribution of single-cut or double-cut deletions'), 
       y = 'log10 Deletion Length')  + 
  theme_classic(base_size = 14) + 
  scale_color_brewer(palette = 'Set1') + 
  scale_fill_brewer(palette = 'Set1')

# facet by treatment or single/double cut deletions 
ggplot2::ggplot(deletions[atCutSite == TRUE][ReadSupport > 5 & freq > 10^-5], 
                aes(x = treatment, y = log10(indelLength))) + 
  geom_violin() + 
  geom_jitter(aes(color = treatment), width = 0.2, alpha = 0.25) + 
  facet_wrap(~ doubleCutEvent, nrow = 1) + 
  labs(title = paste('Length distribution of single-cut or double-cut deletions'), 
       y = 'log10 Deletion Length')  + 
  theme_classic(base_size = 14) + 
  scale_color_brewer(palette = 'Set1') + 
  scale_fill_brewer(palette = 'Set1') 

# Off-target indels 
## count number of off-target indel alignments 
off_target <- do.call(rbind, pbapply::pblapply(split(indels, indels$sample), function(dt) {
  #for each sample, calculate how many of the indels don't overlap the actual target region
  s <- unique(dt$sample)
  target <- as(sampleSheet[sample_name == s]$target_region, "GRanges")
  overlaps <- IRanges::overlapsAny(as(dt, "GRanges"), target)
  #count number of reads on and off-target
  percent_off_target_insertions <- round((1-sum(dt[overlaps][indelType == 'I']$ReadSupport)/sum(dt[indelType == 'I']$ReadSupport))*100, 2)
  percent_off_target_deletions <- round((1-sum(dt[overlaps][indelType == 'D']$ReadSupport)/sum(dt[indelType == 'D']$ReadSupport))*100, 2)
  
  return(data.table("sample" = s, "percent_off_target_insertions" = percent_off_target_insertions, 
                    "percent_off_target_deletions" = percent_off_target_deletions))
}))

off_target$treatment <- ifelse(sampleSheet[match(off_target$sample, sample_name)]$sgRNA_ids == 'none', 
                               'untreated', 'treated')
off_target$target_name <- sampleSheet[match(off_target$sample, sample_name)]$target_name

pdf("percent_off_target_indel_alignments.pdf")
p1 <- ggplot2::ggplot(off_target, aes(x = target_name, y = percent_off_target_insertions)) + 
  geom_boxplot(fill = 'blue') + coord_flip()
p2 <- ggplot2::ggplot(off_target, aes(x = target_name, y = percent_off_target_deletions)) + 
  geom_boxplot(fill = 'red') + coord_flip()

cowplot::plot_grid(p1, p2)

dev.off()

# Number of bases affected by number of deletions on the target region
# considering only single-cut deletions 
subsetRleListByRange <- function(input.rle, input.gr) {
  as.vector(input.rle[[seqnames(input.gr)]])[start(input.gr):end(input.gr)]
}

get_perbase_deletion_diversity <- function(dt, target) {
  #convert deletion table to GRanges
  dt.gr <- as(dt, "GRanges")
  #subset deletions by overlap with the target region
  dt.gr <- IRanges::subsetByOverlaps(dt.gr, target)
  #per base diversity of deletions 
  div <- IRanges::coverage(dt.gr)
  #get a vector of values for each base 
  div <- subsetRleListByRange(div, target)
  div[is.na(div)] <- 0
  return(div)
}

perbase_deletions <- do.call(rbind, pbapply::pblapply(sampleSheet$sample_name, function(s) {
  target <- as(sampleSheet[sample_name == s]$target_region, "GRanges")
  #remove low confidence deletions and only consider those at cutsites
  dt <- deletions[sample == s & atCutSite == TRUE & 
                    ReadSupport > 5 & freq > 10^-5]
  
  results <- data.table("sample" = s, 
                        "bp" = seq(width(target)),
                        "del_single_cut" = rep(0, width(target)), 
                        "del_double_cut" = rep(0, width(target)))
  # number of deletions per base considering only single cut events
  if(nrow(dt[doubleCutEvent == FALSE]) > 0) {
    results$del_single_cut <- get_perbase_deletion_diversity(dt[doubleCutEvent == FALSE], target)
  }
  
  # number of deletions per base considering only double cut events
  if(nrow(dt[doubleCutEvent == TRUE]) > 0) {
    results$del_double_cut <- get_perbase_deletion_diversity(dt[doubleCutEvent == TRUE], target)
  }
  return(results)
}))

perbase_deletions$treatment <- ifelse(sampleSheet[match(perbase_deletions$sample, sample_name)]$sgRNA_ids == 'none', 
                               'untreated', 'treated')
perbase_deletions$target_name <- sampleSheet[match(perbase_deletions$sample, sample_name)]$target_name

pdf("bases_affected_versus_number_of_indels.pdf")
ggplot2::ggplot(perbase_deletions[treatment == 'treated' & del_single_cut > 0], 
                aes(x = del_single_cut)) + 
  geom_histogram(binwidth = 1) 
ggplot2::ggplot(perbase_deletions[treatment == 'treated' & del_double_cut > 0], aes(x = del_double_cut)) + 
  geom_histogram(binwidth = 1) 
dev.off()

# for each guide, look for number of deletions with respect to the guide location 
# extending +/- 1000 bp in each direction 
#inverted version of the sampleGuides object 
#only consider treated samples 
treated_samples <- sampleSheet[sampleSheet$sgRNA_ids != 'none']$sample_name
guideSamples <- do.call(rbind, lapply(treated_samples, function(s) {
  data.frame("sample" = s, "guide" = sampleGuides[[s]], stringsAsFactors = FALSE)
}))

guideSamples <- lapply(split(guideSamples, guideSamples$guide), function(x) {
  unique(x$sample)
})

get_diversity <- function(guideSamples, flanksize, doubleCuts = NULL) {
  #for each guide, find number of deletions overlapping the +/- 500 bp of the cut site
  diversity <- do.call(cbind, pbapply::pblapply(names(guideSamples), function(g) {
    #coordinates of the guide
    gr <- flank(cutSites[cutSites$name == g], flanksize, both = TRUE)
    #get deletions from the samples in which the guide was used 
    #and find those that overlap the specific cut site 
    dt <- deletions[sample %in% guideSamples[[g]] & 
                      ReadSupport > 5 & 
                      freq > 10^-5][get(g) == 1]
    if(!is.null(doubleCuts)) {
      dt <- dt[doubleCutEvent == doubleCuts]
    }
    #now get number of deletions per base around the given cut site
    if(nrow(dt) > 0) {
      dt <- as(dt, "GRanges")
      div <- subsetRleListByRange(input.rle = GenomicRanges::coverage(dt), input.gr = gr)
      div[is.na(div)] <- 0
      # normalize by the number of samples 
      div <- div / length(guideSamples[[g]])
    } else {
      div <- rep(0, flanksize * 2)
    }
    df <- data.frame(div)
    colnames(df) <- g
    return(df)
  }))
  return(diversity)
}

diversity_single <- get_diversity(guideSamples, 250, doubleCuts = FALSE)
diversity_all <- get_diversity(guideSamples, 2000)

pdf("Deletion_diversity_around_cut_sites.pdf")
M <- t(diversity_single)
pheatmap::pheatmap(M[names(sort(apply(M, 1, mean), decreasing = T)),], 
                   cluster_cols = FALSE, cluster_rows = FALSE, 
                   color = colorRampPalette(c("white", "blue", "red"))(max(M)), 
                   fontsize_row = 6,  
                   #filename = 'deletion_diversity.only_single_cuts.pdf', 
                   main = 'Deletion Diversity +/- 250 bp of cut sites\n
                   Only single cuts', 
                   labels_col = c(-250, rep('', 248), 0, rep('', 249), 250), 
                   fontsize_col = 12
                   )
M <- t(diversity_all)
pheatmap::pheatmap(M[names(sort(apply(M, 1, mean), decreasing = T)),], 
                   cluster_cols = FALSE, cluster_rows = FALSE, 
                   color = colorRampPalette(c("white", "blue", "red"))(max(M)), 
                   fontsize_row = 6, 
                   #filename = 'deletion_diversity.single_or_doublecuts.pdf', 
                   main = 'Deletion Diversity +/- 2000 bp of cut sites\n
                   Single or Double Cuts',
                   labels_col = c(-2000, rep('', 1998), 0, rep('', 1999), 2000), 
                   fontsize_col = 12
)
dev.off()


###############################################################################

### make some plots on insertions (including inserted sequences information)

insertions <- getInsertions(settings$`output-dir`, samples = sampleSheet$sample_name)

# borrow coverage information from the indels object

insertions <- merge(insertions, 
                    subset(indels, select = c('seqnames', 'sample', 'start', 'end', 
                                              'coverage', 'atCutSite')), 
      by = c('seqnames', 'sample', 'start', 'end'))

insertions$freq <- insertions$ReadSupport / insertions$coverage
insertions$treatment <- ifelse(sampleSheet[match(insertions$sample, sampleSheet$sample_name),]$sgRNA_ids == 'none', 
                           'untreated', 'treated')

insertions <- cbind(insertions, overlapCutSites(indels = insertions, cutSites = cutSites, extend = 5))

insertions.melt <- reshape2::melt(insertions[atCutSite == TRUE], 
                                  id.vars = which(grepl(pattern = 'sg[0-9]+', colnames(insertions)) == FALSE), 
                                  measure.vars = which(grepl(pattern = 'sg[0-9]+', colnames(insertions)) == TRUE))
#remove the amplicon names from the variable field - TODO: update sampleGuides sgRNA ids to account for this 
insertions.melt$variable <- gsub("^.+?:", '', insertions.melt$variable)

insertions.melt$sampleMatchesGuide <- paste(insertions.melt$sample, insertions.melt$variable, sep = ":") %in% sample_guides


#summarize indel counts by sample and treatment condition
#note: indels from untreated samples are compared to all cut sites, while treated samples are compared for a subset of those cut sites
#      so, it might not be fair to compare these numbers at sample level. It should be more fair to compare at sgRNA level. 
#summarize indel counts at cut sites by guide RNA, treatment condition using different frequency thresholds 
pdf("insertion_counts_persample_perguide_plots.pdf")
insertion_count_plots <- sapply(simplify = FALSE, USE.NAMES = TRUE, X = c(1e-04, 3e-05, 2e-05, 1e-05), 
                                function(freqThreshold) {
                                  
                                  dt <- insertions.melt[ReadSupport > 5 & freq >= freqThreshold]
                                  
                                  dt1 <- dt[sampleMatchesGuide == TRUE, sum(value), by = c('sample', 'treatment')]
                                  
                                  #in case some samples were left out due to freq threshold, fill them up with zero values
                                  left_out <- unique(insertions[sample %in% setdiff(insertions$sample, dt1$sample),c('sample', 'treatment')])
                                  if(nrow(left_out) > 0) {
                                    left_out$V1 <- 0
                                    dt1 <- rbind(dt1, left_out)
                                  }
                                  
                                  p1 <- ggboxplot(dt1, x = "treatment", y = "V1",
                                                  color = "treatment", palette = "jco",
                                                  add = "jitter") + 
                                    geom_text(stat="count", aes(label=paste('n =',..count..)), y = -max(dt1$V1)/50) +
                                    stat_compare_means() +
                                    ggtitle('# unique insertions per sample', 
                                            subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%')) + 
                                    theme(axis.title.y = element_blank())
                                  
                                  dt2 <- dt[sampleMatchesGuide == TRUE, sum(value), by = c('variable', 'treatment')]
                                  
                                  p2 <- ggboxplot(dt2, x = "treatment", y = "V1",
                                                  color = "treatment", palette = "jco",
                                                  add = "jitter") + 
                                    geom_text(stat="count", aes(label=paste('n =',..count..)), y = -max(dt2$V1)/50) +
                                    stat_compare_means() +
                                    ggtitle('# unique insertions per sgRNA', 
                                            subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%')) + 
                                    theme(axis.title.y = element_blank())
                                  
                                  p <- cowplot::plot_grid(p1, p2, labels = c('A', 'B'), ncol = 2)
                                  print(p)
                                  return(p)
                                })
dev.off()


insertion_density_plots <- lapply(X = c(1e-04, 3e-05, 2e-05, 1e-05), 
                                  function(freqThreshold) {
                                    ggplot(insertions[atCutSite == TRUE & ReadSupport > 5 & 
                                                        freq >= freqThreshold], 
                                           aes(x = log10(insertionWidth))) + 
                                      geom_density(aes(fill = treatment), alpha = 0.3) + 
                                      ggtitle('', subtitle = paste('Min Freq =',round(freqThreshold * 100, 3),'%'))
                                  }
)

pdf("insertion_length_density_plot.pdf")
print(cowplot::plot_grid(plotlist = insertion_density_plots, labels = "AUTO"))
dev.off()

# length distribution of insertions

ggplot2::ggplot(insertions[atCutSite == TRUE & freq > 10^-5 & ReadSupport > 5], 
                aes(insertionWi)




