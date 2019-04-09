library(ggplot2)
library(data.table)
library(rtracklayer)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]

# some necessary functions
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
#initial setup
settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)
pipelineOutputDir <- settings$`output-dir`

#define which samples to keep
sampleSheet <- sampleSheet[sample_name %in% c("lin41_DNA_1516",
                                              "lin41_RNApacbio_L1_1516",
                                              "lin41_RNApacbio_L4_1516",
                                              "gen_24C_F2_lin-41_sg26sg27",
                                              "lin41_RNApacbio_L1_2627",
                                              "lin41_RNApacbio_L4_2627",
                                              "lin41_DNA_pool3",
                                              "lin41_RNApacbio_L1_pool3",
                                              "lin41_RNApacbio_L4_pool3",
                                              "lin41_RNA_L1_1516_3end",
                                              "lin41_RNA_L4_1516_3end")]


#define let-7 binding site coordinates on the genome  
let7_sites <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame('seqname' = 'I', 
                  'start' = c(9335255, 9335208), 
                  'end' = c(9335262,9335214), 'strand' = '-', 
                  row.names = c('let7_site1', 'let7_site2'))
)


makePlots <- function(indels, samples) {
  dt <- indels[sample %in% samples & indelType == 'D']

  #now summarize indels by sample 
  df <- data.table::dcast(dt, indelID ~ sample, value.var = 'freq')
  
  # remove the ones don't don't occur in all samples
  df <- df[apply(df, 1 , function(x) sum(is.na(x)) == 0),]

  df <- merge(df, unique(subset(dt, select = c('indelID', 'let7_site1', 'let7_site2', 'let7_either', 'let7_both'))), by = 'indelID')
  
  df$start <-  dt[match(df$indelID, dt$indelID),]$start
  df$end <- dt[match(df$indelID, dt$indelID),]$end
  df$dna_vs_L1 <- df[,get(samples[['dna']])] / df[,get(samples[['L1']])]
  df$dna_vs_L4 <- df[,get(samples[['dna']])] / df[,get(samples[['L4']])]
  df$L4_vs_L1 <- df[,get(samples[['L4']])] / df[,get(samples[['L1']])]
  
  df$let7_overlap <- rep("no-overlap", nrow(df))
  df[which(df$let7_site1),]$let7_overlap <- 'let7_site1'
  df[which(df$let7_site2),]$let7_overlap <- 'let7_site2'
  df[which(df$let7_both),]$let7_overlap <- 'let7_both'
  
  df_text <- df[,quantile(log2(L4_vs_L1))[[4]]+0.2,by = let7_overlap]
  df_text$N <- table(df$let7_overlap)[df_text$let7_overlap]
  
  p1 <- ggplot2::ggplot(df, aes(x = df$let7_overlap, y = log2(L4_vs_L1))) + 
    geom_boxplot(aes(fill = let7_overlap)) + 
    #geom_text(stat = 'count', aes(label = paste("n = ", ..count..)), y = ceiling(min(log2(df$L4_vs_L1)))) + 
    geom_text(data = df_text, 
              aes(label = paste("n =",N), x = let7_overlap, y = V1)) + 
    ggtitle(label = "Relative Deletion Frequency: L4 vs L1 RNA samples",  
                            subtitle = paste("DNA:",samples[['dna']],
                                             "\nRNA L1:",samples[['L1']],
                                             "\nRNA L4:",samples[['L4']])) +
    labs(x = 'Overlap status with let-7 binding sites') + 
    theme_bw(base_size = 14)
  
  return(list(p1))
}

# BEGIN analysis
indels <- do.call(rbind, getIndels(settings$`output-dir`, sampleSheet$sample_name))
indels$freq <- indels$ReadSupport/indels$coverage
indels$indelID <- paste(indels$seqname, indels$start, indels$end, indels$indelType, sep =  ':')

#find overlaps with let7 binding sites on lin-41.

indels.gr <- GenomicRanges::makeGRangesFromDataFrame(indels)
overlaps <- as.data.table(GenomicRanges::findOverlaps(indels.gr, let7_sites))
indels <- cbind(indels, do.call(cbind, lapply(split(overlaps, overlaps$subjectHits), function(x) {
  res <- data.frame(rep(FALSE, nrow(indels)))
  colnames(res)[1] <- names(let7_sites)[unique(x$subjectHits)]
  res[x$queryHits,1] <- TRUE
  return(res)
})))
#find which ones overlap with either
indels$let7_either <- indels$let7_site1 | indels$let7_site2
#find which ones overlap with both
indels$let7_both <- indels$let7_site1 & indels$let7_site2


sg1516_illumina <- makePlots(indels = indels,
                             samples = list('dna' = 'lin41_DNA_1516',
                                            'L1' = 'lin41_RNA_L1_1516_3end',
                                            'L4' = 'lin41_RNA_L4_1516_3end'))

pdf(file = 'sg1516_illumina.pdf')
for(plot in sg1516_illumina) {
  print(plot)
}
dev.off()


sg1516_pacbio <- makePlots(indels = indels[indelType == 'D'], 
                             samples = list('dna' = 'lin41_DNA_1516', 
                                            'L1' = 'lin41_RNApacbio_L1_1516', 
                                            'L4' = 'lin41_RNApacbio_L4_1516'))

pdf(file = 'sg1516_pacbio.pdf')
for(plot in sg1516_pacbio) {
  print(plot)
}
dev.off()


sgPool3_pacbio <- makePlots(indels = indels, 
                           samples = list('dna' = 'lin41_DNA_pool3', 
                                          'L1' = 'lin41_RNApacbio_L1_pool3', 
                                          'L4' = 'lin41_RNApacbio_L4_pool3'))

pdf(file = 'sgPool3_pacbio.pdf')
for(plot in sgPool3_pacbio) {
  print(plot)
}
dev.off()

sg2627_pacbio <- makePlots(indels = indels, 
                           samples = list('dna' = 'gen_24C_F2_lin-41_sg26sg27',
                                          'L1' = 'lin41_RNApacbio_L1_2627',
                                          'L4' = 'lin41_RNApacbio_L4_2627'))

pdf(file = 'sg2627_pacbio.pdf')
for(plot in sg2627_pacbio) {
  print(plot)
}
dev.off()


