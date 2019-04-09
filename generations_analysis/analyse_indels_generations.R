library(ggplot2)
library(data.table)
library(rtracklayer)
library(yaml)

args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]

#initial setup
settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)
cutSites <- rtracklayer::import.bed(con = settings$cutsites)
pipelineOutputDir <- settings$`output-dir`
analysisTable <- read.table(file = 'analysis_table.tsv', 
                            header = TRUE, sep = '\t')

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


#exclude samples from sample sheet that are not in analysis table 

sampleSheet <- sampleSheet[sample_name %in% analysisTable$sample]

#define let-7 binding site coordinates on the genome  
let7_sites <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame('seqname' = 'I', 
                  'start' = c(9335255, 9335208), 
                  'end' = c(9335262,9335214), 'strand' = '-', 
                  row.names = c('let7_site1', 'let7_site2'))
)
  

# match samples to the actual sgRNA guides used in that sample
sampleGuides <- lapply(sampleSheet$sample_name, function(s) {
  sgRNAs <- unlist(strsplit(x = sampleSheet[sampleSheet$sample_name == s,]$sgRNA_ids, 
                   split = ':'))
})
names(sampleGuides) <- as.character(sampleSheet$sample_name)

indels <- do.call(rbind, getIndels(settings$`output-dir`, sampleSheet$sample_name))
indels$indelLength <- indels$end - indels$start + 1
indels$freq <- indels$ReadSupport/indels$coverage
indels$indelID <- paste(indels$seqname, indels$start, indels$end, indels$indelType, sep =  ':')

indels <- merge(indels, analysisTable, by = 'sample')

#print plots for each analysis group for each temperature condition
plots <- lapply(unique(as.character(analysisTable$analysisGroup)), function(ag) {
  lapply(unique(as.character(analysisTable$temperature)), function(t) {
    cat(ag, t, "\n")
    image_filename <- paste0(ag,".",t,".generations")
    pdf(paste0(image_filename,".pdf"))
    samples <- as.character(analysisTable[analysisTable$analysisGroup == ag & analysisTable$temperature == t,]$sample)

    #remove deletions with read support < 5 and min freq < 0.01 at F2
    selectedIndels <- indels[analysisGroup == ag & 
                               indelType == 'D' & 
                               ReadSupport >= 5 & 
                               freq >= 0.00001 & 
                               generation == 'F2']$indelID
    
    dt <- indels[indelID %in% selectedIndels & 
                   sample %in% samples & 
                   temperature == t]
    
    #find overlaps with let7 binding sites
    dt.gr <- GenomicRanges::makeGRangesFromDataFrame(dt)
    overlaps <- as.data.table(GenomicRanges::findOverlaps(dt.gr, let7_sites))
    dt <- cbind(dt, do.call(cbind, lapply(split(overlaps, overlaps$subjectHits), function(x) {
      res <- data.frame(rep(FALSE, nrow(dt)))
      colnames(res)[1] <- names(let7_sites)[unique(x$subjectHits)]
      res[x$queryHits,1] <- TRUE
      return(res)
    })))
    #find which ones overlap with either
    dt$let7_either <- dt$let7_site1 | dt$let7_site2
    #find which ones overlap with both
    dt$let7_both <- dt$let7_site1 & dt$let7_site2
    
    #now summarize indels by sample, 
    df <- dcast(dt, indelID ~ generation, value.var = 'freq', fill = 0)
    
    #prepare data frame for a heatmap 
    M <- as.matrix(subset(df, select = c('F2', 'F3', 'F4', 'F5')))
    rownames(M) <- df$indelID
    
    #annotation data frame
    annoDf <- unique(data.frame(subset(dt[indelID %in% rownames(M)], select = c('indelID', 'let7_site1', 
                                               'let7_site2', 'let7_either', 
                                               'let7_both'))))
    rownames(annoDf) <- annoDf$indelID
    annoDf$indelID <- NULL
    
    ## calculate average counts per indel over generations
    df.count <- dcast(dt, indelID ~  generation, value.var = 'ReadSupport', fill = 0)
    annoDf$log10_read_count <- log10(rowMeans(df.count[match(rownames(annoDf), df.count$indelID),2:5])) #log10(annoDf$count)
    
    #fontsize <- ifelse(ceiling(nrow(mat)/50) >= 5, 1, ceiling(nrow(mat)/50))
    print(pheatmap::pheatmap(mat = M, 
                             cluster_cols = FALSE, 
                             cutree_rows = 5,
                             #fontsize_row = fontsize, 
                             cellwidth = 60, 
                             scale = 'row',
                             annotation_row = annoDf * 1,  
                             annotation_legend = FALSE, 
                             show_rownames = FALSE, 
                             main = paste('Deletion frequency over generations\nsample group:',ag,
                                          '\ntemperature:',t), border_color = NA))
    
    rownames(df) <- df$indelID
    df2 <- merge(df, annoDf, by = 'row.names')
    df2 <- melt(df2, id.vars = c('indelID', 'F2', 'F5'), measure.vars = grep('let7', colnames(df2), value = T))
    print(ggplot(df2, 
           aes(y = log2((F5+10^-5)/(F2+10^-5)))) + geom_boxplot(aes(fill = value)) + 
      facet_wrap(~ variable, nrow = 2))

    dev.off()
  })
})

