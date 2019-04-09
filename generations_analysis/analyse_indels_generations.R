#output folder which contains pigx_crispr results and processed tables
pigxOutDir <- '/data/local/buyar/collaborations/jonathan/results_johnny_student_seminar/pipeline_output/'
#sample_sheet.csv file that was used to produce the pigx output 
sampleSheet <- read.csv('/data/local/buyar/collaborations/jonathan/student_seminar_figures/generations_analysis/sample_sheet.csv', stringsAsFactors = F)
analysisTable <- read.table(file = '/data/local/buyar/collaborations/jonathan/student_seminar_figures/generations_analysis/analysis_table.tsv', 
                            header = TRUE, sep = '\t')
cutSites <- read.table('/data/local/buyar/collaborations/jonathan/analyse_generations/cutSites.tsv',
                       stringsAsFactors = TRUE, header = T, sep = '\t')

#define let-7 binding site coordinates on lin-41 amplicon sequence 
let7_sites <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame('seqname' = 'lin-41', 
                  'start' = c(1391, 1439), 
                  'end' = c(1398, 1445), row.names = c('let7_site1', 'let7_site2'))
)
  
library(ggplot2)
library(ggrepel)
library(data.table)
library(plotly)
library(knitr)
library(crosstalk)
library(parallel)

# match samples to the actual sgRNA guides used in that sample
sampleGuides <- lapply(sampleSheet$sample_name, function(s) {
  sgRNAs <- unlist(strsplit(x = sampleSheet[sampleSheet$sample_name == s,]$sgRNA_ids, 
                   split = ':'))
})
names(sampleGuides) <- as.character(sampleSheet$sample_name)

#get coordinates of all detected indels 
indels <- as.data.table(do.call(rbind, lapply(1:nrow(analysisTable), function(i) {
  sampleName <- analysisTable[i, 'sample']
  amplicon <- analysisTable[i, 'amplicon']
  
  f <- file.path(pigxOutDir, 'indels', amplicon, paste0(sampleName, '.indels.unfiltered.tsv'))
  if(file.exists(f)) {
    dt <- data.table::fread(f)
    return(dt)
  } else {
    stop("Can't open indels.unfiltered.tsv file for sample ",sampleName,
            " at ",f,"\n")
  }
})))

#get coverage stats for the samples
coverageStats <- as.data.table(do.call(rbind, lapply(1:nrow(analysisTable), function(i) {
  sampleName <- analysisTable[i, 'sample']
  amplicon <- analysisTable[i, 'amplicon']
  
  f <- file.path(pigxOutDir, 'indels', amplicon, paste0(sampleName, '.coverageStats.tsv'))
  if(file.exists(f)) {
    dt <- data.table::fread(f)
    return(dt)
  } else {
    stop("Can't open coverageStats file for sample ",sampleName,
            " at ",f,"\n")
  }
})))

#collapse indels by coordinates
indels <- indels[,length(readID), by = c('seqname', 'sample', 'start', 'end', 'indelType')]
colnames(indels)[6] <- 'count'

#for each deletion, extract average coverage across the deletion coordinates
cl <- parallel::makeCluster(10)
parallel::clusterExport(cl = cl, varlist = c('indels', 'coverageStats'))
indels$meanCoverage <- pbapply::pbapply(cl = cl, X = indels, MARGIN = 1, FUN = function(x) {
  require(data.table)
  amp <- x['seqname']
  sampleName <- x['sample']
  st <- as.numeric(x['start'])
  end <- as.numeric(x['end'])
  meanCoverage <- mean(coverageStats[seqname == amp & sample == sampleName & bp >= st & bp <= end]$cov)
  return(meanCoverage)
})
parallel::stopCluster(cl)

indels$indelLength <- indels$end - indels$start + 1
indels$freq <- indels$count/indels$meanCoverage
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
                               count >= 5 & 
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
    
    #now summarize indels by generation & freq 
    # df.freq <- dcast(dt, indelID ~ generation, value.var = 'freq', fill = 0)
    # #generation vs count
    # df.count <- dcast(dt, indelID ~ generation, value.var = 'count', fill = 0)
    # 
    # #find deletions that have adequate evidence in any of the generations
    # #filter by both read count and frequency
    # select_freq <- df.freq[apply(df.freq[,2:5], 1, function(x) sum(x > 0.00001) > 0),]$indelID
    # select_count <- df.count[apply(df.count[,2:5], 1, function(x) sum(x > 5) > 0),]$indelID
    # selectedIndels <- base::intersect(select_freq, select_count)
    # 
    # df <- df.freq[df.freq$indelID %in% selectedIndels,]
    
    #now summarize indels by sample, 
    df <- dcast(dt, indelID ~ generation, value.var = 'freq', fill = 0)
    
    #calculate average counts per indel over generations
    df.count <- dcast(dt, indelID ~  generation, value.var = 'count', fill = 0)
    
    
    #prepare data frame for a heatmap 
    M <- as.matrix(subset(df, select = c('F2', 'F3', 'F4', 'F5')))
    rownames(M) <- df$indelID
    
    #annotation data frame
    annoDf <- unique(data.frame(subset(dt[indelID %in% rownames(M)], select = c('indelID', 'let7_site1', 
                                               'let7_site2', 'let7_either', 
                                               'let7_both'))))
    rownames(annoDf) <- annoDf$indelID
    annoDf$indelID <- NULL
    
    annoDf$log10_read_count <- log10(rowMeans(df.count[match(rownames(annoDf), df.count$indelID),2:5])) #log10(annoDf$count)
    
    #fontsize <- ifelse(ceiling(nrow(mat)/50) >= 5, 1, ceiling(nrow(mat)/50))
    print(pheatmap::pheatmap(mat = M, 
                             cluster_cols = FALSE, 
                             cutree_rows = 4,
                             #fontsize_row = fontsize, 
                             cellwidth = 60, 
                             scale = 'row',
                             annotation_row = annoDf * 1,  
                             annotation_legend = FALSE, 
                             show_rownames = FALSE, 
                             main = paste('Deletion frequency over generations\nsample group:',ag,
                                          '\ntemperature:',t), border_color = NA))
    dev.off()
  })
})


