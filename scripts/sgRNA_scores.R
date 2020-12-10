# make plots correlating sgRNA efficiencies with other metrics for sgRNA 

library(ggplot2)
library(ggrepel)
library(data.table)
library(ggpubr)
library(rtracklayer)
library(yaml)
library(ComplexHeatmap)
library(circlize)

args = commandArgs(trailingOnly=TRUE)


# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]
sgrna_score_file <- args[2]
settings <- yaml::read_yaml(settings_file)

ggplot2::theme_set(theme_pubclean())
# exclude some artifically created samples (e.g. by merging multiple samples into one file)
              
sampleSheet <- data.table::fread(settings$sample_sheet)
excluded <- c(
                # remove samples artifically created by merging multiple treated samples
                grep('lin-41.sg15sg16_sg26sg27_sgPool', sampleSheet$sample_name, value = T),
                # remove any non-F2 samples
                grep('gen_.*?_F[1345]', sampleSheet$sample_name, value = T), 
                #remove F2 samples at 16C
                grep('gen_16C', sampleSheet$sample_name, value = T), 
                #for the gen_plates experiments, keep only those at F2 
                grep('^gen_plates.*plate[35]', sampleSheet$sample_name, value = T),
                
                ##remove all RNA samples
                grep('pacbio', sampleSheet$sample_name, value = T),
                grep('middle', sampleSheet$sample_name, value = T), 
                grep('_3end', sampleSheet$sample_name, value = T),
                
                ##related to RNA experiment
                grep('lin41_DNA', sampleSheet$sample_name, value = T),
                
                ##remove PCR replicates lin-41 pool3
                grep('PCR2', sampleSheet$sample_name, value = T),
                grep('PCR3', sampleSheet$sample_name, value = T),
                
                ##removed samples which were hand-picked for phenotypes and his-72 timecourse experiment 
                grep('picked', sampleSheet$sample_name, value = T),
                grep('timecourse', sampleSheet$sample_name, value = T))

sampleSheet <- sampleSheet[!sample_name %in% excluded]

# get a mapping between guides and samples
sampleGuides <- sapply(simplify = F, unique(sampleSheet$sample_name), function(s) {
  sgRNAs <- unlist(strsplit(x = sampleSheet[sampleSheet$sample_name == s,]$sgRNA_ids, 
                            split = ':'))
  return(sgRNAs)
})

# import sgRNA efficiencies from all processed samples
sgrna_stats <- do.call(rbind, lapply(names(sampleGuides), function(s) {
  g <- sampleGuides[[s]]
  if(g[1] != 'none') {
    dt <- data.table::fread(file.path(settings$`output-dir`, 
                                      'indels', 
                                      s, 
                                      paste0(s, '.sgRNA_efficiency.tsv')))
    # only keep scores for the sgRNAs used in the sample
    dt <- dt[sgRNA %in% g]
    return(dt)
  }
}))

# most sgrnas have multiple efficiency scores because they were used in multiple samples.
# we take the average in such cases. 
# which(table(scores$sgRNA) > 1) 
sgrna_stats <- sgrna_stats[,mean(scores),by = sgRNA]
colnames(sgrna_stats)[2] <- 'sgrna_efficiency'


# import sgrna scores collated from various resources 
# these scores are about the sgrna sequence content and design
# keep the relevant columns 
sgrna_scores_boolean <- data.table::fread(sgrna_score_file)[,c(1, 16:27)] 
sgrna_scores_boolean <- cbind(sgrna_scores_boolean[,1], 
                              apply(sgrna_scores_boolean[,-1], 2, as.logical))
colnames(sgrna_scores_boolean)[1] <- 'sgRNA'
colnames(sgrna_scores_boolean) <- gsub("[\\' ]", "_", colnames(sgrna_scores_boolean))

                                

sgrna_scores_numeric <- data.table::fread(sgrna_score_file)[,c(1, 28:38)] 
sgrna_scores_numeric <- cbind(sgrna_scores_numeric[,1], 
                              apply(sgrna_scores_numeric[,-1], 2, as.numeric))
colnames(sgrna_scores_numeric)[1] <- 'sgRNA'
colnames(sgrna_scores_numeric) <- gsub("[\\' ]", "", colnames(sgrna_scores_numeric))

# now merge the sgRNA efficiency scores we obtained from the pipeline and the metrics obtained from external sources
dt <- merge(sgrna_stats, sgrna_scores_boolean, by = 'sgRNA')
mdt <- melt(dt, measure.vars = 3:ncol(dt))
p <- ggpubr::ggboxplot(mdt, x = 'value', y = 'sgrna_efficiency', add = 'jitter', 
                  facet.by = 'variable', color = 'value') +
  stat_compare_means(method = 'wilcox.test', label.y = 15)
ggsave(filename = 'sgrna_scores.boolean.pdf', plot = p, width = 10, height = 8, units = 'in')


dt <- merge(sgrna_stats, sgrna_scores_numeric, by = 'sgRNA')
dt[is.na(percent_in_injection)]$percent_in_injection <- mean(dt$percent_in_injection, na.rm = T)
mdt <- melt(dt, measure.vars = 3:ncol(dt))
mdt$efficiency_threshold <- "sgrna efficiency >= 0%"
added_row_count <- nrow(mdt[sgrna_efficiency >= 4])
mdt <- rbind(mdt, mdt[sgrna_efficiency >= 4])
mdt[(nrow(mdt)-added_row_count+1):nrow(mdt)]$efficiency_threshold <- 'sgrna efficiency >= 4%'
p <- ggpubr::ggscatter(mdt, x = 'sgrna_efficiency',y = 'value',
                  add = 'reg.line', color = 'efficiency_threshold', 
                  facet.by = 'variable', ncol = 4, scales = 'free') + 
  stat_cor(aes(color = efficiency_threshold), label.x = 4)
ggsave(filename = 'sgrna_scores.numerical.double_regression_line.pdf', plot = p, width = 12, height = 9, units = 'in')

# same as above, with single regression line
p <- ggpubr::ggscatter(mdt, x = 'sgrna_efficiency',y = 'value',
                       add = 'reg.line', color = 'gray',
                       facet.by = 'variable', ncol = 5, scales = 'free') + 
  stat_cor(label.x = 3)
ggsave(filename = 'sgrna_scores.numerical.single_regression_line.pdf', plot = p, width = 15, height = 9, units = 'in')

# random forest model, variable importance calculation
sgRNA_scores <- merge(cbind(sgrna_scores_boolean[,1], 
                            apply(sgrna_scores_boolean[,-1], 2, function(x) {as.factor(as.numeric(x))})), 
                      sgrna_scores_numeric, by = 'sgRNA')

dt <- merge(sgrna_stats, sgRNA_scores, by = 'sgRNA')
dt[is.na(percent_in_injection)]$percent_in_injection <- mean(dt$percent_in_injection, na.rm = T)
df <- data.frame(dt[,-1], row.names = dt$sgRNA)

# one variable has a missing value, just impute that
fit <- ranger::ranger(sgrna_efficiency ~ ., df, num.trees = 1000, importance = 'permutation')
imp <- ranger::importance(fit)
p1 <- ggplot(data.frame('importance' = imp, 'variable' = names(imp)), aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = 'identity') + coord_flip() + labs(x = 'Score', y = 'Variable Importance') +
  theme(text = element_text(size = 14))

# pick top variables, re-build model with top variables
top <- names(sort(imp, decreasing = T)[1:10])
fit2 <- ranger::ranger(sgrna_efficiency ~ ., subset(df, select = c('sgrna_efficiency', top)), 
                       num.trees = 1000, importance = 'permutation')
imp2 <- ranger::importance(fit2)
# show correlation between prediction and actual value
p2 <- ggpubr::ggscatter(data.frame('pred' = ranger::predictions(fit2, df[,-1]), 'actual' = df$sgrna_efficiency), 
                  x = 'actual', y = 'pred', add = 'reg.line', cor.coef = T, cor.coef.coord = c(5, 10), 
                  cor.coef.size = 6) + 
                  labs(title = 'Correlation between actual\nand predicted values\nusing top 10 most important variables')

p <- cowplot::plot_grid(p1, p2)
ggsave(filename = 'sgrna_scores.meta_prediction_randomforest.pdf', plot = p, width = 10, height = 6, units = 'in')



