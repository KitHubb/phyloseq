# https://www.nature.com/articles/s41467-022-28034-z
# https://microbiome.github.io/OMA/differential-abundance.html
# https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html
# https://www.nicholas-ollberding.com/post/identifying-differentially-abundant-features-in-microbiome-data/
# https://www.yanh.org/2021/01/01/microbiome-r/ 




#### package #### 
getwd() # "D:/KSY/Study/Microbiome/R/phyloseq"

p1 <- c("tidyverse", "vegan", "BiocManager", "RColorBrewer", "remotes")
p2 <- c("phyloseq", "mia", "ANCOMBC", "DESeq2", "ComplexHeatmap", "ALDEx2", "Maaslin2" , "microbiomeMarker")
load_package <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    ifelse(p %in% p1, 
           install.packages(p, repos = "http://cran.us.r-project.org/"), 
           BiocManager::install(p))
  }
  library(p, character.only = TRUE, quietly = TRUE)
}
invisible(lapply(c(p1,p2), load_package))
# For ANCOMBC we need development version
# Can be installed with:
# remotes::install_github("FrederickHuangLin/ANCOMBC")


#### phyloseq #### 
otu <- read.table(file = "feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Greengenes format
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("tree.nwk")
# merge the data
ps <- phyloseq(OTU, TAX, SAMPLE,TREE)
ps

saveRDS(ps, "./ps.rds")
ps <- readRDS("ps.rds")
ps

ps <- readRDS("ps.rds") # >>>> ps.rds <<<<##################################################
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
saveRDS(ps.taxa.rel, "./ps_taxa_rel.rds")
ps.taxa.rel <-readEDS("./ps_taxa_rel.rd") #>>>> ps.taxa.rel.rds <<<< 




####  Prevalence Filtering #### 
# https://microbiome.github.io/OMA/differential-abundance.html
# https://microbiome.github.io/OMA/differential-abundance.html#ref-Nearing2022
# applying a 10% threshold for the prevalence of the taxa generally resulted in more robust results



install.packages("devtools")
library(devtools)
devtools::install_github("vmikk/metagMisc")
library(metagMisc)

ps.10 <- phyloseq_filter_prevalence(ps, prev.trh = 0.1, abund.type = "total")


saveRDS(ps.10, "./ps_10.rds")
ps.10 <-readRDS("./ps_10.rds")#>>>> ps.10.rds <<<<###############################################


sample_data(ps.10)$body.site <- as.factor(sample_data(ps.10)$body.site) # factorize for DESeq2
ps.10.taxa <- tax_glom(ps.10, taxrank = 'Species', NArm = FALSE)
saveRDS(ps.10.taxa, "./ps_10_taxa.rds")
ps.10.taxa<-readRDS("./ps_10_taxa.rds") # >>>> ps.10.taxa.rds <<<< ##AFTER factor and species level ##########################
ps.10.taxa

ps.10.taxa.rel <- transform_sample_counts(ps.10.taxa, function(x) x/sum(x)*100) 
saveRDS(ps.10.taxa.rel, "./ps_10_taxa_rel.rds")

#### +) 선발 후 glam or glam 후 선발? #### 
ps # 770
ps.taxa <- tax_glom(ps,  taxrank = 'Species', NArm = FALSE) # 351



ps.glom.10 <- phyloseq_filter_prevalence(ps.taxa, prev.trh = 0.1, abund.type = "total") # 148
ps.glom.10 <- mia::subsetByPrevalentTaxa(ps.taxa,  detection = 0, prevalence = 0.1)


p_t = rownames(otu_table(ps.10.taxa)) # 125
t_p = rownames(otu_table(ps.glom.10)) # 148

otp = setdiff( t_p, p_t) # 29 only t_p
opt = setdiff(p_t,  t_p) # 5 only p_t


## 29 ## 
ps.glom.10 %>% otu_table() %>% subset( t_p %in% otp) %>%  apply( 1, sum) %>% summary()
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     10      63      98     408     171    3650

## 5 ##
ps.10.taxa %>% otu_table() %>% subset( p_t %in% opt) %>%  apply( 1, sum) %>% summary()
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   22.0    59.0   111.0   548.7   984.2  1731.0


# OTU아이디로 subset 
#ps.10.taxa %>% otu_table() %>% subset( p_t %in% opt) %>% 
#  merge_phyloseq(tax_table(ps.10.taxa), sample_data(ps.10.taxa),phy_tree(ps.10.taxa)) %>%
#  otu_table() %>% apply( 1, sum)  %>% summary()






#### 1-1) deseq2 ####

# https://www.yanh.org/2021/01/01/microbiome-r/

ps <- ps.10.taxa # after 10% and tax_glom()
otu_table(ps) <- otu_table(ps.10.taxa) + 1

ps.sub <- subset_samples(ps, body.site %in% c("gut", "tongue")) # pairwise comparison between gut and tongue

ds = phyloseq_to_deseq2(ps.sub, ~ body.site)
ds = DESeq(ds, test="Wald",fitType="parametric" )
alpha = 0.05 
res = results(ds, alpha=alpha)

res_j = res[order(res$padj, na.last=NA), ]
#res_p = res[order(res$pvalue, na.last=NA), ]

#taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
#taxa_sig = rownames(res[res$padj > 0.05, ]) # select bottom 20 with lowest p.adj values
#taxa_sig = rownames(res[res$padj<0.05,])
#taxa_sig_p = rownames(res_p[res_p$pvalue<alpha,])
res_sig_j = res_j[res_j$pvalue<alpha,]
'log2 fold change (MLE): body.site tongue vs gut 
Wald test p-value: body.site tongue vs gut 
DataFrame with 61 rows and 6 columns
                                  baseMean log2FoldChange     lfcSE      stat       pvalue         padj
                                 <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
4b5eeb300368260019c1fbc7a3c718fc  1638.529       -8.90983  0.278312  -32.0138 7.00011e-225 5.04008e-223
0305a4993ecf2d8ef4149fdfc7592603   309.235       -9.35755  0.671690  -13.9313  4.08564e-44  1.47083e-42'
taxa_sig_j = rownames(res_sig_j)

#ps.taxa.rel.sig <- prune_taxa(taxa_sig_p, ps.taxa.rel)
#ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.sub)), ps.taxa.rel.sig)# Only keep gut and tongue samples


ps.taxa.rel.sig <- prune_taxa(taxa_sig_j, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.sub)), ps.taxa.rel.sig)# Only keep gut and tongue samples




#### 1-2) deseq2-heatmap #### 
matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
deseq_otu <- rownames(matrix)

# Define the annotation color for columns and rows
annotation_col = data.frame(
  `Body site` = as.factor(metadata_sub$body.site), 
  check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"]) )
rownames(annotation_row) = rownames(matrix)

# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  `Body site` = c(gut = "purple", tongue = "yellow"),
  Phylum = phylum_col)

# heat- map
ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row,
                         annotation_colors = ann_colors)
# save as 1000x650






#### 1-3) deseq2-genus~log2FC-of-the-significant-OTUs ####
# https://micca.readthedocs.io/en/latest/phyloseq.html#otu-differential-abundance-testing-with-deseq2

res_sig_j= cbind(as(res_sig_j, "data.frame"), as(tax_table(ps)[rownames(res_sig_j), ], "matrix")) # 추출 
res_sig_j = data.frame(res_sig_j)
ggplot(res_sig_j, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


#### 2-1)ANCOM-BC -1 #### 
# https://www.yanh.org/2021/01/01/microbiome-r/

# pairwise comparison between gut and tongue
ps.10.taxa<-readRDS("./ps_10_taxa.rds") 

ps <- ps.10.taxa
ps.sub <- subset_samples(ps, body.site %in% c("gut", "tongue")) 

ps
# ancombc
out <- ancombc(phyloseq = ps.sub, 
               formula = "body.site", 
               p_adj_method = "fdr", # same as DESeq, edgeR
               group = "body.site",
               struc_zero = TRUE, #Group에서 Zeros 데이터 detect
               neg_lb = TRUE, # taxon에서 zero 분류 
               conserve = TRUE,  # It is recommended if the sample size is small and/or the number of differentially abundant taxa is believed to be large
               global = TRUE,
               lib_cut = 1000,
               alpha = 0.05, 
               tol = 1e-5
               ) 



#res.or_p <- rownames(res$q_val[,"body.sitetongue"])[base::order(res$q_val[,"body.sitetongue"])]
taxa_sig <- out$res$p_val %>%arrange(body.sitetongue) %>% 
  filter( body.sitetongue <.05 ) %>% rownames()


# select the bottom 20 with lowest p values
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.sub)), ps.taxa.rel.sig)


#### 2-2) ANCOM-BC -1 heatmap #### 
matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig)) 
ancom_bc_otu <- rownames(matrix)
'# for pheatmap
# remove sum = 0
apply(matrix, 1, sum)>0  # Unclassified  Mogibacterium = 0
matrix <- matrix[apply(matrix, 1, sum)>0, ]
# matrix1 <- data_frame[rownames(matrix) %in% rows, ] 
apply(matrix, 1, sum) 

colnames(matrix) # 18
# [1] "L1S105" "L1S140" "L1S208" "L1S257" "L1S281" "L1S57"  "L1S76"  "L1S8"   "L5S104"
# [10] "L5S155" "L5S174" "L5S203" "L5S222" "L5S240" "L6S20"  "L6S68"  "L6S93"
metadata_sub <- metadata_sub[colnames(matrix), ] # 17
rownames(metadata_sub) # OK'



# Define the annotation color for columns and rows
annotation_col = data.frame(
  `Body site` = as.factor(metadata_sub$body.site), 
  check.names = FALSE )
rownames(annotation_col) = rownames(metadata_sub)


'taxx = tax_table(ps.taxa.rel.sig)
taxx = taxx[!(taxx[,"Species"]  == "Unclassified  Mogibacterium"), ]'
annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"]) )
rownames(annotation_row) = rownames(matrix)


# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  `Body site` = c(`gut` = "purple", tongue = "yellow"),
  Phylum = phylum_col )

ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors)

# Error in hclust(get_dist(submat, distance), method = method) : 
# NA/NaN/Inf in foreign function call (arg 10)



#### 3-1) MaAsLin2 ####

#### 3-2) MaAsLin2 heatmap#### 
#### ALDEx2 #### 
#### Wilcoxon #### 
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
# Comparing taxon abundance in two groups

obj$data$tax_abund <- calc_taxon_abund(ps, "otu_rarefied",
                                       cols = sample_data$SampleID)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

set.seed(1)
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs, # number of OTUs
          node_color = log2_median_ratio, # difference between groups
          node_color_interval = c(-10, 10), # symmetric interval
          node_color_range = c("cyan", "gray", "magenta"), # diverging colors
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median counts")

#### MaAsLin2#### 
#### Lefse#### 

# https://rpubs.com/mohsen/lefse_analysis_cleaned_prevalence_phyloseq
# phyloseq MUST be
  # 1. Removing NA phylum sequences.
  # 2. Removing sequences that occur in not more than 3 samples .
    # How to ?
  # 3. Removing the phyla **Cyanobacteria/Chloroplast, Deferribacteres, and Tenericutes** which are comprised of low-prevalence features.



# Tutorial # 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/microbiomeMarker/inst/doc/microbiomeMarker-vignette.html
# https://rpubs.com/mohsen/lefse_analysis_cleaned_prevalence_phyloseq


install.packages("rlang")
install.packages('rlang')
library(rlang)
sessionInfo() # rlang_1.0.5 

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("microbiomeMarker")
library(microbiomeMarker) #

## Tutorial #### 
data(kostic_crc)
kostic_crc
mm_test <- run_lefse( kostic_crc,
                      wilcoxon_cutoff = 0.01,
                      norm = "CPM", 
                      group = "DIAGNOSIS",
                      kw_cutoff = 0.01,
                      multigrp_strat = TRUE,
                      lda_cutoff = 4
  )
# if error  =>  '.rs.restartR()'

plot_ef_bar(mm_test)

mm_lefse <- run_lefse(
  kostic_crc,
  wilcoxon_cutoff = 0.01,
  group = "DIAGNOSIS",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)


samp <- sample_data(mm_lefse)
str(samp$DIAGNOSIS)


plot_cladogram(mm_lefse, color = c(Healthy = "darkgreen", Tumor = "red")) + 
  theme(plot.margin = margin(0, 0, 0, 0))

plot_ef_bar(mm_lefse)

plot_ef_dot(mm_lefse)

head(marker_table(mm_lefse))

#### phyloseq ####

ps.10.taxa # 10% prevalenced 
ps.10.taxa.rel # relative abundance
sample_data(ps.10.taxa.rel)

ps.10.taxa.sub<-phyloseq::subset_samples(ps.10.taxa, body.site %in% c("gut","tongue"))
# ound more than one class "phylo" in cache; using the first, from namespace 'phyloseq' Also defined by ‘tidytree
# -> ignore this massege 


ps.10.taxa.sub

# run_lefse 
?run_lefse
lef_out_01_cut2 <- run_lefse(ps.10.taxa.sub, group = "body.site", 
                   norm = "CPM", # normalization 
                   # CPM  : pre-sample normalization of the sum of the values to 1e+06
                   # 만약 rAB면 -> none하면 됨
                   taxa_rank = "all" ,
                   kw_cutoff = 0.05, 
                   wilcoxon_cutoff = 0.05,
                   bootstrap_n = 999,
                   lda_cutoff = 4)

# Visualization 
plot_cladogram(lef_out_01_cut2, color = c("darkgreen", "red")) # 'all' 로 하니까 된다..
plot_ef_bar(lef_out_01_cut2)
plot_ef_dot(lef_out_01_cut2)

plot_sl_roc(lef_out_01_cut2, group = "body.site")
p_pht <- plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")& theme_bw()

# table 
lef_table <- marker_table(lef_out_01_cut2)
lef_table <- lef_table[lef_table$ef_lda >= 4.0 & lef_table$pvalue< 0.05& lef_table$padj< 0.05,1:2 ]
lef_table

ss <- grep('s__', lef_table$feature)
lef_table_ss <- lef_table[ss,]
lef_table_ss
#### corncob #### 
#### limma voom #### 
#### metagenomeSeq #### 
#### t-test #### 
#### LinDA #### 
#### edgeR #### 