ps.taxa.rel <- readRDS("ps_taxa_rel.rds") ###################################################
# relab_genera = transform_sample_counts(physeq1, function(x) x / sum(x) * 100) 

library(phyloseq)
library(vegan)
library(tidyverse)

ps <- readRDS("./ps.rds") ###################################################
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x) )

# 
ps.taxa.rel
#### 01. PERMANOVA ####
# https://scienceparkstudygroup.github.io/microbiome-lesson/06-beta-diversity/index.html

set.seed(1782) # set seed for analysis reproducibility

#ps.taxa.rel.filt = filter_taxa(ps.taxa.rel, function(x) sum(x > 2) > (0.11 * length(x)), TRUE) 

# https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.taxa.rel))
# Calculate bray curtis distance matrix
bray <- phyloseq::distance(ps.taxa.rel, method = "bray")

# test 
permanova <- adonis2(bray~body.site, data=sampledf, permutations=9999, method = "bray")
'Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 9999

adonis2(formula = bray2 ~ body.site, data = sampledf, permutations = 9999, method = "bray")
          Df SumOfSqs      R2      F Pr(>F)    
body.site  3   5.2330 0.42437 7.3723  1e-04 ***
Residual  30   7.0982 0.57563                  
Total     33  12.3312 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1'

permanova2 <- adonis2(bray~sampledf$subject, data=sampledf, permutations=9999, method = "bray")
permanova2


permanova3 <- adonis2(bray~sampledf$reported.antibiotic.usage, data=sampledf, permutations=9999, method = "bray")
permanova3


permanova4 <- adonis2(bray~sampledf$body.site*sampledf$year, data=sampledf, permutations=9999, method = "bray")
permanova4

permanova5 <- adonis2(bray~sampledf, data=sampledf, permutations=9999, method = "bray")
permanova5

#### 02. ANOVA ####
otu = t(otu.filt)
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
anova_result <- aov(otu ~ body.site, samp)
anova_result
summary(anova_result)


bray <- phyloseq::distance(ps.taxa.rel, method = "bray")
sampledf <- data.frame(sample_data(ps.taxa.rel))
Anova <- anova(betadisper(bray, sampledf$body.site))
?betadisper
'Analysis of Variance Table

Response: Distances
Df  Sum Sq  Mean Sq F value    Pr(>F)    
Groups     3 0.28883 0.096278  9.0096 0.0002083 ***
  Residuals 30 0.32059 0.010686                      
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1'

Anova



#### 03. Tukey’s Honest Significant Difference (HSD) test XXX  #### 사후검정 
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html


#install.packages("agricolae")
library(agricolae)

tukey_result <- HSD.test(anova_result, "body.site", group = TRUE)
'Error in quantile.default(x) : 
  만약 'na.rm'이 FALSE라면 결측값들과 NaN은 허용되지 않습니다
In addition: Warning message:
In tapply.stat(junto[, 1], junto[, 2], stat = "median") :
  NAs introduced by coercion'

print(tukey_result)




#### 03. wilcox ####


#### 03. ####
#### 03. ####
#### 03. ####
#### 03. ####

#### +)  PCoA ####

bray = ordinate(ps.taxa.rel, method = "PCoA", distance = 'bray')

plot_ordination(ps.taxa.rel, bray, color = "body.site") + 
  geom_point(mapping = aes(size = day, 
                           shape = factor(subject))) +
  ggtitle("PCoA: Bray-Curtis")

##### Bray-Curtis PCoA Axis-1 vs. Day ##### 
# Join sample data and ordination axes together in one data.table
bray_dt = data.frame(bray$vectors, keep.rownames = TRUE)

df.1 = cbind(bray_dt[1:3], sample_data(ps.taxa.rel))
df.1

# Axis 1
ggplot(df.1, aes(x = day,  y = Axis.1, color = factor(subject))) +
  geom_point(size = 5) + 
  geom_path() +
  facet_wrap(~body.site) + ggtitle("Bray-Curtis PCoA Axis-1 vs. Day")

# Axis 2
ggplot(df.1, aes(x = day,  y = Axis.2, color = factor(subject))) +
  geom_point(size = 5) + 
  geom_path() +
  facet_wrap(~body.site) + ggtitle("Bray-Curtis PCoA Axis-1 vs. Day")


p.bray = ordinate(ps.rel, method = "PCoA", distance = 'bray')

plot_ordination(ps.rel, p.bray, color = "body.site") + 
  geom_point(mapping = aes(size = day, 
                           shape = factor(subject))) +
  theme_bw() + 
  ggtitle("PCoA: Bray-Curtis")


plot_ordination(ps.rel, p.bray, type="samples", color="body.site", shape = "body.site") +
  geom_point(size=3) + 
  stat_ellipse() + theme_bw() + 
  ggtitle("PCoA: Bray-Curtis")



p.jaccard = ordinate(ps.rel, method = "PCoA", distance = 'jaccard')
plot_ordination(ps.rel, p.jaccard, color = "body.site")


?ordinate
N.bray = ordinate(ps.rel, method = "NMDS", distance = 'bray')
plot_ordination(ps.rel, N.bray, color = "body.site")




