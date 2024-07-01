

library(phyloseq)
library(ggplot2)
library(vegan)


ps <- readRDS("./ps.rds")
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x) )
ps.rel.2 <- subset_samples(ps.rel, body.site %in% c("gut", "tongue"))
bray_dist <- distance(ps.rel.2, method = "bray")
metadata <- sample_data(ps.rel.2)

set.seed(42)
result <- anosim(bray_dist , metadata$body.site, permutations = 9999)
result
# Call:
#   anosim(x = bray_dist, grouping = metadata$body.site, permutations = 9999) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R:     1 
# Significance: 3e-04 
# 
# Permutation: free
# Number of permutations: 9999


set.seed(42)
result2 <- adonis2(bray_dist ~ metadata$body.site,  permutations = 9999 )
result2
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = bray_dist ~ metadata$body.site, permutations = 9999)
# Df SumOfSqs      R2      F Pr(>F)    
# metadata$body.site  1   3.1095 0.58703 21.322  1e-04 ***
#   Residual           15   2.1875 0.41297                  
# Total              16   5.2969 1.00000                  
# ---

permanova <- adonis2(t(as.matrix(otu_table(ps.rel.2))) ~ body.site,
                    data = data.frame(metadata), permutations=9999, method = "bray")


coef <- coefficients(permanova)["body.site",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")






result2$`Pr(>F)`[1]


set.seed(42)
anoV <- anova(betadisper(bray_dist, metadata$body.site),  permutations = 9999)
anoV
# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq   Mean Sq F value   Pr(>F)   
# Groups     1 0.026279 0.0262786  9.0194 0.008914 **
#   Residuals 15 0.043704 0.0029136                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anoV$`Pr(>F)`[1]
# 0.008913996