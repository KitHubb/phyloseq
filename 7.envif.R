library(phyloseq)
library(ggplot2)
library(vegan)


ps <- readRDS("./ps.rds")
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x) )
# ps.rel.2 <- subset_samples(ps.rel, body.site %in% c("gut", "tongue"))
meta <- sample_data(ps.rel)



# https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#test-of-significance
# https://sw1.github.io/teaching/microbiome.html
# https://medium.com/@saurav12das/cca-plot-using-ggplot2-125159f13bbd
PCoA <- ordinate(ps.rel.2, method = "PCoA", distance = "uunifrac")
scree.cap <- plot_scree(PCoA, "Scree Plot for MCs in Constrained Analysis of Principal Coordinates (CAPSCALE)")
print(scree.cap)

# https://bioinformatics.stackexchange.com/questions/16331/how-to-plot-a-pcoa-biplot-with-otu-loadings-as-arrows
# https://forum.qiime2.org/t/how-to-make-pcoa-biplot-in-r-using-q2-deicode-ordination/8377/4









# PCoA 수행
PCoA <- ordinate(ps.rel, method = "PCoA", distance = "bray")

# PCoA biplot 시각화
df <- PCoA$vectors
df2 <- merge(meta, df, by = "row.names")
# PCoA plot
pcoa_plot <- ggplot(data = df2, aes(x = Axis.1, y = Axis.2, color = body.site)) +
  geom_point(size = 3) +
  geom_segment(data = biplot_data, 
               aes(xend = Axis.1*0.2, yend = Axis.2*0.2, x = 0, y = 0,
               arrow = arrow(length = unit(0.25, "cm")), # 화살표
               colour = "red", alpha = 0.8) +
  
  theme_minimal()
pcoa_plot

# geom_segment(data = vec.sp.df, # only line
#              aes(x = 0, xend = MDS1, y = 0, yend = MDS2) ,      #node=0 / length = NMDS /
#              arrow = arrow(length = unit(0.25, "cm")), # 화살표
#              colour = "grey")+
  
  
# 환경변수를 포함한 biplot
biplot_data <- data.frame(PCoA = df, body.site = sample_data(ps.rel)$body.site)
pcoa_biplot <- pcoa_plot +
  geom_segment(data = biplot_data, aes(xend = PC1*0.2, yend = PC2*0.2, x = 0, y = 0, 
                                       arrow = arrow(length = unit(0.2, "inches"))),
               color = "red", alpha = 0.8) +
  geom_text(data = biplot_data, aes(x = PC1*0.3, y = PC
                                    