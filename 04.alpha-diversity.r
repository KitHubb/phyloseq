library(phyloseq)
library(ggplot2)
library(vegan)


# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html

ps <- readRDS("./ps.rds") ###################################################




##### 01. basic  #####
?estimate_richness
######################################################################## Not rarefy
summary(sample_sums(ps))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 897    1838    4010    4524    7013    9820 

library(ape)
library(picante)

alpha_df <- estimate_richness(ps,
                              measures = c("Observed", "Chao1", "ACE", "Shannon",
                                           "Simpson", "InvSimpson","Fisher"))
OTU = data.frame(otu_table(ps))

alpha_df$FaithPD <- picante::pd(samp = t(OTU), tree = phy_tree(ps), include.root = F)[,1]
alpha_df

alpha_df %>% head
# Observed Chao1 se.chao1 ACE   se.ACE  Shannon   Simpson InvSimpson    Fisher   FaithPD
# L1S105       63    63        0  63 3.454925 2.682108 0.8707597   7.737527  9.370957  6.455550
# L1S140       65    65        0  65 3.022989 2.660947 0.8518507   6.749945  9.864796  6.258540
# L1S208       85    85        0  85 3.766727 3.121034 0.8999369   9.993691 13.229121  7.764582
# L1S257       81    81        0  81 3.813556 3.262504 0.9261295  13.537212 13.078800  7.560627
# L1S281       72    72        0  72 3.162278 3.189387 0.9082814  10.902920 11.295130  6.379477
# L1S57        70    70        0  70 3.044902 2.905920 0.8661380   7.470381 10.400133  5.943171



meta <- sample_data(ps)
alpha_df2 <- merge(meta, alpha_df, by = "row.names")

ggplot(alpha_df2, aes(x=body.site, y = Shannon)) +
  geom_boxplot(aes(fill = body.site), alpha = 0.6) + theme_classic() +
  scale_fill_manual(values = c("red", "blue", "orange", "yellow")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))



library(reshape2)
alpha_df3 <- alpha_df2[, c(1,3, 10, 11, 13, 15:19)]
alpha_df3
melt <- reshape2::melt(alpha_df3, c("Row.names","body.site"))


ggplot(melt, aes(x=body.site, y = value)) +
  facet_grid(~variable)+
  geom_boxplot(aes(fill = body.site), alpha = 0.6) + theme_classic() +
  scale_fill_manual(values = c("red", "blue", "orange", "yellow")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))







######################################################################## Rarefy
summary(sample_sums(ps))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 897    1838    4010    4524    7013    9820 

ps.rar <- rarefy_even_depth(ps)
# You set `rngseed` to FALSE. Make sure you've set & recorded
#  the random seed of your session for reproducibility.
# See `?set.seed`
# 
# ...
# 137OTUs were removed because they are no longer 
# present in any sample after random subsampling
# 
# ...

head(sample_sums(ps.rar)) 
# L1S105 L1S140 L1S208 L1S257 L1S281  L1S57 
#    897    897    897    897    897    897






alpha_df.rar <- estimate_richness(ps.rar,
                              measures = c("Observed", "Chao1", "ACE", "Shannon",
                                           "Simpson", "InvSimpson","Fisher"))
OTU.rar = data.frame(otu_table(ps.rar))

alpha_df.rar$FaithPD <- picante::pd(samp = t(OTU.rar), tree = phy_tree(ps.rar), include.root = F)[,1]
alpha_df.rar


meta.rar <- sample_data(ps.rar)
alpha_df2.rar <- merge(meta.rar, alpha_df.rar, by = "row.names")



library(reshape2)
alpha_df3.rar <- alpha_df2.rar[, c(1,3, 10, 11, 13, 15:19)]
alpha_df3.rar
melt.rar <- reshape2::melt(alpha_df3.rar, c("Row.names","body.site"))


ggplot(melt.rar, aes(x=body.site, y = value)) +
  facet_grid(~variable)+
  geom_boxplot(aes(fill = body.site), alpha = 0.6) + theme_classic() +
  scale_fill_manual(values = c("red", "blue", "orange", "yellow")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))




alpha_df3 %>% summary


# Row.names              body.site    Observed          Chao1             ACE            Shannon         Simpson         InvSimpson         Fisher          FaithPD      
# Length:34          gut       :8   Min.   : 20.00   Min.   : 20.00   Min.   : 20.00   Min.   :1.993   Min.   :0.7331   Min.   : 3.747   Min.   : 3.337   Min.   : 3.553  
# Class :AsIs        left palm :8   1st Qu.: 40.75   1st Qu.: 40.75   1st Qu.: 40.75   1st Qu.:2.472   1st Qu.:0.8462   1st Qu.: 6.506   1st Qu.: 6.264   1st Qu.: 5.457  
# Mode  :character   right palm:9   Median : 53.50   Median : 53.50   Median : 53.50   Median :2.972   Median :0.8911   Median : 9.180   Median : 9.535   Median : 7.130  
#                    tongue    :9   Mean   : 68.44   Mean   : 68.44   Mean   : 68.44   Mean   :2.960   Mean   :0.8844   Mean   :12.349   Mean   :12.085   Mean   : 9.948  
#                                   3rd Qu.: 84.00   3rd Qu.: 84.00   3rd Qu.: 84.00   3rd Qu.:3.247   3rd Qu.:0.9237   3rd Qu.:13.105   3rd Qu.:13.192   3rd Qu.:12.717  
#                                   Max.   :230.00   Max.   :230.00   Max.   :230.00   Max.   :4.474   Max.   :0.9817   Max.   :54.767   Max.   :42.237   Max.   :37.713  
alpha_df3.rar %>% summary
# Row.names              body.site    Observed          Chao1             ACE            Shannon         Simpson         InvSimpson         Fisher          FaithPD      
# Length:34          gut       :8   Min.   : 20.00   Min.   : 20.00   Min.   : 20.29   Min.   :1.988   Min.   :0.7343   Min.   : 3.763   Min.   : 3.627   Min.   : 3.553  
# Class :AsIs        left palm :8   1st Qu.: 36.00   1st Qu.: 37.08   1st Qu.: 37.96   1st Qu.:2.395   1st Qu.:0.8415   1st Qu.: 6.309   1st Qu.: 7.515   1st Qu.: 4.904  
# Mode  :character   right palm:9   Median : 48.50   Median : 54.50   Median : 52.94   Median :2.970   Median :0.8941   Median : 9.459   Median :10.993   Median : 5.806  
#                    tongue    :9   Mean   : 57.71   Mean   : 63.64   Mean   : 64.71   Mean   :2.916   Mean   :0.8835   Mean   :12.159   Mean   :14.948   Mean   : 7.871  
#                                   3rd Qu.: 71.75   3rd Qu.: 80.25   3rd Qu.: 80.99   3rd Qu.:3.230   3rd Qu.:0.9244   3rd Qu.:13.236   3rd Qu.:18.368   3rd Qu.: 9.816  
#                                   Max.   :158.00   Max.   :181.65   Max.   :195.63   Max.   :4.342   Max.   :0.9794   Max.   :48.503   Max.   :55.621   Max.   :20.250  

alpha_df3$rarefy <- "No"
melt<- reshape2::melt(alpha_df3, c("Row.names","body.site", "rarefy"))

alpha_df3.rar$rarefy <- "Yes"
melt.rar <- reshape2::melt(alpha_df3.rar, c("Row.names","body.site", "rarefy"))
melt.rar

# 
# colnames(melt)[4] <- "Not"
# melt
# melt_total <- cbind(melt,`Yes` = melt.rar$value)
# melt_total %>% head


melt_total <- rbind(melt, melt.rar)
melt_total %>% tail
melt_total %>% head


gut <- melt_total[melt_total$body.site == "gut", ]
ggplot(gut, aes(x=rarefy , y = value)) +
  facet_grid(~variable)+
  geom_boxplot(aes(fill = rarefy), alpha = 0.6) + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))


##### correlation ##### 
alpha_df3
pairs(alpha_df3[,c(3:10)])
cor(alpha_df3[,c(3:10)])
# Observed     Chao1       ACE   Shannon   Simpson InvSimpson    Fisher   FaithPD
# Observed   1.0000000 1.0000000 1.0000000 0.8462887 0.5787187  0.6361918 0.9837859 0.9186595
# Chao1      1.0000000 1.0000000 1.0000000 0.8462887 0.5787187  0.6361918 0.9837859 0.9186595
# ACE        1.0000000 1.0000000 1.0000000 0.8462887 0.5787187  0.6361918 0.9837859 0.9186595
# Shannon    0.8462887 0.8462887 0.8462887 1.0000000 0.8708922  0.8495185 0.8961197 0.8082042
# Simpson    0.5787187 0.5787187 0.5787187 0.8708922 1.0000000  0.7405269 0.6345334 0.5694444
# InvSimpson 0.6361918 0.6361918 0.6361918 0.8495185 0.7405269  1.0000000 0.7294284 0.6744517
# Fisher     0.9837859 0.9837859 0.9837859 0.8961197 0.6345334  0.7294284 1.0000000 0.9506943
# FaithPD    0.9186595 0.9186595 0.9186595 0.8082042 0.5694444  0.6744517 0.9506943 1.0000000
##### 01-2 . statistic test #####

alpha_sh <- alpha_df2[, c("body.site", "Shannon")]
sh_GT <- alpha_sh[alpha_sh$body.site %in% c("gut", "tongue"), ]
####################################################################### t.test

t.test(alpha_sh[alpha_sh$body.site %in% c("gut"), "Shannon"], 
       alpha_sh[alpha_sh$body.site %in% c("tongue"), "Shannon"], paired = F)



#          Welch Two Sample t-test
# 
# t = 2.5999, df = 12.121, p-value = 0.02307
# alternative hypothesis: true difference in means is not equal to 0
# 
# 95 percent confidence interval:
#   0.06596118 0.74388176
#
# sample estimates:
#   mean of x mean of y 
#   2.814283  2.409362 

####################################################################### wilcox.test
wilcox.test(alpha_sh[alpha_sh$body.site %in% c("gut"), "Shannon"],
            alpha_sh[alpha_sh$body.site %in% c("tongue"), "Shannon"], paired = F)


# Wilcoxon rank sum exact test
# 
# W = 59, p-value = 0.0274
# alternative hypothesis: true location shift is not equal to 0





####################################################################### kruskal wallis

?kruskal.test
kruskal.test(Shannon~body.site, data = alpha_sh)

#        Kruskal-Wallis rank sum test
# 
# data:  Shannon by body.site
# Kruskal-Wallis chi-squared = 14.336, df = 3, p-value = 0.002482


kruskal.test(Shannon~body.site, data = alpha_sh,)
library(FSA)
# install.packages("rcompanion")
library(rcompanion)
PT <- dunnTest(Shannon ~ body.site, data = alpha_sh, method="bh") # require the FSA package

# Comparison          Z      P.unadj       P.adj
# 1        gut - left palm -1.7071279 0.0877982798 0.131697420
# 2       gut - right palm -0.9988624 0.3178613394 0.381433607
# 3 left palm - right palm  0.7577577 0.4485960474 0.448596047
# 4           gut - tongue  1.7795825 0.0751443146 0.150288629
# 5     left palm - tongue  3.5362026 0.0004059232 0.002435539
# 6    right palm - tongue  2.8639555 0.0041838683 0.012551605


PT2 <- PT$res  
cldList(comparison = PT2$Comparison, p.value = PT2$P.adj, threshold  = 0.05) # require the rcompanion package   
#        Group Letter MonoLetter
# 1       gut     ab         ab
# 2  leftpalm      a         a 
# 3 rightpalm      a         a 
# 4    tongue      b          b
####################################################################### ANOVA

aov_test <- aov(Shannon ~ body.site, data = alpha_df2)
summary(aov_test)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
# body.site    3  5.308  1.7693   6.616 0.00145 **
# Residuals   30  8.023  0.2674                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# post-hoc test
library(agricolae)
hsd_test <- TukeyHSD(aov_test) # require the agricolae package  
# left palm-gut         0.6239103 -0.07917265  1.3269932 0.0962135
# right palm-gut        0.3995687 -0.28370511  1.0828424 0.3992897
# tongue-gut           -0.4049215 -1.08819525  0.2783523 0.3876978
# right palm-left palm -0.2243416 -0.90761537  0.4589322 0.8086121
# tongue-left palm     -1.0288317 -1.71210550 -0.3455580 0.0015937
# tongue-right palm    -0.8044901 -1.46736306 -0.1416172 0.0126417

hsd_res <- HSD.test(aov_test, "reported.antibiotic.usage", group=T)$reported.antibiotic.usage
hsd_res


##### 02. Custom Alpha Diversity Graphic #####
# Store as a new data variable
alphadt = data.frame(pAlpha$data)
# Subset to just the Shannon index
alphadt <- alphadt[alphadt$variable == "Shannon",]
# Order by Days
alphadt <- alphadt[order(alphadt$day),]
alphadt <- alphadt[,-12]
str(alphadt)
# Define the plot


alphadt
ggplot(data = alphadt, 
       mapping = aes(x = day, y = value, color = body.site, shape = reported.antibiotic.usage)) + 
  geom_point(size = 5) + 
  geom_path() +
  geom_point(data = alphadt[alphadt$reported.antibiotic.usage == "Yes", ], 
             size = 8, alpha = 0.35) +
  ylab("Shannon Index") +
  ggtitle("Shannon Index in Moving Pictures Dataset")

