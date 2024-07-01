ps <- readRDS("ps.rds") ###################################################


library(phyloseq)
library(reshape2)
library(plyr)
library(tidyverse)


g_abundance = ps %>% tax_glom(taxrank = "Phylum") %>%  # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)  
head(g_abundance)


df_movpic <- ddply(g_abundance, c("Phylum", "body.site", "subject"), summarise,mean=mean(Abundance),sd=sd(Abundance, na.rm=TRUE),n=length(Abundance),se=sd/sqrt(n))
df_movpic


ggplot(df_movpic, aes(x=Phylum, y=mean, fill=body.site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=position_dodge(0.9),width=0.2)+
  theme_bw() +  ylab("Mean proportion")+
  theme(axis.text.x = element_text(angle=60, size=10, hjust=1, vjust=1.0))+
  facet_grid(~body.site)


ggplot(df_movpic, aes(x=Phylum, y=mean, fill=body.site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=position_dodge(0.9),width=0.2)+
  theme_bw() + 
  theme(title=element_text(size=30), axis.text.x = element_text(angle=60, size=20, hjust=1, vjust=1.0), axis.text.y=element_text(size=20), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30),legend.text=element_text(size=20), strip.text.x=element_text(size=30, face="bold")) +
  ggtitle("") + xlab("") + ylab("Mean proportion")
  




g_t <- subset(g_abundance, body.site=="gut" |  body.site=="tongue")
df_g_t <- ddply(g_t, c("Phylum", "body.site", "subject"), summarise,mean=mean(Abundance),sd=sd(Abundance, na.rm=TRUE),n=length(Abundance),se=sd/sqrt(n))

ggplot(df_g_t, aes(x=Phylum, y=mean, fill=body.site)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = c("red", "blue")) +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=position_dodge(0.9),width=0.2)+
  theme_bw() +  ylab("Mean proportion")+
  theme(axis.text.x = element_text(angle=60, size=10, hjust=1, vjust=1.0))



ggplot(df_g_t, aes(x=Phylum, y=mean, fill=body.site)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~ body.site, scales = "free_x", 1) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=position_dodge(0.9),width=0.2)+
  theme_bw() +  ylab("Mean proportion")+
  coord_flip()+
  theme(axis.text.x = element_text(angle=60, size=10, hjust=1, vjust=1.0))
