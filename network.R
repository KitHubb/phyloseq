library(phyloseq)
library(ggplot2)
library(vegan)

ps <- readRDS("./ps.rds") 
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x) )

plot_bar(ps.rel, fill="Phylum")
plot_bar(ps, fill="Phylum") + 
  geom_bar(position = "fill", stat='identity')


library(dplyr)
tax_table(ps.rel) %>% data.frame() %>% select(Phylum) %>% View()


sample_data(ps.rel)
?plot_net
plot_net(ps.rel, maxdist = 0.4,  
         color = "body.site",
         shape="subject",
         point_size = 5,
         laymeth = "circle")


plot_net(ps.rel, type = "taxa", 
        distance = "bray", maxdist=0.1, 
        color="Phylum")
?plot_net

ps.rel.g<- tax_glom(ps.rel, taxrank = "Genus")
ps.rel.g.top30 <- tax_
plot_net(ps.rel.g, type = "taxa", 
         distance = "bray", maxdist=0.1, 
         color="Genus")
