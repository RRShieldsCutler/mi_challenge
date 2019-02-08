## Alpha div
library(ggplot2)

alpha <- read.delim('../../data/prok/alpha_diversity/alphadiv.txt', header=T, row.names = 1, check.names = F, sep = '\t')
alpha$sampleid <- rownames(alpha)

map <- read.delim('../../data/metadata_vlavage.txt', header=T, check.names = F)
map$rel.day <- map$exp_day - 1

alpha.meta <- merge(alpha, map, by='sampleid')

ggplot(alpha.meta, aes(x=rel.day, y=faith_pd, group=mouse, color=treatment)) +
  geom_line() +
  theme_classic() + labs(x="experiment day relative to first sensitization", y="Faith's phylogenetic diversity") +
  scale_color_manual(name='Treatment', values=c('red','blue')) +
  theme(axis.text = element_text(color='black'))
