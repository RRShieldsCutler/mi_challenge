library(ggplot2)


unifrac <- read.delim('../../data/prok/beta_diversity/unweighted_unifrac_distance_matrix.txt', header = T, row.names = 1, check.names = F)
map <- read.delim('../../data/metadata_vlavage.txt', header=T, check.names = F)
mouse.ids <- as.character(unique(map$mouse))
mouse.x <-sort(as.character(map[map$mouse==mouse.ids[1],]$sampleid))
self.distance <- data.frame(matrix(nrow = length(mouse.x), ncol = 4))
self.distance$X1