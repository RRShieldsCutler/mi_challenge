library("gplots")
library("RColorBrewer")
library(beeswarm)
library(reshape2)
library(vegan)
library(ggplot2)

otu.L3 <- read.delim('../../data/prok/otu_table_L3.txt', header=T, check.names = F, row.names = 1, skip = 1)
map <- read.delim('../../data/metadata_vlavage.txt', header=T, check.names = F)

mouse.ids <- as.character(unique(map$mouse))

rownames(map) <- map$sampleid
samp.ids <- intersect(colnames(otu.L3), rownames(map))
samp.ids <- sort(samp.ids)
otu.L3 <-otu.L3[,samp.ids]
map <- map[samp.ids,]

for (m in 1:length(mouse.ids)) {
  m.samples <- as.character(map[map$mouse==mouse.ids[m],]$sampleid)
  otu.t <- otu.L3[,m.samples]
  map.t <- map[map$mouse==mouse.ids[m],]
  
  dates.ax <- sort(map.t$exp_day)
  
  otu.n <- sweep(otu.t,2,colSums(otu.t),'/');                # Normalize to relative abundance
  
  otu.m <- sweep(sqrt(otu.n), 2, colSums(sqrt(otu.n)), '/')
  meanAb <- apply(otu.m, 1, FUN=function(xx) tapply(xx, map.t$treatment, mean)) # group mean
  
  
  ranked = order(apply(meanAb, 2, max),decreasing=T)
  otu.m = otu.m[ranked, ]
  
  rownames(otu.m) <- lapply(X = rownames(otu.m), FUN = function(xx) strsplit(as.character(xx), ';', fixed = T)[[1]][3])
  Taxa = rownames(otu.m)
  lim = 20
  if (nrow(otu.m) > lim) Taxa[lim:nrow(otu.m)] = "Other"
  otu.m = rowsum(otu.m, Taxa)
  byAbundance = rownames(otu.m)[order(rowMeans(otu.m), decreasing=T)]
  #byAbundance = gsub(';','.',byAbundance)
  
  otu.m <- data.frame(t(otu.m), check.names=F)      # flip table
  otu.m$sampleid <- rownames(otu.m)  # add a column for the sample IDs
  rownames(map.t) <- map.t$sampleid      # add a column for the sample IDs
  
  # The following separates taxa abundances per sample, then splices in the column of interest
  otu.m <- melt(otu.m, id.vars = "sampleid", variable.name = "Taxa", value.name = "RelativeAbundance")
  otu.m <- merge(otu.m, map.t[,c("sampleid","exp_day")], by="sampleid")
  otu.m$Taxa <- factor(otu.m$Taxa, levels=byAbundance, ordered=T)
  otu.m.bar <- otu.m
  otu.m.bar$exp_day <- factor(as.character(otu.m.bar$exp_day), levels = as.character(dates.ax), ordered = T)
  # otu.m.bar$Taxa <- lapply(X = otu.m.bar$Taxa, FUN = function(xx) strsplit(xx, ';', fixed = T)[[1]][length(strsplit(xx, ';', fixed = T))])
  ## Plot according to Lifestyle, sorted by abundance
  # png(paste0("TaxaSummary_L",bT[L],".png"),width = 12,height=8, units = "in", res = '300') # Make room for legend
  # ggplot(otu.m.bar, aes(x = samp_rel_hct, y = RelativeAbundance, fill = Taxa)) +
  #   geom_bar(stat ="identity", position="fill") + labs(x="date relative to transplant",y="root relative abundance") +
  #   guides(fill=guide_legend(ncol=1)) +
  #   scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
  #                              "green4",
  #                              "#6A3D9A", # purple
  #                              "#FF7F00", # orange
  #                              "black","gold1",
  #                              "skyblue2","#FB9A99", # lt pink
  #                              "palegreen2",
  #                              "#CAB2D6", # lt purple
  #                              "#FDBF6F", # lt orange
  #                              "gray70", "khaki2",
  #                              "maroon","orchid1","deeppink1","blue1","steelblue4",
  #                              "darkturquoise","green1","yellow4","yellow3",
  #                              "darkorange4","brown")) +
  #   theme_classic() + theme(axis.text = element_text(color="black"))
  # ggsave(filename = paste0("../../results/individualplots/patient97/Taxa_barplot_summary_L",bT[L],".png"), height = 7, width = 5.5, dpi = 300)
  
  otu.m.area <- otu.m
  otu.m.area$exp_day_n <- as.numeric(otu.m$exp_day)
  ggplot(otu.m.area, aes(x = exp_day_n, y = RelativeAbundance, fill = Taxa)) +
    geom_area(stat ="identity") + labs(x="experiment day",y="root relative abundance") +
    guides(fill=guide_legend(ncol=1)) + scale_x_continuous() +
    scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                               "green4",
                               "#6A3D9A", # purple
                               "#FF7F00", # orange
                               "black","gold1",
                               "skyblue2","#FB9A99", # lt pink
                               "palegreen2",
                               "#CAB2D6", # lt purple
                               "#FDBF6F", # lt orange
                               "gray70", "khaki2",
                               "maroon","orchid1","deeppink1","blue1","steelblue4",
                               "darkturquoise","green1","yellow4","yellow3",
                               "darkorange4","brown")) +
    theme_classic() + theme(axis.text = element_text(color="black"))
  ggsave(filename = paste0("../../results/mouse_plots_by_animal_",mouse.ids[m],"_L3.png"), height = 6, width = 7, dpi = 300)
}
