library(ggplot2)


unifrac <- read.delim('../../data/prok/beta_diversity/unweighted_unifrac_distance_matrix.txt', header = T, row.names = 1, check.names = F)
unifrac <- unifrac[sort(rownames(unifrac)),sort(rownames(unifrac))]
map <- read.delim('../../data/metadata_vlavage.txt', header=T, check.names = F)
map$rel_exp_day <- map$exp_day - 2
rownames(map) <- map$sampleid
map <- map[rownames(unifrac),]
mouse.ids <- as.character(unique(map$mouse))
# mx <- 1
self.distances <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(self.distance) <- c("dist_to_0", "mouse", "exp_day", "trx")
for (mx in 1:length(mouse.ids)) {
  mouse.x <-sort(as.character(map[map$mouse==mouse.ids[mx],]$sampleid))
  self.dist.mx <- data.frame(matrix(nrow = length(mouse.x), ncol = 4))
  colnames(self.dist.mx) <- c("dist_to_0", "mouse", "exp_day", "trx")
  self.dist.mx$mouse <- mouse.ids[mx]
  self.dist.mx$exp_day <- as.numeric(map[map$mouse==mouse.ids[mx],]$rel_exp_day)
  self.dist.mx$trx <- as.character(map[map$mouse==mouse.ids[mx],]$treatment)
  self.dist.mx$dist_to_0 <- as.numeric(unifrac[mouse.x, mouse.x[1]])
  self.distances <- rbind(self.distances, self.dist.mx)
}

ggplot(self.distances, aes(x=rel_exp_day, y=dist_to_0, group=mouse, color=trx)) +
  geom_line() +
  theme_classic() + labs(x="experiment day", y="distance to pre-sensitization") +
  scale_color_manual(values=c('red','blue'))
ggsave('../../results/distance_self_unweighted_unifrac.png', height = 4, width = 5, dpi=300)


ggplot(self.distances, aes(x=rel_exp_day, y=dist_to_0, group=trx)) +
  geom_smooth(aes(color=trx)) +
  theme_classic() + labs(x="experiment day", y="distance to pre-sensitization") +
  scale_color_manual(values=c('red','blue'))
ggsave('../../results/distance_self_splines_unweighted_unifrac.png', height = 4, width = 5, dpi=300)
  


### Now distance relative to prev. day

unifrac <- read.delim('../../data/prok/beta_diversity/unweighted_unifrac_distance_matrix.txt', header = T, row.names = 1, check.names = F)
unifrac <- unifrac[sort(rownames(unifrac)),sort(rownames(unifrac))]
map <- read.delim('../../data/metadata_vlavage.txt', header=T, check.names = F)
map$rel_exp_day <- map$exp_day - 2
rownames(map) <- map$sampleid
map <- map[rownames(unifrac),]
mouse.ids <- as.character(unique(map$mouse))

# mx <- 1
self.distances <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(self.distances) <- c("dist", "mouse", "rel_exp_day", "trx")
for (mx in 1:length(mouse.ids)) {
  mouse.x <-sort(as.character(map[map$mouse==mouse.ids[mx],]$sampleid))
  if (length(mouse.x) < 4) next
  self.dist.mx <- data.frame(matrix(nrow = (length(mouse.x)-1), ncol = 4))
  colnames(self.dist.mx) <- c("dist", "mouse", "rel_exp_day", "trx")
  self.dist.mx$mouse <- (mouse.ids[mx])[1]
  self.dist.mx$rel_exp_day <- as.numeric(map[map$mouse==mouse.ids[mx],]$rel_exp_day)[2:length(mouse.x)]
  self.dist.mx$trx <- as.character(map[map$mouse==mouse.ids[mx],]$treatment)[1]
  mx.distances <- diag(as.matrix(unifrac[mouse.x[1:(length(mouse.x)-1)], mouse.x[2:length(mouse.x)]]))
  dist.1 <- mx.distances[1]
  mx.distances <- unlist(lapply(X = mx.distances, FUN = function(xx) xx - dist.1))
  self.dist.mx$dist <- as.numeric(mx.distances)
  self.distances <- rbind(self.distances, self.dist.mx)
}

ggplot(self.distances, aes(x=rel_exp_day, y=dist, group=mouse, color=trx)) +
  geom_line() +
  theme_classic() + labs(x="experiment day relative to first sensitization", y="relative distance to previous sample") +
  scale_color_manual(name='Treatment', values=c('red','blue'))
ggsave('../../results/distance_self_unweighted_unifrac.png', height = 4, width = 5, dpi=300)


ggplot(self.distances, aes(x=rel_exp_day, y=dist, group=trx)) +
  geom_smooth(aes(color=trx)) +
  theme_classic() + labs(x="experiment day relative to first sensitization", y="relative distance to previous sample") +
  scale_color_manual(name='Treatment', values=c('red','blue'))
ggsave('../../results/distance_self_splines_unweighted_unifrac.png', height = 4, width = 5, dpi=300)
