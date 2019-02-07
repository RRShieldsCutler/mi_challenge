map <- read.delim('~/Documents/Research/chatterjea/data/metadata_vlavage.txt',
                  header = T, sep = '\t')

otdf <- read.delim('~/Documents/Research/chatterjea/data/vlav6b_ggotu.txt',
                   row.names = 1, header = T, sep = '\t')
sum(otdf)
dim(otdf)

# otdf <- otdf[-grep('Chloroplast', rownames(otdf)),]
otdf <- otdf[rowSums(otdf) > 2, ]
otdf <- otdf[rowMeans(otdf > 0) >= 0.05, ]

depths <- colSums(otdf)
sort(depths)[1:20]

otdf <- otdf[, colSums(otdf)>2280]
otdf <- otdf[, !names(otdf) %in% c("PBS.control")]

map.g <- map[map$sampleid %in% colnames(otdf), ]

otdf.w <- tibble::rownames_to_column(otdf, var='OTU ID')

write.table(otdf.w, '~/Documents/Research/chatterjea/data/vlav6b_filt_ggotu.txt', row.names = F, sep='\t',quote = F)
write.table(map.g, '~/Documents/Research/chatterjea/data/metadata_gg.txt', row.names = F, sep='\t',quote = F)


#######for prok

otdf <- read.delim('~/Documents/Research/chatterjea/data/vlav6b_prokotu.txt',
                   row.names = 1, header = T, sep = '\t')
sum(otdf)
dim(otdf)
# otdf <- otdf[-grep('Chloroplast', rownames(otdf)),]
otdf <- otdf[rowSums(otdf) >= 2, ]
otdf <- otdf[rowMeans(otdf > 0) >= 0.05, ]

depths <- colSums(otdf)
sort(depths)[1:41]

otdf <- otdf[, colSums(otdf)>195]
otdf <- otdf[, !names(otdf) %in% c("BLANK.B06")]
otdf <- otdf[, !names(otdf) %in% c("BLANK.C06")]
otdf <- otdf[, !names(otdf) %in% c("PBS.control")]

map.p <- map[map$sampleid %in% colnames(otdf), ]

otdf.w <- tibble::rownames_to_column(otdf, var='OTU ID')

write.table(otdf.w, '~/Documents/Research/chatterjea/data/vlav6b_filt_prokotu.txt', row.names = F, sep='\t',quote = F)
write.table(map.p, '~/Documents/Research/chatterjea/data/metadata_prok.txt', row.names = F, sep='\t',quote = F)
