## ------------------------------------------------------------------------
library(gapmap)
data("sample_tcga")

## ------------------------------------------------------------------------
#transpose
dataTable <- t(sample_tcga)
#calculate the correlation based distance
row_dist <- as.dist(1-cor(t(dataTable), method = "pearson"))
col_dist <- as.dist(1-cor(dataTable, method = "pearson"))
#hierarchical clustering
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")
col_d <- as.dendrogram(col_hc)
row_d <- as.dendrogram(row_hc)

## ----, fig.height=13, fig.width=20---------------------------------------
gapmap(m = as.matrix(dataTable), d_row = rev(row_d), d_col = col_d, mode = "quantitative", mapping="exponential", 
       ratio = 0.3, verbose=FALSE, scale = 0.5, label_size=2, v_ratio= c(0.1,0.8,0.1), h_ratio=c(0.1,0.8,0.1))

## ----, fig.height=13, fig.width=20---------------------------------------
library(dendsort)
gapmap(m = as.matrix(dataTable), d_row = rev(dendsort(row_d, type = "average")), d_col = dendsort(col_d, type = "average"),  
       mode = "quantitative", mapping="exponential", ratio = 0.3, verbose=FALSE, scale = 0.5, v_ratio= c(0.1,0.8,0.1), 
       h_ratio=c(0.1,0.8,0.1), label_size=2, show_legend=TRUE)

## ----, eval=FALSE--------------------------------------------------------
#  png("myplot.png", width= 2000, height= 1300, units="px")
#  gapmap(m = as.matrix(dataTable), d_row = rev(dendsort(row_d, type = "average")), d_col = dendsort(col_d, type = "average"),
#         mode = "quantitative", mapping="exponential", ratio = 0.3, verbose=FALSE, scale = 0.5, v_ratio= c(0.1,0.8,0.1),
#         h_ratio=c(0.1,0.8,0.1), label_size=2, show_legend=TRUE)
#  dev.off()

## ----, eval=FALSE--------------------------------------------------------
#  pdf("myplot.pdf", width= 20, height= 13)
#  gapmap(m = as.matrix(dataTable), d_row = rev(dendsort(row_d, type = "average")), d_col = dendsort(col_d, type = "average"),
#         mode = "quantitative", mapping="exponential", ratio = 0.3, verbose=FALSE, scale = 0.5, v_ratio= c(0.1,0.8,0.1),
#         h_ratio=c(0.1,0.8,0.1), label_size=2, show_legend=TRUE)
#  dev.off()

