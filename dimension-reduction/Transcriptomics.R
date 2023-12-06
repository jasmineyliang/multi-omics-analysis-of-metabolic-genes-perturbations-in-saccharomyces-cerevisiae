library("Rtsne")
library(data.table)
library("factoextra")
library("FactoMineR")

df <- read.table(file="G:\\My Drive\\pythonProject\\COMS570FinalProject\\TranscriptomicsData.txt", header=FALSE, sep="\t")
cond = read.table(file="G:\\My Drive\\pythonProject\\COMS570FinalProject\\Conds.txt", header=FALSE, sep="\t")

#dataP <- transpose(df)
dataP <- t(df)
dataPP <- t(df)
install.packages("janitor")

library(janitor)
dataPP<- remove_constant(dataPP)
dataPP<- remove_constant(dataPP, na.rm= TRUE)
dataPP.pr <- prcomp(dataPP, center = TRUE, scale = TRUE)
summary(dataPP.pr)

res.pca <- PCA(dataPP, graph = TRUE)
eig.val <- get_eigenvalue(res.pca)

# screeplot
fviz_eig(res.pca, addlabels = TRUE,ylim = c(0, 80),ncp = 10)
cumpro <- cumsum(dataPP.pr$sdev^2 / sum(dataPP.pr$sdev^2))

plot(cumpro[0:30], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 5, col="blue", lty=5)
abline(v = 10, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

# ELLIPSES ----------------------------------------------------------------
# plots PCs 1 and 2 labelled
# plot with ellipses

plot(dataPP.pr$x[,1],dataPP.pr$x[,2], xlab="PC1 (58.8%)", ylab = "PC2 (10.1%)", main = "PC1 / PC2 - plot")
fviz_pca_ind(dataPP.pr, geom.ind = "point", pointshape = 21,
             pointsize = 5,
             fill.ind = t(cond),
             col.ind = "black",
             palette =c('#e6194b', '#3cb44b', '#ffe119', '#000000', '#f58231',
                             '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                             '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                             '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
                             '#ffffff', '#000000'),
         #   palette = 'GSEA',
         #   addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Condictions") +
  #  ggtitle("2D PCA-plot from 30 feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))

------------------------------------------------------------------------
library(Rtsne)
tsne_out1 <- Rtsne(dataP,dims = 2, initial_dims = 30,perplexity = 3, 
                  theta = 0.5, 
                  check_duplicates = TRUE,
                  pca = TRUE, partial_pca = FALSE, max_iter = 5000,
                  verbose = getOption("verbose", FALSE), is_distance = FALSE,
                  Y_init = NULL, pca_center = TRUE, pca_scale = FALSE,
                  normalize = TRUE, stop_lying_iter = 250, 
                  mom_switch_iter = 250,
                  momentum = 0.5, final_momentum = 0.8, eta = 200,
                  exaggeration_factor = 12)



colors <- c('#e6194b', '#3cb44b', '#ffe119', '#000000', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
        '#ffffff', '#000000')
d_tsne_1 = as.data.frame(tsne_out1$Y)  
plot(d_tsne_1,
  #   ylim = c(-200,200),
  #   xlim = c(-200,200),
     
     col=colors[as.factor(cond)],
     cex = 2,
     pch = 16,
    # asp=1,
    ) # Plot the result

legend("bottomright", 
  #     inset = c(0, 2),
       cex = 1.2,
       legend = levels(as.factor(cond)), col =colors, 
       pch = 16, 
       bty = "n"
       )
#legend("topright", legend = as.factor(cond), col = 1:16, pch = 19, bty = "n")

--------------------------------------------------------------------------------------------------------------
  library(M3C)
  umap(t(dataP),
     labels=as.factor(cond),
     controlscale=TRUE,
     scale=3,
     dotsize = 7,
     textlabelsize= 2,
     #colvec=c('skyblue')
     colvec =c('#e6194b', '#3cb44b', '#ffe119', '#000000', '#f58231',
               '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
               '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
               '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
               '#ffffff', '#000000'),
     #text = c('1')
)
