library(mixOmics) # import the mixOmics library

tdata <- read.delim('filtered_transcriptomic.txt', header = TRUE, row.names = 1, sep = "\t")
pdata <- read.delim('proteome.txt', header = TRUE, row.names = 2, sep = "\t")
metadata <- read.delim('metadata.txt', header = TRUE, sep = "\t")

tdata <- tdata[,c(5:52)]
pdata <- pdata[,c(7:54)]
row.names(pdata) <- paste(row.names(pdata),".p",sep="")

tdata <- t(tdata)
pdata <- t(pdata)

pca.t <- pca(tdata, ncomp = 10, center = TRUE, scale = TRUE)
pca.p <- pca(pdata, ncomp = 10, center = TRUE, scale = TRUE)

plot(pca.t)
plot(pca.p)

plotIndiv(pca.t, comp = c(1, 2), 
          group = metadata[, 2], 
          ind.names = metadata[, 2], 
          legend = TRUE, title = 'Transcriptome, PCA comp 1 - 2')
plotIndiv(pca.p, comp = c(1, 2), 
          group = metadata[, 2], 
          ind.names = metadata[, 2], 
          legend = TRUE, title = 'Proteome, PCA comp 1 - 2')

spls.res <- spls(X = tdata, Y = pdata, ncomp = 5, mode = 'regression')

perf.spls.res <- perf(spls.res, validation = 'Mfold',
                        folds = 10, nrepeat = 5)

plot(perf.spls.res, criterion = 'Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 

tune.spls.tp <- tune.spls(tdata, pdata, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor') 
plot(tune.spls.tp)         # use the correlation measure for tuning

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.tp$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.tp$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

final.spls.tp <- spls(tdata, pdata, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression")

plotIndiv(final.spls.tp, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = metadata[, 2], # colour by time group
          col.per.group = color.jet(16), 
          legend = TRUE, legend.title = 'Group')

plotIndiv(final.spls.tp, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = metadata[, 2], # colour by time group
          col.per.group = color.jet(16), 
          legend = TRUE, legend.title = 'Group')

plotIndiv(final.spls.tp, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = metadata[, 2], # colour by time group
          col.per.group = color.jet(16), 
          legend = TRUE, legend.title = 'Group')

plotArrow(final.spls.tp, ind.names = FALSE,
          group = metadata[, 2],
          col.per.group = color.jet(16), 
          legend.title = 'Transcriptome.Proteome')

plotVar(final.spls.tp, cex = c(3,4), var.names = c(FALSE, TRUE))

color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
network.res <- network(final.spls.tp, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Transcriptome-Proteome Network Plot')

network.resM <- network.res$M


twoOmicsGraph <- matrix(ncol = 5, nrow = length(network.resM[,1])*length(network.resM[1,]))
count <- 1
for (i in 1:length(network.resM[,1])) {
  for (j in 1:length(network.resM[1,])){
    if (network.resM[i,j]>0){
      twoOmicsGraph[count,1] <- paste(row.names(network.resM)[i],".t",sep="")
      twoOmicsGraph[count,2] <- 0
      twoOmicsGraph[count,3] <- colnames(network.resM)[j]
      twoOmicsGraph[count,4] <- 0
      twoOmicsGraph[count,5] <- network.resM[i,j]
      count <- count+1
    }
  }
}

twoOmicsGraph <- na.omit(twoOmicsGraph)
twoOmicsGraph[,2] <- gsub(".t","",as.character(twoOmicsGraph[,1]))
twoOmicsGraph[,4] <- gsub(".p","",as.character(twoOmicsGraph[,3]))
colnames(twoOmicsGraph) <- c("Transcriptome.t","Transcriptome","Proteome.p","Proteome","Score")
write.table(twoOmicsGraph, 'twoOmicsGraph_v2.txt', quote=FALSE, row.names=FALSE, col.names = TRUE, sep="\t")








