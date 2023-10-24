library(phytools)
library(rphylopic)
setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

addSil <- function (img, x = NULL, y = NULL, ysize = 1, xsize = 1, alpha = 0.7, color = NULL) {                                                                         
	img <- rphylopic:::recolor_phylopic(img, alpha, color)
	graphics::rasterImage(img, x - xsize/2, y - ysize/2, x + xsize/2, y + ysize/2, interpolate = TRUE)                         
} 

addLoadings <- function(Axis, Scores, Loadings, Side=3, Thresh = 0.4, Buffer=0.2, Cex=0.5)	{
	# par("usr") = x1, x2, y1, y2
	Keeps <- which(abs(Loadings[,Axis]) >= Thresh)
	ToMark <- Loadings[Keeps,Axis]
	ToMark <- ToMark * max(abs(Scores[,Axis]))
	Names <- gsub("\\_", " ", rownames(Loadings)[Keeps])
	if (Side == 1)	{
		Ax <- 3
		silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[Ax]+Buffer, Names[x], adj = 0, srt = 35, xpd=NA, cex=Cex))
	}
	if (Side == 2)	{
		Ax <- 1
		silent <- sapply(1:length(ToMark), function(x) text(par("usr")[Ax]+Buffer, ToMark[x], Names[x], adj = 0, srt = 35, xpd=NA, cex=Cex))
	}
	if (Side == 3)	{
		Ax <- 4
		silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[Ax]+Buffer, Names[x], adj = 0, srt = 35, xpd=NA, cex=Cex))
	}
	if (Side == 4)	{
		Ax <- 2
		silent <- sapply(1:length(ToMark), function(x) text(par("usr")[Ax]+Buffer, ToMark[x], Names[x], adj = 0, srt = 35, xpd=NA, cex=Cex))
	}

}

Genera <- sapply(rownames(Scores$S), function(x) unlist(strsplit(x, " "))[1])


setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
pdf("pPCA_NoLoadings.pdf", height=5, width=10)
#par(mfrow=c(1, 2), las=1, mgp=c(1.5, 0.5, 0), tck=-0.01, mar=c(3,3,4,4), bty="l")
par(mfrow=c(1, 2), las=1, mgp=c(1.5, 0.5, 0), tck=-0.01, mar=c(4,4,1,1), bty="l", xpd=NA)
Xchoice <- 1
Ychoice <- 2
X <- Scores$S[,Xchoice]
Y <- Scores$S[,Ychoice]
PCH <- c(16, 5)
plot(X, Y, pch=PCH[ModelDat$flight + 1], xlab=paste("pPC", Xchoice, sep=""), ylab=paste("pPC ", Ychoice, sep=""))
#addLoadings(Xchoice, Scores$S, Scores$L)
#addLoadings(Ychoice, Scores$S, Scores$L, Side = 4)
legend("topleft", legend=c("flightless", "volant"), pch=PCH, bty="n")

### Add silhouettes
Genus <- "Struthio"
Sil <- image_data("df5979fc-7cdd-41c7-ac47-7107ba3f7e1a", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=9, xsize=15, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Aptenodytes"
Sil <- image_data("21c50828-58d8-42df-9219-2c1b0fb57c99", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=7, xsize=15, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Rhea"
Sil <- image_data("83097107-2510-44e4-b899-273362c7ca62", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=6, xsize=18, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Dryocopus"
Sil <- image_data("ddd5783c-ded5-48f2-a07d-cc37c83b227b", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=10, xsize=15, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Columba"
Sil <- image_data("93124ce1-73c9-4278-bfd4-a29886bc7ca7", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=6, xsize=18, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Cathartes"
Sil <- image_data("72a6434a-3ed4-43f0-b6ef-b3146c2f07ae", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=7.5, xsize=18+9, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Casuarius"
Sil <- image_data("6fe53db2-831d-47c4-84db-34b7b429cc40", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=6, xsize=18, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Rollandia"
Sil <- image_data("deba1d91-daa8-40a6-8d48-7a9f295bc662", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus)[1],Xchoice], y=Scores$S[which(Genera==Genus)[1],Ychoice], ysize=6, xsize=18, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus)[1],Xchoice], y=Scores$S[which(Genera==Genus)[1],Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Xchoice <- 3
Ychoice <- 4
X <- Scores$S[,Xchoice]
Y <- Scores$S[,Ychoice]
plot(X, Y, pch=PCH[ModelDat$flight + 1], xlab=paste("pPC", Xchoice, sep=""), ylab=paste("pPC ", Ychoice, sep=""))
#addLoadings(Xchoice, Scores$S, Scores$L)
#addLoadings(Ychoice, Scores$S, Scores$L, Side = 4)

### Add silhouettes
Genus <- "Struthio"
Sil <- image_data("df5979fc-7cdd-41c7-ac47-7107ba3f7e1a", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=6, xsize=4, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Aptenodytes"
Sil <- image_data("21c50828-58d8-42df-9219-2c1b0fb57c99", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=5, xsize=4*(5/6), alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Dryocopus"
Sil <- image_data("ddd5783c-ded5-48f2-a07d-cc37c83b227b", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=5, xsize=3, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Casuarius"
Sil <- image_data("6fe53db2-831d-47c4-84db-34b7b429cc40", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=6, xsize=4, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Phoenicopterus"
Sil <- image_data("a1244226-f2c2-41dc-b113-f1c6545958ce", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=7, xsize=3, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Caprimulgus"
Sil <- image_data("19899811-afff-4ba8-8084-cc7fe762da7e", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=3, xsize=6, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])

Genus <- "Mesitornis"
Sil <- image_data("3baf6f47-7496-45cc-b21f-3f0bee99d65b", size = 512)[[1]]
addSil(Sil, x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], ysize=4, xsize=5, alpha=0.7, color="red")
points(x=Scores$S[which(Genera==Genus),Xchoice], y=Scores$S[which(Genera==Genus),Ychoice], pch=PCH[ModelDat$flight[which(Genera==Genus)] + 1])


dev.off()