library(HDInterval)

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
barb_ang <- read.csv("barb_ang_MM.csv")
barb_ang <- barb_ang[floor(0.5*nrow(barb_ang)):nrow(barb_ang),]
barb_dens <- read.csv("barb_dens_MM.csv")
barb_dens <- barb_dens[floor(0.5*nrow(barb_dens)):nrow(barb_dens),]
barb_len <- read.csv("barb_len_MM.csv")
barb_len <- barb_len[floor(0.5*nrow(barb_len)):nrow(barb_len),]
barbule_dens <- read.csv("barbule_dens_MM.csv")
barbule_dens <- barbule_dens[floor(0.5*nrow(barbule_dens)):nrow(barbule_dens),]
wing <- read.csv("wing_MM.csv")
wing <- wing[floor(0.5*nrow(wing)):nrow(wing),]

options(scipen = -3)

Obj <- barb_ang
Table <- rbind(apply(Obj, 2, function(x) signif(mean(x), digits = 2)), apply(Obj, 2, function(x) signif(hdi(x), digits=2)))
Table <- Table[,c("b0","b1","b2","b3","s1","s2","lam")]
# intercept comparison
Table <- cbind(Table, c(signif(mean(Obj$b0 - Obj$b2), digits = 2), signif(hdi(Obj$b0 - Obj$b2), digits = 2))) 
# slope comparison
Table <- cbind(Table, c(signif(mean(Obj$b1 - Obj$b3), digits = 2), signif(hdi(Obj$b1 - Obj$b3), digits = 2))) 
# rate comparison
Table <- cbind(Table, c(signif(mean(log(Obj$s1) - log(Obj$s2)), digits = 2), signif(hdi(log(Obj$s1) - log(Obj$s2)), digits =2))) 
sum_ang <- paste(Table[1,], " (", Table[2,], ", ", Table[3,], ")", sep = "")


Obj <- barb_dens
Table <- rbind(apply(Obj, 2, function(x) signif(mean(x), digits = 2)), apply(Obj, 2, function(x) signif(hdi(x), digits=2)))
Table <- Table[,c("b0","b1","b2","b3","s1","s2","lam")]
# intercept comparison
Table <- cbind(Table, c(signif(mean(Obj$b0 - Obj$b2), digits = 2), signif(hdi(Obj$b0 - Obj$b2), digits = 2))) 
# slope comparison
Table <- cbind(Table, c(signif(mean(Obj$b1 - Obj$b3), digits = 2), signif(hdi(Obj$b1 - Obj$b3), digits = 2))) 
# rate comparison
Table <- cbind(Table, c(signif(mean(log(Obj$s1) - log(Obj$s2)), digits = 2), signif(hdi(log(Obj$s1) - log(Obj$s2)), digits =2))) 
sum_dens <- paste(Table[1,], " (", Table[2,], ", ", Table[3,], ")", sep = "")

Obj <- barb_len
Table <- rbind(apply(Obj, 2, function(x) signif(mean(x), digits = 2)), apply(Obj, 2, function(x) signif(hdi(x), digits=2)))
Table <- Table[,c("b0","b1","b2","b3","s1","s2","lam")]
# intercept comparison
Table <- cbind(Table, c(signif(mean(Obj$b0 - Obj$b2), digits = 2), signif(hdi(Obj$b0 - Obj$b2), digits = 2))) 
# slope comparison
Table <- cbind(Table, c(signif(mean(Obj$b1 - Obj$b3), digits = 2), signif(hdi(Obj$b1 - Obj$b3), digits = 2))) 
# rate comparison
Table <- cbind(Table, c(signif(mean(log(Obj$s1) - log(Obj$s2)), digits = 2), signif(hdi(log(Obj$s1) - log(Obj$s2)), digits =2))) 
sum_len <- paste(Table[1,], " (", Table[2,], ", ", Table[3,], ")", sep = "")

Obj <- wing
Table <- rbind(apply(Obj, 2, function(x) signif(mean(x), digits = 2)), apply(Obj, 2, function(x) signif(hdi(x), digits=2)))
Table <- Table[,c("b0","b1","b2","b3","s1","s2","lam")]
# intercept comparison
Table <- cbind(Table, c(signif(mean(Obj$b0 - Obj$b2), digits = 2), signif(hdi(Obj$b0 - Obj$b2), digits = 2))) 
# slope comparison
Table <- cbind(Table, c(signif(mean(Obj$b1 - Obj$b3), digits = 2), signif(hdi(Obj$b1 - Obj$b3), digits = 2))) 
# rate comparison
Table <- cbind(Table, c(signif(mean(log(Obj$s1) - log(Obj$s2)), digits = 2), signif(hdi(log(Obj$s1) - log(Obj$s2)), digits =2))) 
sum_wing <- paste(Table[1,], " (", Table[2,], ", ", Table[3,], ")", sep = "")


sum_table <- rbind(sum_dens, sum_len, sum_ang, sum_wing)
sum_table <- cbind(rownames(sum_table), sum_table)
colnames(sum_table) <- c("traits", "b0","b1","b2","b3","s1","s2","lam", "b0_to_b2", "b1_to_b3", "s1_to_s2")

write.table(sum_table, "mm_summary_table.csv", sep=";", quote=F, row.names=F)


