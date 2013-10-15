### R code from vignette source 'diveRsity.Rnw'

###################################################
### code chunk number 1: diveRsity.Rnw:66-67 (eval = FALSE)
###################################################
## divOnline()


###################################################
### code chunk number 2: diveRsity.Rnw:83-84 (eval = FALSE)
###################################################
## citation("diveRsity")


###################################################
### code chunk number 3: diveRsity.Rnw:124-125 (eval = FALSE)
###################################################
## install.packages("diveRsity")


###################################################
### code chunk number 4: diveRsity.Rnw:130-138 (eval = FALSE)
###################################################
## # install and load devtools
## install.packages("devtools")
## 
## library("devtools")
## 
## # download and install diveRsity
## 
## install_github(name = "diveRsity", subdir = "kkeenan02")


###################################################
### code chunk number 5: diveRsity.Rnw:157-158 (eval = FALSE)
###################################################
## install.packages("package_name")


###################################################
### code chunk number 6: diveRsity.Rnw:166-167 (eval = FALSE)
###################################################
## library("diveRsity")


###################################################
### code chunk number 7: diveRsity.Rnw:171-185 (eval = FALSE)
###################################################
## ?divPart
## ?inCalc
## ?readGenepop
## ?corPlot
## ?difPlot
## ?chiCalc
## ?divOnline
## ?divBasic
## ?fstOnly
## ?divRatio
## ?microPlexer
## ?arp2gen
## ?divMigrate
## ?haploDiv


###################################################
### code chunk number 8: diveRsity.Rnw:421-424 (eval = FALSE)
###################################################
## divPart(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
##         WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
##         bootstraps = 0, plot = FALSE, parallel = FALSE)


###################################################
### code chunk number 9: diveRsity.Rnw:480-486
###################################################
#data(Test_data,package='diveRsity')
#library(diveRsity)
#x<-capture.output(res<-divPart(Test_data,"outtt",3,T,T,T,T,3,F,T))

load("./div_res.RData")
res <- div_results


###################################################
### code chunk number 10: diveRsity.Rnw:496-498
###################################################
options(width=50)
res$standard[c(1:10,nrow(res$standard)),] 


###################################################
### code chunk number 11: diveRsity.Rnw:528-530
###################################################
options(width=50)
res$estimate[c(1:10,nrow(res$estimate)),]


###################################################
### code chunk number 12: diveRsity.Rnw:567-583
###################################################
noquote(names(res$pairwise)[1])
noquote(res$pairwise[[1]][1:4,1:4])
noquote(names(res$pairwise)[2])
noquote(res$pairwise[[2]][1:4,1:4])
noquote(names(res$pairwise)[3])
noquote(res$pairwise[[3]][1:4,1:4])
noquote(names(res$pairwise)[4])
noquote(res$pairwise[[4]][1:4,1:4])
noquote(names(res$pairwise)[5])
noquote(res$pairwise[[5]][1:4,1:4])
noquote(names(res$pairwise)[6])
noquote(res$pairwise[[6]][1:4,1:4])
noquote(names(res$pairwise)[7])
noquote(res$pairwise[[7]][1:4,1:4])
noquote(names(res$pairwise)[8])
noquote(res$pairwise[[8]][1:4,1:4])


###################################################
### code chunk number 13: diveRsity.Rnw:596-613
###################################################
noquote(names(res$bs_locus)[1])
res$bs_locus[[1]][c(1:3,nrow(res$bs_locus[[1]])),]
noquote(names(res$bs_locus)[2])
res$bs_locus[[2]][c(1:3,nrow(res$bs_locus[[2]])),]
noquote(names(res$bs_locus)[3])
res$bs_locus[[3]][c(1:3,nrow(res$bs_locus[[3]])),]
noquote(names(res$bs_locus)[4])
res$bs_locus[[4]][c(1:3,nrow(res$bs_locus[[4]])),]
noquote(names(res$bs_locus)[5])
res$bs_locus[[5]][c(1:3,nrow(res$bs_locus[[5]])),]
noquote(names(res$bs_locus)[6])
res$bs_locus[[6]][c(1:3,nrow(res$bs_locus[[6]])),]
noquote(names(res$bs_locus)[7])
res$bs_locus[[7]][c(1:3,nrow(res$bs_locus[[7]])),]
noquote(names(res$bs_locus)[8])
res$bs_locus[[8]][c(1:3,nrow(res$bs_locus[[8]])),]
#$


###################################################
### code chunk number 14: diveRsity.Rnw:628-645
###################################################
noquote(names(res$bs_pairwise)[1])
res$bs_pairwise[[1]][c(1:3,nrow(res$bs_pairwise[[1]])),]
noquote(names(res$bs_pairwise)[2])
res$bs_pairwise[[2]][c(1:3,nrow(res$bs_pairwise[[2]])),]
noquote(names(res$bs_pairwise)[3])
res$bs_pairwise[[3]][c(1:3,nrow(res$bs_pairwise[[3]])),]
noquote(names(res$bs_pairwise)[4])
res$bs_pairwise[[4]][c(1:3,nrow(res$bs_pairwise[[4]])),]
noquote(names(res$bs_pairwise)[5])
res$bs_pairwise[[5]][c(1:3,nrow(res$bs_pairwise[[5]])),]
noquote(names(res$bs_pairwise)[6])
res$bs_pairwise[[6]][c(1:3,nrow(res$bs_pairwise[[6]])),]
noquote(names(res$bs_pairwise)[7])
res$bs_pairwise[[7]][c(1:3,nrow(res$bs_pairwise[[7]])),]
noquote(names(res$bs_pairwise)[8])
res$bs_pairwise[[8]][c(1:3,nrow(res$bs_pairwise[[8]])),]
#$


###################################################
### code chunk number 15: diveRsity.Rnw:656-659 (eval = FALSE)
###################################################
## inCalc(infile, outfile = NULL, gp = 3, bs_locus = FALSE,
##        bs_pairwise = FALSE, bootstraps = 0, plot = FALSE,
##        parallel = FALSE)


###################################################
### code chunk number 16: diveRsity.Rnw:708-711
###################################################

load("./in_res.RData")
res_in <- in_results


###################################################
### code chunk number 17: diveRsity.Rnw:723-724
###################################################
noquote(res_in$Allele_In[1:10,c(1:5,ncol(res_in$Allele_In))])


###################################################
### code chunk number 18: diveRsity.Rnw:742-744
###################################################
res_in$l_bootstrap[1:10,]
#$


###################################################
### code chunk number 19: diveRsity.Rnw:763-773
###################################################
noquote(names(res_in$PW_bootstrap)[1])
res_in$PW_bootstrap[[1]][1:5,]
noquote(names(res_in$PW_bootstrap)[2])
res_in$PW_bootstrap[[2]][1:5,]
noquote(names(res_in$PW_bootstrap)[3])
res_in$PW_bootstrap[[3]][1:5,]
noquote(names(res_in$PW_bootstrap)[4])
res_in$PW_bootstrap[[4]][1:5,]
noquote(names(res_in$PW_bootstrap)[5])
res_in$PW_bootstrap[[5]][1:5,]


###################################################
### code chunk number 20: diveRsity.Rnw:780-781 (eval = FALSE)
###################################################
## readGenepop(infile = NULL, gp = 3, bootstrap = FALSE)


###################################################
### code chunk number 21: diveRsity.Rnw:841-842 (eval = FALSE)
###################################################
## corPlot(x,y)


###################################################
### code chunk number 22: diveRsity.Rnw:876-877 (eval = FALSE)
###################################################
## difPlot(x, outfile = NULL, interactive = FALSE)


###################################################
### code chunk number 23: diveRsity.Rnw:921-922 (eval = FALSE)
###################################################
## chiCalc(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL)


###################################################
### code chunk number 24: diveRsity.Rnw:949-950 (eval = FALSE)
###################################################
## divOnline()


###################################################
### code chunk number 25: diveRsity.Rnw:957-959 (eval = FALSE)
###################################################
## fstOnly(infile = NULL, outfile = NULL, gp = 3, bs_locus = FALSE,
##         bs_pairwise = FALSE, bootstraps = 0, parallel = FALSE)


###################################################
### code chunk number 26: diveRsity.Rnw:993-995 (eval = FALSE)
###################################################
## divRatio(infile = NULL, outfile = NULL, gp = 3, pop_stats = NULL,
##          refPos = NULL, bootstraps = 1000, parallel = FALSE)


###################################################
### code chunk number 27: diveRsity.Rnw:1048-1050 (eval = FALSE)
###################################################
## bigDivPart(infile = NULL, outfile = NULL, WC_Fst = FALSE,
##            format = NULL)


###################################################
### code chunk number 28: diveRsity.Rnw:1077-1078 (eval = FALSE)
###################################################
## microPlexer()


###################################################
### code chunk number 29: diveRsity.Rnw:1098-1099 (eval = FALSE)
###################################################
## arp2gen(infile)


###################################################
### code chunk number 30: diveRsity.Rnw:1107-1108 (eval = FALSE)
###################################################
## divMigrate(infile = NULL, stat = c("gst", "d_jost"))


###################################################
### code chunk number 31: diveRsity.Rnw:1116-1118 (eval = FALSE)
###################################################
## haploDiv(infile = NULL, outfile = NULL, pairwise = FALSE,
##          bootstraps = 0)


###################################################
### code chunk number 32: diveRsity.Rnw:1131-1132
###################################################
library("diveRsity")


###################################################
### code chunk number 33: diveRsity.Rnw:1142-1143 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 34: diveRsity.Rnw:1150-1151
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 35: diveRsity.Rnw:1160-1165 (eval = FALSE)
###################################################
## div_results <- divPart(infile = Test_data, outfile = "Test", 
##                          gp = 3, pairwise = TRUE, 
##                          WC_Fst = TRUE, bs_locus = TRUE, 
##                          bs_pairwise = TRUE, bootstraps = 100, 
##                          plot = FALSE, parallel = TRUE)


###################################################
### code chunk number 36: diveRsity.Rnw:1168-1170
###################################################

load("./div_res.RData")


###################################################
### code chunk number 37: diveRsity.Rnw:1180-1181
###################################################
names(div_results)


###################################################
### code chunk number 38: diveRsity.Rnw:1186-1187
###################################################
typeof(div_results$bs_locus)


###################################################
### code chunk number 39: diveRsity.Rnw:1190-1191
###################################################
#$


###################################################
### code chunk number 40: diveRsity.Rnw:1195-1196
###################################################
names(div_results$bs_locus)


###################################################
### code chunk number 41: diveRsity.Rnw:1198-1199
###################################################
#$


###################################################
### code chunk number 42: diveRsity.Rnw:1206-1207 (eval = FALSE)
###################################################
## mymatrix[5, 1]


###################################################
### code chunk number 43: diveRsity.Rnw:1213-1214
###################################################
div_results$bs_locus$Gst[1:10, ]


###################################################
### code chunk number 44: diveRsity.Rnw:1218-1219 (eval = FALSE)
###################################################
## div_results$bs_locus$Gst[ ,1]


###################################################
### code chunk number 45: diveRsity.Rnw:1231-1232 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 46: diveRsity.Rnw:1239-1240
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 47: diveRsity.Rnw:1249-1253 (eval = FALSE)
###################################################
## in_results <- inCalc (infile = Test_data, outfile = "Test", 
##                          gp = 3, bs_locus = TRUE, 
##                          bs_pairwise = TRUE, bootstraps = 100, 
##                          plot = FALSE, parallel = TRUE)


###################################################
### code chunk number 48: diveRsity.Rnw:1256-1258
###################################################
load("./in_res.RData")



###################################################
### code chunk number 49: diveRsity.Rnw:1271-1272
###################################################
names(in_results)


###################################################
### code chunk number 50: diveRsity.Rnw:1277-1278
###################################################
typeof(in_results$PW_bootstrap)


###################################################
### code chunk number 51: diveRsity.Rnw:1281-1282
###################################################
#$


###################################################
### code chunk number 52: diveRsity.Rnw:1287-1288
###################################################
names(in_results$PW_bootstrap)


###################################################
### code chunk number 53: diveRsity.Rnw:1290-1291
###################################################
#$


###################################################
### code chunk number 54: diveRsity.Rnw:1296-1297 (eval = FALSE)
###################################################
## mymatrix[5, 1]


###################################################
### code chunk number 55: diveRsity.Rnw:1303-1304
###################################################
in_results$PW_bootstrap[["pop1, vs. pop2,"]][1:3, ]


###################################################
### code chunk number 56: diveRsity.Rnw:1306-1307
###################################################
#$


###################################################
### code chunk number 57: diveRsity.Rnw:1311-1312 (eval = FALSE)
###################################################
## in_results$PW_bootstrap[["pop1, vs. pop2,"]][ ,1]


###################################################
### code chunk number 58: diveRsity.Rnw:1323-1324 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 59: diveRsity.Rnw:1331-1332
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 60: diveRsity.Rnw:1340-1342
###################################################
gp_res <- readGenepop(infile = Test_data, gp = 3,
                      bootstrap = FALSE)


###################################################
### code chunk number 61: diveRsity.Rnw:1350-1351
###################################################
names(gp_res)


###################################################
### code chunk number 62: diveRsity.Rnw:1366-1373
###################################################
locus18_pop1 <- c(gp_res$pop_alleles[[1]][[1]][,18], 
                  gp_res$pop_alleles[[2]][[1]][,18])
# sort alleles by size
allele_sort <- order(locus18_pop1, decreasing = FALSE)
#plot
plot(locus18_pop1[allele_sort], ylab = "allele size", col="blue",
     pch = 16)


###################################################
### code chunk number 63: diveRsity.Rnw:1389-1424
###################################################
# Define a results matrix with 37 columns (loci) and
# 1000 rows (bootstraps)to record allele number per locus

num_all <- matrix(rep(0,(37*10)), ncol = 37)

# Now using readGenepop we can fill the matrix
bs<-10
for(i in 1:bs){
    # first produce a bootstrap file
    
    x <- readGenepop(infile = Test_data, gp = 3,
                     bootstrap = TRUE)
                     
    # Now record the number of alleles at each locus
    
    num_all[i, ] <- x$nalleles                      
}

# Now we can use this data to calculate the mean
# number of alleles per locus as well at their
# 95% confidence intervals

mean_num <- colMeans(num_all)
lower<-vector()
upper<-vector()
for(i in 1:ncol(num_all)){
    lower[i] <- mean_num[i] - (1.96 * sd(num_all[,i]))  
    upper[i] <- mean_num[i] + (1.96 * sd(num_all[,i]))
}

# Now we can create a data frame of these results

bs_res <- data.frame(mean_num, lower, upper)

bs_res[1:10,]


###################################################
### code chunk number 64: diveRsity.Rnw:1449-1507 (eval = FALSE)
###################################################
##  # load the diveRsity package 
## library("diveRsity")
## 
## # We can specify the names of our simulation folders in two ways 
## 
## # manually 
## fold_names <- paste("sim", 1:10, sep = "")
## 
## or
## 
## # automatically (when there is only a single level below the
## # top directory)
## fold_names <- list.dirs(full.names = TRUE, recursive = FALSE) 
## 
## # Now we can determine the names of all genepop files in each folder
## 
## file_names <- lapply(fold_names, function(x){
##   files <- dir(path = paste(x, "/", sep = ""), pattern = "*.gen",
##                full.names = TRUE)
##   return(files) })
## 
## # file_names will be a list of length 10. Each element will contain
## # the names of all .gen files within the respective simulation folder
## 
## # Before we are ready to run the main analyses, we should set up 
## # the parallel environment
## 
## # load the doParallel package
## library("doParallel")
## # set up a cluster of 10 CPUs (one for each batch of files)
## cl <- makeCluster(10)
## # Export the 'divPart' function to the cluster cores 
## clusterExport(cl, "divPart", envir = environment())
## 
## # Now we can run the main analyses
## 
## results <- parLapply(cl, file_names, function(x){
##   sim_res <- sapply(x, function(y){
##     out <- divPart(infile = y, gp = 3, WC_Fst = TRUE)
##     return(out$estimate[nrow(out$estimate), 4:7])
##   })
##   return(t(sim_res))  # transpose sim_res
## })
## 
## # This will generate a list (of length 10), with each element 
## # containing a matrix of 1000 rows (1 per file) and 4 columns 
## # (1 for each diversity statistic) 
## 
## # Example of output for simulation 1
## # G_st_est        G_hed_st_est    D_Jost_est      Fst_WC
## # 0.3905          0.8938          0.8256          0.4010
## # 0.5519          0.8719          0.6986          0.6031
## # 0.5924          0.8880          0.7092          0.6096
## # ...             ...             ...             ...
## # ...             ...             ...             ...
## 
## # these results could then be piped to further analyses or 
## # visualisation tools


###################################################
### code chunk number 65: diveRsity.Rnw:1522-1523
###################################################
print(sessionInfo(), locale = FALSE)


