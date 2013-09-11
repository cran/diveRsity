
## @knitr eval=FALSE
## divOnline()


## @knitr eval=FALSE
## citation("diveRsity")


## @knitr eval=FALSE
## install.packages("diveRsity")


## @knitr eval=FALSE
## install.packages("package_name")


## @knitr eval=FALSE
## library("diveRsity")


## @knitr eval = FALSE
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


## @knitr eval=FALSE
## divPart(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
##         WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE,
##         bootstraps = 0, plot = FALSE, parallel = FALSE)


## @knitr echo=FALSE, quiet=TRUE, results='hide'
#data(Test_data,package='diveRsity')
#library(diveRsity)
#x<-capture.output(res<-divPart(Test_data,"outtt",3,T,T,T,T,3,F,T))
load("./div_res.RData")
res <- div_results


## @knitr echo=FALSE
options(width=50)
res$standard[c(1:10,nrow(res$standard)),] 


## @knitr echo=FALSE
options(width=50)
res$estimate[c(1:10,nrow(res$estimate)),]


## @knitr echo=FALSE
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


## @knitr echo=FALSE
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


## @knitr echo=FALSE
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


## @knitr eval=FALSE
## inCalc(infile, outfile = NULL, gp = 3, bs_locus = FALSE,
##         bs_pairwise = FALSE, bootstraps = 0, plot = FALSE
##         parallel = FALSE)


## @knitr echo=FALSE, results='hide'
load("./in_res.RData")
res_in <- in_results


## @knitr echo=FALSE
noquote(res_in$Allele_In[1:10,c(1:5,ncol(res_in$Allele_In))])


## @knitr echo=FALSE
res_in$l_bootstrap[1:10,]
#$


## @knitr echo=FALSE
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


## @knitr eval=FALSE
## readGenepop(infile = NULL, gp = 3, bootstrap = FALSE)


## @knitr eval=FALSE
## corPlot(x,y)


## @knitr eval=FALSE
## difPlot(x, outfile = NULL, interactive = FALSE)


## @knitr eval=FALSE
## chiCalc(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL)


## @knitr eval=FALSE
## divOnline()


## @knitr eval=FALSE
## fstOnly(infile = NULL, outfile = NULL, gp = 3, bs_locus = FALSE,
##         bs_pairwise = FALSE, bootstraps = 0, parallel = FALSE)


## @knitr eval=FALSE
## divRatio(infile = NULL, outfile = NULL, gp = 3, pop_stats = NULL,
##          refPos = NULL, bootstraps = 1000, parallel = FALSE)


## @knitr eval=FALSE
## bigDivPart(infile = NULL, outfile = NULL, WC_Fst = FALSE,
##            format = NULL)


## @knitr eval=FALSE
## microPlexer()


## @knitr echo=FALSE
library("diveRsity")


## @knitr eval=FALSE
## setwd("mypath")


## @knitr 
data(Test_data, package = "diveRsity")


## @knitr eval=FALSE
## div_results <- divPart(infile = Test_data, outfile = "Test",
##                          gp = 3, pairwise = TRUE,
##                          WC_Fst = TRUE, bs_locus = TRUE,
##                          bs_pairwise = TRUE, bootstraps = 100,
##                          plot = FALSE, parallel = TRUE)


## @knitr echo=FALSE
load("./div_res.RData")


## @knitr 
names(div_results)


## @knitr 
typeof(div_results$bs_locus)


## @knitr echo=FALSE
#$


## @knitr 
names(div_results$bs_locus)


## @knitr echo=FALSE
#$


## @knitr eval=FALSE
## mymatrix[5, 1]


## @knitr 
div_results$bs_locus$Gst[1:10, ]


## @knitr eval=FALSE
## div_results$bs_locus$Gst[ ,1]


## @knitr eval=FALSE
## setwd("mypath")


## @knitr 
data(Test_data, package = "diveRsity")


## @knitr eval=FALSE
## in_results <- inCalc (infile = Test_data, outfile = "Test",
##                          gp = 3, bs_locus = TRUE,
##                          bs_pairwise = TRUE, bootstraps = 100,
##                          plot = FALSE, parallel = TRUE)


## @knitr echo=FALSE
load("./in_res.RData")


## @knitr 
names(in_results)


## @knitr 
typeof(in_results$PW_bootstrap)


## @knitr echo=FALSE
#$


## @knitr 
names(in_results$PW_bootstrap)


## @knitr echo=FALSE
#$


## @knitr eval=FALSE
## mymatrix[5, 1]


## @knitr 
in_results$PW_bootstrap[["pop1, vs. pop2,"]][1:3, ]


## @knitr echo=FALSE
#$


## @knitr eval=FALSE
## in_results$PW_bootstrap[["pop1, vs. pop2,"]][ ,1]


## @knitr eval=FALSE
## setwd("mypath")


## @knitr 
data(Test_data, package = "diveRsity")


## @knitr 
gp_res <- readGenepop(infile = Test_data, gp = 3,
                      bootstrap = FALSE)


## @knitr 
names(gp_res)


## @knitr 
locus18_pop1 <- c(gp_res$pop_alleles[[1]][[1]][,18], 
                  gp_res$pop_alleles[[2]][[1]][,18])
# sort alleles by size
allele_sort <- order(locus18_pop1, decreasing = FALSE)
#plot
plot(locus18_pop1[allele_sort], ylab = "allele size", col="blue",
     pch = 16)


## @knitr 
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


## @knitr eval = FALSE
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
## G_st_est        G_hed_st_est    D_Jost_est      Fst_WC
## 0.3905          0.8938          0.8256          0.4010
## 0.5519          0.8719          0.6986          0.6031
## 0.5924          0.8880          0.7092          0.6096
## ...             ...             ...             ...
## ...             ...             ...             ...
## 
## # these results could then be piped to further analyses or
## # visualisation tools


