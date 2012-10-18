################################################################################
################################################################################
##                              diveRsity v1.3.2                              ##  
##                            by Kevin Keenan QUB                             ##  
##            An R package for the calculation of differentiation             ##
##              statistics and locus informativeness statistics               ##  
##                V 1.2.0 and up allows parallel computations                 ##  
##                       copyright Kevin Keenan 2012                          ##  
################################################################################
################################################################################

# div.part, a wrapper function for the calculation of differentiation stats.
div.part<-function(infile, outfile = NULL, gp = 3, WC_Fst = FALSE, 
                   bs_locus = FALSE, bs_pairwise = FALSE, 
                   bootstraps = 0, Plot = FALSE, parallel = FALSE){
  ############################ Argument definitions ############################
  D<-infile
  on<-outfile
  gp<-gp
  fst<-WC_Fst
  bstrps<-bootstraps
  bsls<-bs_locus
  bspw<-bs_pairwise
  plt<-Plot
  para<-parallel
  
  ##############################################################################
  if(bsls==T && bstrps<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else if (bspw==T && bstrps<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else {
    #Use pre.div to calculate the standard global and locus stats
    pre_data_in<-list(infile=D,gp=gp,bootstrap="FALSE",locs="TRUE",fst=fst)
    accDat<-pre.divLowMemory(pre_data_in)
    # create a directory for output
    if(!is.null(on)){
      suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                             "-[diveRsity]","/",sep="")))
    }
    of=paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
    wd<-getwd()
    write_res<-is.element("xlsx",installed.packages()[,1])
    plot_res<-is.element("sendplot",installed.packages()[,1])
    if(Sys.info()["sysname"][[1]]=="Linux"){
      para_pack_inst<-is.element(c("snow","doSNOW","foreach","iterators"),
                                 installed.packages()[,1]) 
    } else {    
      para_pack_inst<-is.element(c("parallel","doParallel","foreach","iterators"),
                                 installed.packages()[,1])
    }
    
    if(plt == TRUE && is.null(on)){
      writeWarn <- paste("", "[NOTE]",
                         "Your results can't be plotted as you have not",
                         "provided an argument for 'outfile'.",
                         "Analysis completed", sep="\n")
      cat(noquote(writeWarn))
    }
    para_pack<-any(para_pack_inst==FALSE)
    if(write_res==F){
      Warning1<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please install the package 'xlsx' if you would like your", 
                       "results written to an Excel workbook.",
                       "Alternatively, your result will automatically be written",
                       "to .txt files.",
                       "___________________________________________________________",
                       "To install 'xlsx' use:",
                       "> install.packages('xlsx', dependencies=TRUE)",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning1))
    } 
    if(plot_res==F && plt==T){
      Warning2<-{paste(" "," "," ",
                       "[NOTE]  ",
                       "___________________________________________________________",
                       "Please install the package 'sendplot' to plot your results.",
                       "Use:",
                       "> install.packages('sendplot', dependencies = TRUE)",
                       "See:",
                       "> ?install.packages - for usage details",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning2))
    }
    if(fst==TRUE){
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est","Fis_WC","Fst_WC","Fit_WC")
    } else {
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est")
    }
    
    ############################################################################
    # output file
    # multilocus stats vector
    # pre output table for global locus stats
    #standard
    pre_ot1<-cbind(accDat$locus_names,round(as.numeric(accDat$hst),4),
                   round(as.numeric(accDat$dst),4),
                   round(as.numeric(accDat$gst),4),
                   round(as.numeric(accDat$gst_hedrick),4),
                   round(as.numeric(accDat$djost),4))
    # Add global multi locus stats to output table
    ot1<-rbind(pre_ot1,c("Global","","",accDat[[18]],accDat[[20]],accDat[[21]]))
    colnames(ot1)<-c("loci","H_st","D_st","G_st","G_hed_st","D_jost")
    #Estimated
    pre_ot2<-cbind(accDat$locus_names,
                   round(as.numeric(accDat$locus_harmonic_N),4),
                   round(as.numeric(accDat$ht_est),4),
                   round(as.numeric(accDat$dst_est),4),
                   round(as.numeric(accDat$gst_est),4),
                   round(as.numeric(accDat$gst_est_hedrick),4),
                   round(as.numeric(accDat$djost_est),4))
    ot2<-rbind(pre_ot2,c("Global","","","",accDat[[24]],accDat[[27]],
                         accDat[[28]]))
    colnames(ot2)<-c("loci","Harmonic_N","H_st_est","D_st_est","G_st_est",
                     "G_hed_st_est","D_Jost_est")
    if(fst==TRUE){
      ot2<-cbind(ot2, accDat$fstats[,1:3])
    }
    if(fst==TRUE){
      plot_data321<-c("Overall","","","",accDat$gst_est_all,
                      accDat$gst_est_all_hedrick,
                      accDat$djost_est_all,
                      as.numeric(accDat$fstats[(accDat$nloci+1),2]))
    } else {
      plot_data321<-c("Overall","","","",accDat$gst_est_all,
                      accDat$gst_est_all_hedrick,
                      accDat$djost_est_all)
    }
    if (!is.null(on)){
      if(write_res==TRUE){
        # write data to excel
        # Load dependencies
        require("xlsx")
        # standard stats
        write.xlsx(ot1,file=paste(of,"[div.part].xlsx",sep=""),
                   sheetName="Standard_stats",col.names=T,
                   row.names=F,append=F)
        # Estimated stats
        write.xlsx(ot2,file=paste(of,"[div.part].xlsx",sep=""),
                   sheetName="Estimated_stats",col.names=T,
                   row.names=F,append=T)
      } else {
        # text file alternatives
        std<-file(paste(of,"Standard-stats[div.part].txt",sep=""), "w")
        cat(paste(colnames(ot1),sep=""),"\n",sep="\t",file=std)
        for(i in 1:nrow(ot1)){
          cat(ot1[i,],"\n",file=std,sep="\t")
        }
        close(std)
        est<-file(paste(of,"Estimated-stats[div.part].txt",sep=""),"w")
        cat(paste(colnames(ot2),sep=""),"\n",sep="\t",file=est)
        for(i in 1:nrow(ot2)){
          cat(ot2[i,],"\n",file=est,sep="\t")
        }
        close(est)
      }
    }
    ot1out<-ot1[,-1]
    ot2out<-ot2[,-1]
    
    ot1out<-matrix(as.numeric(ot1[,2:6]),ncol=5)
    rownames(ot1out)<-ot1[,1]
    colnames(ot1out)<-colnames(ot1)[-1]
    
    ot2out<-matrix(as.numeric(ot2[,-1]),ncol=(ncol(ot2)-1))
    rownames(ot2out)<-ot2[,1]
    colnames(ot2out)<-colnames(ot2)[-1]
    if (para == TRUE && para_pack == TRUE){
      if(Sys.info()["sysname"][[1]]=="Linux"){
        Warning3<-{paste(" "," ",
                         "[NOTE]",
                         "___________________________________________________________",
                         "Please make sure the packages 'doSNOW', 'snow', 'foreach'",
                         " and 'iterators' are installed. These are required to run",
                         " your analysis in parallel.",
                         "Your analysis will be run sequentially!",
                         "___________________________________________________________",
                         "To install these use:",
                         "> install.packages()",
                         "See:",
                         "> ?install.packages - for usage details.",
                         "___________________________________________________________",
                         sep="\n")
        }
      } else {
        Warning3<-{paste(" "," ",
                         "[NOTE]",
                         "___________________________________________________________",
                         "Please make sure the packages 'parallel', 'doParallel',",
                         "'foreach' and 'iterators' are installed. These are required",
                         " to run your analysis in parallel.",
                         "Your analysis will be run sequentially!",
                         "___________________________________________________________",
                         "To install these use:",
                         "> install.packages()",
                         "See:",
                         "> ?install.packages - for usage details.",
                         "___________________________________________________________",
                         sep="\n")
        }
      }
      cat(noquote(Warning3))
    }
    
    ############################################################################
    ############################ Bootstrapper ##################################
    ############################################################################
    if (para == TRUE && para_pack == FALSE) {
      #count cores
      if(any((.packages())=="snow") && any((.packages())=="doSNOW")){
        detach(package:doSNOW)
        detach(package:snow)
        library(parallel)
        cores<-detectCores(logical=TRUE)
        if(any((.packages())=="doParallel")){
          detach(package:doParallel)
          detach(package:parallel)
        } else {
          detach(package:parallel)
        }
        library(doSNOW)
        library(snow)
      } else {
        library(parallel)
        cores<-detectCores(logical=TRUE)
        if(any((.packages())=="doParallel")){
          detach(package:doParallel)
          detach(package:parallel)
        } else {
          detach(package:parallel)
        }
      }
      if(Sys.info()["sysname"][[1]]=="Linux"){
        library("doSNOW")
        #make clusters
        cl<-makeCluster(cores,type="SOCK")
        registerDoSNOW(cl)         # Make sure 'cl' is closed
      } else if (Sys.info()["sysname"][[1]]!="Linux"){
        library("doParallel")
        cl<-makeCluster(cores)
        registerDoParallel(cl)
      }
      note<-paste( "[NOTE]  ",
                   "Cores successfully registered for parallel computations...",
                   " ",
                   sep="\n")
      cat(noquote(note))
    }
    
    # Used only if bootstraps is greater than zero
    if(bsls==T){
      
      if (para == TRUE && para_pack == FALSE) {
        
        #vectorize prallele#
        gp_inls<-list(infile=D,gp=gp,bootstrap="TRUE",locs="TRUE",fst=fst)
        gp_in<-list()
        for(i in 1:bstrps){
          gp_in[[i]]<-gp_inls
        }
        
        # calculate stats from readGenepop objects
        bs_loc<-parLapply(cl,gp_in, pre.divLowMemory)
        rm(gp_in)                          ###
        z<-gc(reset=T, verbose=FALSE)        ## tidy up
        rm(z)                              ###
        #vectorize data extraction#
        if(fst==TRUE){
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4),
              as.numeric(bs_loc[[x]]$fstats[(accDat$nloci+1),1:3]))
          }))
        }else{
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4))
          }))
        }
        bs_std<-lapply(1:accDat$nloci, function(x){
          do.call("rbind",lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x],4),
              round(bs_loc[[y]]$gst_hedrick[x],4),
              round(bs_loc[[y]]$djost[x],4))}))
        })
        if(fst==TRUE){
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4),
                as.numeric(bs_loc[[y]]$fstats[x,1:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)                  ###
        z<-gc(reset=T)                ### tidy up
        rm(z)                       ###
        
      } else {
        #vectorize non-parallel#
        
        gp_inls<-list(infile=D,gp=gp,bootstrap="TRUE",locs="TRUE",fst=fst)
        gp_in<-list()
        for(i in 1:bstrps){
          gp_in[[i]]<-gp_inls
        }
        # calculate stats from readGenepop objects
        bs_loc<-lapply(gp_in,pre.divLowMemory)
        rm(gp_in)                          ###
        z<-gc(reset=T, verbose=FALSE)        ## tidy up
        rm(z)                              ###
        if(fst==TRUE){
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4),
              as.numeric(bs_loc[[x]]$fstats[(accDat$nloci+1),1:3]))
          }))
        }else{
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4))
          }))
        }
        bs_std<-lapply(1:accDat$nloci, function(x){
          do.call("rbind",lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x],4),
              round(bs_loc[[y]]$gst_hedrick[x],4),
              round(bs_loc[[y]]$djost[x],4))}))
        })
        if(fst==TRUE){
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4),
                as.numeric(bs_loc[[y]]$fstats[x,1:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)
        z<-gc(reset=T)
        rm(z)
        
      }
      
      
      #vectorize#
      if(fst==TRUE){
        bs_res<-lapply(1:9,function(x){matrix(ncol=3, nrow=(accDat$nloci+1))})
      } else {
        bs_res<-lapply(1:6,function(x){matrix(ncol=3, nrow=(accDat$nloci+1))})
      }
      bs_join<-cbind(bs_std, bs_est)
      ciCalc<-function(x){
        res<-lapply(x, function(y){
          ci<-function(x){
            sd(na.omit(x))*1.96
          }
          apply(y,2,ci)
        })
        return(c(res[[1]],res[[2]]))
      }
      ci<-function(x){
        sd(na.omit(x))*1.96
      }
      bs_cis<-t(apply(bs_join, 1, ciCalc))
      bs_cis<-rbind(bs_cis, apply(bs_glb,2,ci))
      if(fst==TRUE){
        for(i in 1:9){
          if(i <= 3){
            bs_res[[i]][,1]<-as.numeric(ot1[,(i+3)])
            bs_res[[i]][,2]<-round(as.numeric(ot1[,(i+3)])-bs_cis[,i],4)
            bs_res[[i]][,3]<-round(as.numeric(ot1[,(i+3)])+bs_cis[,i],4)
          } else {
            bs_res[[i]][,1]<-as.numeric(ot2[,(i+1)])
            bs_res[[i]][,2]<-round(as.numeric(ot2[,(i+1)])-bs_cis[,i],4)
            bs_res[[i]][,3]<-round(as.numeric(ot2[,(i+1)])+bs_cis[,i],4)
          }
          bs_res[[i]][is.na(bs_res[[i]])]<-0
        }
      } else {
        for(i in 1:6){
          if(i <= 3){
            bs_res[[i]][,1]<-as.numeric(ot1[,(i+3)])
            bs_res[[i]][,2]<-round(as.numeric(ot1[,(i+3)])-bs_cis[,i],4)
            bs_res[[i]][,3]<-round(as.numeric(ot1[,(i+3)])+bs_cis[,i],4)
          } else {
            bs_res[[i]][,1]<-as.numeric(ot2[,(i+1)])
            bs_res[[i]][,2]<-round(as.numeric(ot2[,(i+1)])-bs_cis[,i],4)
            bs_res[[i]][,3]<-round(as.numeric(ot2[,(i+1)])+bs_cis[,i],4)
          }
          bs_res[[i]][is.na(bs_res[[i]])]<-0
        }
      }
      
      names(bs_res)<-namer
      
      bs_res1<-bs_res
      if(fst==FALSE){
        for(i in 1:6){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names,"global"),
                                       c("Actual","Lower_CI","Upper_CI"))
        }
      } else if (fst==TRUE){
        for(i in 1:9){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names,"global"),
                                       c("Actual","Lower_CI","Upper_CI"))
        }
      }
      # bs results output object header
      hdr<-matrix(c("locus","Actual","Lower_95%CI","Upper_95%CI"),ncol=4)
      bs_out<-matrix(rbind(hdr,c(names(bs_res)[1],"","",""),
                           cbind(c(accDat$locus_names,"Overall"),
                                 bs_res[[1]])),ncol=4)
      if(fst==TRUE){
        for(i in 2:9){
          bs_out<-matrix(rbind(bs_out,c(names(bs_res)[i],"","",""),
                               cbind(c(accDat$locus_names,"Global"),
                                     bs_res[[i]])),ncol=4)
        }
      } else {
        for(i in 2:6){
          bs_out<-matrix(rbind(bs_out,c(names(bs_res)[i],"","",""),
                               cbind(c(accDat$locus_names,"Global"),
                                     bs_res[[i]])),ncol=4)
        }
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(bs_out,file=paste(of,"[div.part].xlsx",sep=""),
                     sheetName="Locus_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          bts<-file(paste(of,"Locus-bootstrap[div.part].txt",sep=""), "w")
          cat(paste(colnames(bs_out),sep=""),"\n",sep="\t",file=bts)
          for(i in 1:nrow(bs_out)){
            cat(bs_out[i,],"\n",file=bts,sep="\t")
          }
          close(bts)
        }
      }
    }
    zzz<-gc()
    rm(zzz)
    if(plot_res==TRUE && plt==TRUE && bsls==TRUE){
      
      #vectorize#
      sorter<-function(x){
        z<-order(x[1:accDat$nloci,1],decreasing=F)
        if(length(z) >= 200){
          z<-z[(length(z)-150):length(z)]
        }
        return(z)
      }
      lso123<-lapply(bs_res, sorter)
      
      #
      names(lso123)<-namer
      plot.call_loci<-list()
      plot.extras_loci<-list()
      xy.labels_loci<-list()
      y.pos_loci<-list()
      x.pos_loci=1:accDat$nloci
      direct=of
      fn_pre_loci<-list()
      #Plot Gst_Nei
      plot.call_loci[[1]]=c("plot(bs_res[[4]][lso123[[4]],1],
                          ylim=c(0,(max(bs_res[[4]][,3])+
                            min(bs_res[[4]][,3]))),xaxt='n',
                            ylab=names(bs_res)[4],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_loci[[1]]=c("points(bs_res[[4]][lso123[[4]],1],
                            pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[4]][lso123[[4]],2],
                              1:accDat$nloci,bs_res[[4]][lso123[[4]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[4]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
  
      xy.labels_loci[[1]]=data.frame(Locus_name=accDat$locus_names[lso123[[4]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[4]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[4]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[4]],1],4))
      
      y.pos_loci[[1]]=bs_res[[4]][lso123[[4]],1]
      fn_pre_loci[[1]]<-names(bs_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_loci[[2]]=c("plot(bs_res[[5]][lso123[[5]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(bs_res)[5],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_loci[[2]]=c("points(bs_res[[5]][lso123[[5]],1],
                            pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[5]][lso123[[5]],2],
                              1:accDat$nloci,bs_res[[5]][lso123[[5]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[5]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
  
      xy.labels_loci[[2]]=data.frame(Locus_name=accDat$locus_names[lso123[[5]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[5]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[5]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[5]],1],4))
      
      y.pos_loci[[2]]=bs_res[[5]][lso123[[5]],1]
      fn_pre_loci[[2]]<-names(bs_res)[5]
      
      
      # Plot D_jost
      plot.call_loci[[3]]=c("plot(bs_res[[6]][lso123[[6]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(bs_res)[6],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_loci[[3]]=c("points(bs_res[[6]][lso123[[6]],1],
                            pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[6]][lso123[[6]],2],
                              1:accDat$nloci,bs_res[[6]][lso123[[6]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[6]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
  
      xy.labels_loci[[3]]=data.frame(Locus_name=accDat$locus_names[lso123[[6]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[6]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[6]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[6]],1],4))
      
      y.pos_loci[[3]]=bs_res[[6]][lso123[[6]],1]
      fn_pre_loci[[3]]<-names(bs_res)[6]
      
      #plot(Fst)
      if(fst==TRUE){
        plot.call_loci[[4]]=c("plot(bs_res[[8]][lso123[[8]],1],
                          ylim=c(0,(max(bs_res[[8]][,3])+
                              min(bs_res[[8]][,3]))),xaxt='n',
                              ylab=names(bs_res)[8],type='n',
                              xlab='Loci \n (Hover over a point to see locus data)',
                              cex.lab=1.5,cex.axis=1.3,las=1)")
        
        plot.extras_loci[[4]]=c("points(bs_res[[8]][lso123[[8]],1],
                            pch=15,col='black',cex=1);
                                arrows(1:accDat$nloci,bs_res[[8]][lso123[[8]],2],
                                1:accDat$nloci,bs_res[[8]][lso123[[8]],3],code=3,
                                angle=90,length=0.05,lwd=0.1);
                                abline(h=c(0,bs_res[[8]][(accDat$nloci+1),2]),
                                lwd=1,lty=c(1,2),col=c('black','red'))")
  
        xy.labels_loci[[4]]=data.frame(Locus_name=accDat$locus_names[lso123[[8]]],
                                       Gst_Nei=round(bs_res[[4]][lso123[[8]],1],4),
                                       Gst_Hedrick=round(bs_res[[5]][lso123[[8]],1],4),
                                       D_jost=round(bs_res[[6]][lso123[[8]],1],4),
                                       Fst_WC=round(bs_res[[8]][lso123[[8]],1],4))
        
        y.pos_loci[[4]]=bs_res[[8]][lso123[[8]],1]
        fn_pre_loci[[4]]<-names(bs_res)[8]
      }
    }
    ############################################################################
    ################################## Pairwise ################################
    ############################################################################
    # population pair combinations
    pw<-combn(accDat$npops,2)
    pwmat<-pw+1
    #pw data creator
    ind_vectors<-lapply(1:accDat$npops, function(x){
      rep(x, accDat$pop_sizes[[x]])}
    )
    #      
    pre_data<-matrix(rep("",((accDat$nloci+1)*(accDat$nloci+1))),
                     ncol=(accDat$nloci+1))
    pre_data[1,]<-rep("",(accDat$nloci+1))
    #
    for(i in 2:(accDat$nloci+1)){
      pre_data[i,1]<-accDat$locus_names[(i-1)]
    }
    #
    pw_data<-list()
    for (i in 1:ncol(pw)){
      pw_data[[i]]<-data.frame(rbind(pre_data,
                                     c("POP",as.vector(rep("",accDat$nloci))),
                                     cbind(ind_vectors[[pw[1,i]]],
                                           matrix(noquote(accDat$pop_list
                                                          [[pw[1,i]]]),
                                                  ncol=accDat$nloci)),
                                     c("POP",as.vector(rep("",accDat$nloci))),
                                     cbind(ind_vectors[[pw[2,i]]],
                                           matrix(noquote(accDat$pop_list
                                                          [[pw[2,i]]]),
                                                  ncol=accDat$nloci))))
    }
    true_stat_gp_in<-list()
    if(fst==TRUE){
      pw_glb<-matrix(rep(0,(9*(ncol(pw)))),ncol=9)
    } else {
      pw_glb<-matrix(rep(0,(6*(ncol(pw)))),ncol=6)
    }
    for (i in 1:ncol(pw)){
      true_stat_gp_in[[i]]<-list(infile=pw_data[[i]],gp=gp,bootstrap="FALSE",
                                 locs="FALSE",fst=fst)
    }
    if (para == TRUE && para_pack == FALSE) {
      true_stat<-parLapply(cl,true_stat_gp_in, pre.divLowMemory)
      # close core connections if not needed further
      if (bspw==FALSE){
        stopCluster(cl)
      }
    } else {
      true_stat<-lapply(true_stat_gp_in, pre.divLowMemory)
    }
    for(i in 1:ncol(pw)){
      if(fst==TRUE){
        pw_glb[i,]<-c(true_stat[[i]]$gst_all,true_stat[[i]]$gst_all_hedrick,
                      true_stat[[i]]$djost_all,true_stat[[i]]$gst_est_all,
                      true_stat[[i]]$gst_est_all_hedrick,
                      true_stat[[i]]$djost_est_all,
                      as.numeric(true_stat[[i]]$fstat))
      } else {
        pw_glb[i,]<-c(true_stat[[i]]$gst_all,true_stat[[i]]$gst_all_hedrick,
                      true_stat[[i]]$djost_all,true_stat[[i]]$gst_est_all,
                      true_stat[[i]]$gst_est_all_hedrick,
                      true_stat[[i]]$djost_est_all)
      }
      
      true_stat[[i]]<-0
    }
    if(fst==TRUE){
      pwMatList<-lapply(1:9, function(x){
        noquote(matrix(rep("--",((accDat$npops+1)^2)),ncol=(accDat$npops+1),
                       nrow=(accDat$npops+1)))
      })
    } else {
      pwMatList<-lapply(1:6, function(x){
        noquote(matrix(rep("--",((accDat$npops+1)^2)),ncol=(accDat$npops+1),
                       nrow=(accDat$npops+1)))
      })
    }
    if(fst==TRUE){
      pwMatListOut<-lapply(1:9, function(x){
        noquote(matrix(rep(NA,((accDat$npops)^2)),ncol=(accDat$npops),
                       nrow=(accDat$npops)))
      })
    } else {
      pwMatListOut<-lapply(1:6, function(x){
        noquote(matrix(rep(NA,((accDat$npops)^2)),ncol=(accDat$npops),
                       nrow=(accDat$npops)))
      })
    }
    names(pwMatList)<-namer
    names(pwMatListOut)<-namer
    #write pw res to matrices
    pnames<-c("", accDat$pop_names)
    pnamesOut<-accDat$pop_names
    if(fst==TRUE){
      for(i in 1:9){
        for(j in 1:ncol(pw)){
          pwMatList[[i]][pwmat[2,j],pwmat[1,j]]<-pw_glb[j,i]
          pwMatList[[i]][pwmat[1,j],pwmat[2,j]]<-""
          pwMatListOut[[i]][pw[2,j],pw[1,j]]<-pw_glb[j,i]
          #pwMatListOut[[i]][pw[1,j],pw[2,j]]<-""
        }
        pwMatList[[i]][1,]<-pnames
        pwMatList[[i]][,1]<-pnames
        dimnames(pwMatListOut[[i]])<-list(pnamesOut,pnamesOut)
      }
    } else {
      for(i in 1:6){
        for(j in 1:ncol(pw)){
          pwMatList[[i]][pwmat[2,j],pwmat[1,j]]<-pw_glb[j,i]
          pwMatList[[i]][pwmat[1,j],pwmat[2,j]]<-""
          pwMatListOut[[i]][pw[2,j],pw[1,j]]<-pw_glb[j,i]
          #pwMatListOut[[i]][pw[1,j],pw[2,j]]<-""
        }
        pwMatList[[i]][1,]<-pnames
        pwMatList[[i]][,1]<-pnames
        dimnames(pwMatListOut[[i]])<-list(pnamesOut,pnamesOut)
      }
    }
    
    
    # write object create
    #pnames list
    
    pwWrite<-pwMatList[[1]]
    pwWrite<-rbind(c(names(pwMatList)[1],rep("",accDat$npops)),pwWrite, 
                   rep("",(accDat$npops+1)))
    if(fst==TRUE){
      for(i in 2:9){
        pwWrite<-rbind(pwWrite,c(names(pwMatList)[i],rep("",accDat$npops)),
                       pwMatList[[i]],rep("",(accDat$npops+1)))
      }
    } else {
      for(i in 2:6){
        pwWrite<-rbind(pwWrite,c(names(pwMatList)[i],rep("",accDat$npops)),
                       pwMatList[[i]],rep("",(accDat$npops+1)))
      }
    }
    if(!is.null(on)){
      if(write_res==TRUE){
        # write data to excel
        # Load dependencies
        # pw stats
        write.xlsx(pwWrite,file=paste(of,"[div.part].xlsx",sep=""),
                   sheetName="Pairwise-stats",col.names=F,
                   row.names=F,append=T)
      } else {
        # text file alternatives
        pw_outer<-file(paste(of,"Pairwise-stats[div.part].txt",sep=""), "w")
        for(i in 1:nrow(pwWrite)){
          cat(pwWrite[i,],"\n",file=pw_outer,sep="\t")
        }
        close(std)
      }
    }
    #cleanup
    rm("pwWrite")
    ##
    zzz<-gc()
    rm(zzz)
    #Bootstrap
    if(bspw==TRUE){
      
      # Bootstrap results data object 
      # bs_pw_glb = bootstrap pairwise global stats
      if(fst==TRUE){
        bs_pw_glb<-matrix(rep(0,(9*bstrps)),ncol=9,nrow=bstrps)
      } else {
        bs_pw_glb<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
      }
      # output results data object
      # pw_res = pairwise results
      if(fst==TRUE){
        pw_res<-lapply(1:9, function(x){
          matrix(nrow=ncol(pw), ncol=3)
        })
      } else {
        pw_res<-lapply(1:6, function(x){
          matrix(nrow=ncol(pw), ncol=3)
        })
      }
      #
      #
      
      #parallel processing option
      if (para == TRUE && para_pack == FALSE) {
        #create a readGenepop list
        bs_pw_glb<-list()
        data_res<-list()
        bs_pw_para<-list()
        for(i in 1:ncol(pw)){
          input<-list(infile=pw_data[[i]],gp=gp,bootstrap="TRUE",
                      locs="FALSE",fst=fst)
          pw_inlist<-list()
          for(j in 1:bstrps){
            pw_inlist[[j]]<-input
          }
          if(fst==TRUE){
            bs_pw_glb[[i]]<-matrix(rep(0,(9*bstrps)),ncol=9,nrow=bstrps)
          } else {
            bs_pw_glb[[i]]<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
          }
          bs_pw_para<-parLapply(cl,pw_inlist,pre.divLowMemory)
          for(j in 1:bstrps){
            if(fst==TRUE){
              bs_pw_glb[[i]][j,]<-c(bs_pw_para[[j]]$gst_all,
                                    bs_pw_para[[j]]$gst_all_hedrick,
                                    bs_pw_para[[j]]$djost_all,
                                    bs_pw_para[[j]]$gst_est_all,
                                    bs_pw_para[[j]]$gst_est_all_hedrick,
                                    bs_pw_para[[j]]$djost_est_all,
                                    as.numeric(bs_pw_para[[j]]$fstats))
            } else {
              bs_pw_glb[[i]][j,]<-c(bs_pw_para[[j]]$gst_all,
                                    bs_pw_para[[j]]$gst_all_hedrick,
                                    bs_pw_para[[j]]$djost_all,
                                    bs_pw_para[[j]]$gst_est_all,
                                    bs_pw_para[[j]]$gst_est_all_hedrick,
                                    bs_pw_para[[j]]$djost_est_all)
            }
          }
        }
        #
        # confidence interval calculator function
        ci<-function(x){
          ci_raw<-NULL
          for(i in 1:ncol(x)){
            ci_raw[i]<- 1.96*(sd(x[,i]))
          }
          return(ci_raw)
        }
        # Calculate confidence interval  
        cis<-parLapply(cl,bs_pw_glb,ci)
        #stopCluster(cl)
        
        for(i in 1:ncol(pw)){
          if(fst==TRUE){
            for(j in 1:9){
              pw_res[[j]][i,1]<-pw_glb[i,j]
              pw_res[[j]][i,2]<-round((pw_glb[i,j]-cis[[i]][j]),4)
              pw_res[[j]][i,3]<-round((pw_glb[i,j]+cis[[i]][j]),4)
              pw_res[[j]][is.na(pw_res[[j]])]<-0
            }
          } else {
            for(j in 1:6){
              pw_res[[j]][i,1]<-pw_glb[i,j]
              pw_res[[j]][i,2]<-round((pw_glb[i,j]-cis[[i]][j]),4)
              pw_res[[j]][i,3]<-round((pw_glb[i,j]+cis[[i]][j]),4)
              pw_res[[j]][is.na(pw_res[[j]])]<-0
            }
          }
        }
        stopCluster(cl)
      } else {
        #sequential vectorized
        pw_inlist<-list()
        for(i in 1:ncol(pw)){
          input<-list(infile=pw_data[[i]],gp=gp,bootstrap="TRUE",locs="FALSE",
                      fst=fst)
          pw_inlist[[i]]<-list()
          for(j in 1:bstrps){
            pw_inlist[[i]][[j]]<-input
          }
        }
        bs_pw_glb<-list()
        for(i in 1:ncol(pw)){
          if(fst==TRUE){
            bs_pw_glb[[i]]<-matrix(rep(0,(9*bstrps)),ncol=9,nrow=bstrps)
          } else {
            bs_pw_glb[[i]]<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
          }
        }
        #create a readGenepop list
        bs_pw_glb<-list()
        data_res<-list()
        bs_pw_para<-list()
        for(i in 1:ncol(pw)){
          input<-list(infile=pw_data[[i]],gp=gp,bootstrap="TRUE",
                      locs="FALSE",fst=fst)
          pw_inlist<-list()
          for(j in 1:bstrps){
            pw_inlist[[j]]<-input
          }
          if(fst==TRUE){
            bs_pw_glb[[i]]<-matrix(rep(0,(9*bstrps)),ncol=9,nrow=bstrps)
          } else {
            bs_pw_glb[[i]]<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
          }
          bs_pw_para<-lapply(pw_inlist,pre.divLowMemory)
          for(j in 1:bstrps){
            if(fst==TRUE){
              bs_pw_glb[[i]][j,]<-c(bs_pw_para[[j]]$gst_all,
                                    bs_pw_para[[j]]$gst_all_hedrick,
                                    bs_pw_para[[j]]$djost_all,
                                    bs_pw_para[[j]]$gst_est_all,
                                    bs_pw_para[[j]]$gst_est_all_hedrick,
                                    bs_pw_para[[j]]$djost_est_all,
                                    as.numeric(bs_pw_para[[j]]$fstat))
            } else {
              bs_pw_glb[[i]][j,]<-c(bs_pw_para[[j]]$gst_all,
                                    bs_pw_para[[j]]$gst_all_hedrick,
                                    bs_pw_para[[j]]$djost_all,
                                    bs_pw_para[[j]]$gst_est_all,
                                    bs_pw_para[[j]]$gst_est_all_hedrick,
                                    bs_pw_para[[j]]$djost_est_all)
            }
          }
        } 
        # confidence interval calculator function
        ci<-function(x){
          ci_raw<-NULL
          for(i in 1:ncol(x)){
            ci_raw[i]<- 1.96*(sd(na.omit(x[,i])))
          }
          return(ci_raw)
        }
        # Calculate confidence interval
        cis<-lapply(bs_pw_glb,ci)
        for(i in 1:ncol(pw)){
          if(fst==TRUE){
            for(j in 1:9){
              pw_res[[j]][i,1]<-pw_glb[i,j]
              pw_res[[j]][i,2]<-round((pw_glb[i,j]-cis[[i]][j]),4)
              pw_res[[j]][i,3]<-round((pw_glb[i,j]+cis[[i]][j]),4)
            }
          } else {
            for(j in 1:6){
              pw_res[[j]][i,1]<-pw_glb[i,j]
              pw_res[[j]][i,2]<-round((pw_glb[i,j]-cis[[i]][j]),4)
              pw_res[[j]][i,3]<-round((pw_glb[i,j]+cis[[i]][j]),4)
            }
          }
        }
        #
      }
      #
      # pairwise comparisons
      # pw_names = pairwise population names
      pw_nms<-paste(accDat$pop_names[pw[1,]],
                    accDat$pop_names[pw[2,]],sep=" vs. ")
      #
      pw_nms1<-paste(pw[1,],pw[2,],sep=" vs. ")
      #
      names(pw_res)<-namer
      #
      pw_res1<-pw_res
      if(fst==TRUE){
        for(i in 1:9){
          dimnames(pw_res1[[i]])<-list(pw_nms,c("Actual","Lower_CI","Upper_CI"))
        }
      } else {
        for(i in 1:6){
          dimnames(pw_res1[[i]])<-list(pw_nms,c("Actual","Lower_CI","Upper_CI"))
        }
      }
      # bs results output object header
      hdr<-matrix(c("Pairwise","Actual","Lower_95%CI","Upper_95%CI"),ncol=4)
      pw_bs_out<-matrix(rbind(hdr,c(names(pw_res)[1],"","",""),
                              cbind(pw_nms,pw_res[[1]])),ncol=4)
      if(fst==TRUE){
        for(i in 2:9){
          pw_bs_out<-matrix(rbind(pw_bs_out,c(names(pw_res)[i],"","",""),
                                  cbind(pw_nms,pw_res[[i]])),ncol=4)
        }
      } else {
        for(i in 2:6){
          pw_bs_out<-matrix(rbind(pw_bs_out,c(names(pw_res)[i],"","",""),
                                  cbind(pw_nms,pw_res[[i]])),ncol=4)
        }
      }
      if(is.null(on)){
        if(write_res==TRUE){
          write.xlsx(pw_bs_out,file=paste(of,"[div.part].xlsx",sep=""),
                     sheetName="Pairwise_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          pw_bts<-file(paste(of,"Pairwise-bootstrap[div.part].txt",sep=""), "w")
          cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
          for(i in 1:nrow(pw_bs_out)){
            cat(pw_bs_out[i,],"\n",file=pw_bts,sep="\t")
          }
          close(pw_bts)
        }
      }
    }
    zzz<-gc()
    rm(zzz)
    ############################################################################
    #pw plotter
    if(plot_res==TRUE && plt==TRUE && bspw==TRUE){
      pwso<-list()
      for(i in 1:length(pw_res)){
        pwso[[i]]<-order(pw_res[[i]][,1],decreasing=F)
        if(length(pwso[[i]]) >= 100){
          pwso[[i]]<-pwso[[i]][(length(pwso[[i]])-99):length(pwso[[i]])]
        }
      }
      names(pwso)<-namer
      # define plot parameters 
      plot.call_pw<-list()
      plot.extras_pw<-list()
      xy.labels_pw<-list()
      y.pos_pw<-list()
      x.pos_pw=1:length(pwso[[i]])
      fn_pre_pw<-list()
      direct=of
      #Plot Gst_Nei
      plot.call_pw[[1]]=c("plot(pw_res[[4]][pwso[[4]],1],
                            ylim=c(0,(max(pw_res[[4]][,3])+
                          min(pw_res[[4]][,3]))),xaxt='n',
                          ylab=names(pw_res)[4],type='n',
                          xlab='Pairwise comparisons 
                              \n (Hover over a point to see pairwise info.)',
                        cex.lab=1.2,cex.axis=1.3,las=1)")

      plot.extras_pw[[1]]=c("points(pw_res[[4]][pwso[[4]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],2],
                            1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[5]),
                            lwd=1,lty=2,col='red')")
  
      xy.labels_pw[[1]]=data.frame(pairwise_name=pw_nms[pwso[[4]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[4]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[4]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[4]],1],4))
      
      y.pos_pw[[1]]=pw_res[[4]][pwso[[4]],1]
      fn_pre_pw[[1]]<-names(pw_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_pw[[2]]=c("plot(pw_res[[5]][pwso[[5]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(pw_res)[5],type='n',
                          xlab='Pairwise comparisons
                              \n (Hover over a point to see pairwise info.)',
                        cex.lab=1.2,cex.axis=1.3,las=1)")

      plot.extras_pw[[2]]=c("points(pw_res[[5]][pwso[[5]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[5]]),pw_res[[5]][pwso[[5]],2],
                            1:length(pwso[[5]]),pw_res[[5]][pwso[[5]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[6]),
                            lwd=1,lty=2,col='red')")
  
      xy.labels_pw[[2]]=data.frame(pairwise_name=pw_nms[pwso[[5]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[5]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[5]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[5]],1],4))
      
      y.pos_pw[[2]]=pw_res[[5]][pwso[[5]],1]
      fn_pre_pw[[2]]<-names(pw_res)[5]
      
      
      # Plot D_jost
      plot.call_pw[[3]]=c("plot(pw_res[[6]][pwso[[6]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(pw_res)[6],type='n',
                          xlab='Pairwise comparisons 
                             \n (Hover over a point to see pairwise info.)',
                        cex.lab=1.2,cex.axis=1.3,las=1)")

      plot.extras_pw[[3]]=c("points(pw_res[[6]][pwso[[6]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[6]]),pw_res[[6]][pwso[[6]],2],
                            1:length(pwso[[6]]),pw_res[[6]][pwso[[6]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[7]),
                            lwd=1,lty=2,col='red')")
    
      xy.labels_pw[[3]]=data.frame(pairwise_name=pw_nms[pwso[[6]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[6]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[6]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[6]],1],4))
      
      y.pos_pw[[3]]=pw_res[[6]][pwso[[6]],1]
      fn_pre_pw[[3]]<-names(pw_res)[6]
      #plot(Fst_WC)
      if(fst==TRUE){
        plot.call_pw[[4]]=c("plot(pw_res[[8]][pwso[[8]],1],
                        ylim=c(0,(max(pw_res[[8]][,3])+
                            min(pw_res[[8]][,3]))),xaxt='n',ylab=names(pw_res)[8],type='n',
                            xlab='Pairwise comparisons 
                        \n (Hover over a point to see pairwise info.)',
                        cex.lab=1.2,cex.axis=1.3,las=1)")

        plot.extras_pw[[4]]=c("points(pw_res[[8]][pwso[[8]],1],
                          pch=15,col='black',cex=1);
                              arrows(1:length(pwso[[8]]),pw_res[[8]][pwso[[8]],2],
                              1:length(pwso[[8]]),pw_res[[8]][pwso[[8]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=as.numeric(plot_data321[7]),
                              lwd=1,lty=2,col='red')")
    
        xy.labels_pw[[4]]=data.frame(pairwise_name=pw_nms[pwso[[8]]],
                                     Gst_Nei=round(pw_res[[4]][pwso[[8]],1],4),
                                     Gst_Hedrick=round(pw_res[[5]][pwso[[8]],1],4),
                                     D_jost=round(pw_res[[6]][pwso[[8]],1],4),
                                     Fst_WC=round(pw_res[[8]][pwso[[8]],1],4))
        
        y.pos_pw[[4]]=pw_res[[8]][pwso[[8]],1]
        fn_pre_pw[[4]]<-names(pw_res)[8]
      }
    }
    ############################### Bootstrap end ################################
    
    
    ################################# Plot resuts ################################
    #make necessary data available
    if(plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==TRUE){
      pl<-list(bs_res=bs_res,
               pw_res=pw_res,
               accDat=accDat,
               lso123=lso123,
               pwso=pwso,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="TRUE",
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               x.pos_pw=x.pos_pw,
               pw=pw,
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==FALSE){
      pl<-list(bs_res=bs_res,
               accDat=accDat,
               lso123=lso123,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="FALSE",
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==FALSE && bspw==TRUE){
      pl<-list(pw_res=pw_res,
               accDat=accDat,
               pwso=pwso,
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               x.pos_pw=x.pos_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               direct=direct,
               plot_loci="FALSE",
               plot_pw="TRUE",
               pw=pw,plot_data321=plot_data321,
               fst=fst)
    }
    if(!is.null(on)){
      if (plt==TRUE && plot_res==TRUE){
        suppressWarnings(plotter(x=pl,img="1000x600"))
      }
    }
    zzz<-gc()
    rm(zzz)
    
    
    #############################################################################
    #Data for output
    if(bspw==T && bsls==T){
      list(standard=ot1out,
           estimate=ot2out,
           pairwise=pwMatListOut,
           bs_locus=bs_res1,
           bs_pairwise=pw_res1)
    } else if(bspw==T && bsls==F){
      list(standard=ot1out,
           estimate=ot2out,
           pairwise=pwMatListOut,
           bs_pairwise=pw_res1)
    } else if(bspw==F && bsls==T){
      list(standard=ot1out,
           estimate=ot2out,
           pairwise=pwMatListOut,
           bs_locus=bs_res1)
    } else if(bspw==F && bsls==F){
      list(standard=ot1out,
           estimate=ot2out,
           pairwise=pwMatListOut)
    }
  }
}

################################################################################
# div.part end                                                                 #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# readGenepop, a function for the generation of basic population parameters    #
################################################################################
readGenepop<- function (x) {
  gp=x$gp
  infile=x$infile
  bootstrap=x$bootstrap
  locs=x$locs
  if(typeof(infile)=="list"){
    data1=infile 
  } else if (typeof(infile)=="character"){
    no_col <- max(count.fields(infile))
    suppressWarnings(flth<-length(readLines(infile)))
    data1<-matrix(NA,ncol=no_col,nrow=flth)
    datt<-list()
    for(i in 1:flth){
      datt[[i]]<-scan(infile,skip=(i-1),nlines=1,fill=T, 
                      what=list(rep("\t",no_col)),quiet = T)
      datt[[i]]<-unlist(datt[[i]])
      if(is.element(",", datt[[i]])==TRUE){
        datt[[i]]<-datt[[i]][-(which(datt[[i]] == ","))]
      }
      if(length(datt[[i]]) < no_col){
        datt[[i]]<-c(datt[[i]],rep(NA, no_col-length(datt[[i]])))
      }
      data1[i,]<-datt[[i]]
    }
    data1<-as.data.frame(data1)
  }
  data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
  raw_data<-data1
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
  pop_weights<- 1/pop_sizes
  
  n_harmonic<-npops/sum(pop_weights)
  
  N<-pop_sizes
  
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                   2:(nloci+1)])
  }
  
  
  
  if (gp==3) {
    plMake<-function(x){
      return(matrix(sprintf("%06g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
    }
  } else if (gp==2) {
    plMake<-function(x){
      return(matrix(sprintf("%04g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
    }
  }
  pop_list<-lapply(pop_list, plMake)
  
  
  for(i in 1:npops){
    pop_list[[i]][pop_list[[i]]=="    NA"]<-NA
  }
  
  
  if(bootstrap == T){
    bs<-function(x){
      return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
    }
    pop_list<-lapply(pop_list, bs)
  }  
  
  ###vectorize loci_pop_sizes#####################################################
  
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
  pre_loci_pop_sizes<-lapply(pop_list,lps)#
  pls<-matrix(ncol=nloci,nrow=npops)#
  for(i in 1:length(pre_loci_pop_sizes)){#
    pls[i,]<-pre_loci_pop_sizes[[i]]#
  }#
  #convert pls to loci_pop_sizes format
  loci_pop_sizes<-split(pls,col(pls))
  
  
  #vectorized loci_pop_weights##################################################
  
  pre_loc_weights<- 1/pls
  loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
  loci_harm_N<-npops/colSums(pre_loc_weights)
  
  #end vectorized loci_pop_weights##############################################
  
  ###vectorize pop_alleles########################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles####################################################
  
  #vectorize allele_names#######################################################
  
  alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res<-list()
    for(i in 1:ncol(x[[1]])){
      res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    }
    return(res)
  }
  
  allele_names<-lapply(pop_alleles,alln)
  
  
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  
  #all_alleles vectorized#######################################################
  
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  
  #end all_alleles vectorized###################################################
  
  aa<-all_alleles
  aa<-lapply(aa, FUN=`list`, npops)
  afMatrix<-function(x){
    np<-x[[2]]
    z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
    rownames(z)<-x[[1]]
    return(z)
  }
  allele_freq<-lapply(aa,afMatrix)
  
  
  #combine pop_alleles
  parbind<-function(x){
    rbind(x[[1]],x[[2]])
  }
  pa1<-lapply(pop_alleles, parbind)
  #create a function to tabulate the occurance of each allele
  afTab<-function(x){
    apply(x,2,table)
  }
  actab<-lapply(pa1, afTab)
  
  afs<-function(x){
    afsint<-function(y){
      length(na.omit(y))/2
    }
    apply(x,2,afsint)
  }
  indtyppop<-lapply(pa1,afs)
  #calculate allele frequencies
  afCalcpop<-lapply(1:length(actab), function(x){
    lapply(1:length(actab[[x]]),function(y){
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
    })
  })
  #assign allele freqs to frequency matrices
  obs_count<-allele_freq
  for(i in 1:npops){
    for(j in 1:nloci){
      allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
    }
  }
  
  
  
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:npops){
    for(j in 1:nloci){
      indtyp[[j]][i]<-indtyppop[[i]][j]
    }
  }
  
  if(bootstrap==T){
    ind_vectors<-list()
    for(i in 1:npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                     ncol=(nloci+1))
    pre_data[1,]<-c("Title",rep("\t",nloci))
    for(i in 2:(nloci+1)){
      pre_data[i,1]<-loci_names[(i-1)]
    }
    pop_data<-list()
    for(i in 1:npops){
      pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                  cbind(ind_vectors[[i]],pop_list[[i]])),
                            ncol=(nloci+1))
    }
    bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
    for(i in 2:npops){
      bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
    }
    bs_data_file<-data.frame(bs_data_file)
  }
  nalleles<-vector()
  for(i in 1:nloci){
    nalleles[i]<- nrow(allele_freq[[i]])
  }
  ##############################################################################
  if(bootstrap==T){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         bs_file=bs_data_file,
         obs_allele_num=obs_count)
  } else if(bootstrap==F){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         obs_allele_num=obs_count)
  }
}
################################################################################
# readGenepop end                                                              #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# pre.div, a function to calculate Gst G'st and Djost                          #
################################################################################
pre.div<-function(x){
  data1<-x #x is a readGenepop out object
  locs<-x$locs
  ##############################################################################
  # create 'easy use' objects from data1 (readGenepop output)
  # pl = pop_list
  pl<-data1$pop_list
  # np = npops
  np<-data1$npops
  # nl = nloci
  nl<-data1$nloci
  # ps = pop sizes
  ps<-data1$pop_sizes
  # pa = pop alleles
  pa<-data1$pop_alleles
  # ant = allele names total
  ant<-data1$all_alleles
  # af = allele frequencies
  af<-data1$allele_freq
  # lnharm = locus harmonic sample size
  lnharm<-round(as.numeric(data1$loci_harm_N),4)
  # ln = locus names
  ln<-data1$loci_names
  # pn = population names
  pn<-data1$pop_names
  # ntpl = number (of individuals) typed per locus
  nt<-data1$indtyp
  ##############################################################################
  #observed heterozygosity count vectorize######################################
  
  ohcFUN<-function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data<-lapply(pa, ohcFUN)
  ohcConvert<-function(x){
    matrix(unlist(x),nrow=length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  
  
  #end observed heterozygosity count vectorize##################################
  #exhmf & exhtf vectorize######################################################
  
  square<-function(x){x^2}
  exhmf<-lapply(af, square)
  exhtf<-do.call("rbind",lapply(exhmf,function(x){
    1-colSums(x)  
  }))
  
  #end exhmf & exhtf vectorize##################################################
  #mean frequency vectorize#####################################################
  
  mf<-lapply(af,function(x){
    rowSums(x)/np  
  })
  mexhmf<-lapply(mf,square)
  ht<-sapply(mexhmf, function(x){
    1-sum(x)
  })
  ht[ht=="NaN"]<-NA
  
  #end mean frequency vectorize#################################################
  
  
  ###end locus stats legacy code
  #locus stats vectorize########################################################
  
  hs<-round(rowSums(exhtf)/np,4)
  hs_est<-round(hs*((2*lnharm)/((2*lnharm)-1)),4)
  ht_est<-round((ht + (hs_est/(2*lnharm*np))),4)
  hst<-(ht-hs)/(1-hs)
  dst<-ht-hs
  gst<-dst/ht
  djost<-((ht-hs)/(1-hs))*(np/(np-1))
  hst_est<-(ht_est-hs_est)/(1-hs_est)
  dst_est<-ht_est- hs_est
  gst_est<-(ht_est-hs_est)/ht_est
  gst_max<-((np-1)*(1-hs))/(np-1+hs)
  gst_est_max<-(((np-1)*(1-hs_est))/(np-1+hs_est))
  gst_hedrick<-gst/gst_max
  gst_est_hedrick<-gst_est/gst_est_max
  djost_est<-(np/(np-1))*((ht_est-hs_est)/(1 - hs_est))
  
  #end locus stats vectorize####################################################
  # Across all loci stats #
  ht_mean<-round(mean(na.omit(ht)),4)
  hs_mean<-round(mean(hs),4)
  gst_all<-round((ht_mean-hs_mean)/ht_mean,4)
  gst_all_max<-round(((np-1)*(1-hs_mean))/(np-1+hs_mean),4)
  gst_all_hedrick<-round(gst_all/gst_all_max,4)
  djost_all<-round(((ht_mean-hs_mean)/(1-hs_mean))*(np/(np-1)),4)
  ##############################################################################
  # Across all loci estimated stats #
  hs_est_mean<-round(mean(hs_est),4)
  ht_est_mean<-round(mean(na.omit(ht_est)),4)
  gst_est_all<-round((ht_est_mean-hs_est_mean)/ht_est_mean,4)
  gst_est_all_max<-round((((np-1)*(1-hs_est_mean))/(np-1+hs_est_mean)),4)
  gst_est_all_hedrick<-round(gst_est_all/gst_est_all_max,4)
  #djost_est_all<-round((np/(np-1))*((ht_est_mean-hs_est_mean)/
  #(1 - hs_est_mean)),4)
  djost_est_all<-round(1/((1/mean(na.omit(djost_est))+(var(na.omit(djost_est))*
    ((1/mean(na.omit(djost_est)))^3)))),4)
  ##############################################################################
  if(locs==T){  
    list(hs=hs,
         hst=hst,
         dst=dst,
         gst=gst,
         djost=djost,
         hs_est=hs_est,
         ht_est=ht_est,
         hst_est=hst_est,
         dst_est=dst_est,
         gst_est=gst_est,
         djost_est=djost_est,
         gst_max=gst_max,
         gst_est_max=gst_est_max,
         gst_hedrick=gst_hedrick,
         gst_est_hedrick=gst_est_hedrick,
         ht_mean=ht_mean,
         hs_mean=hs_mean,
         gst_all=gst_all,
         gst_all_max=gst_all_max,
         gst_all_hedrick=gst_all_hedrick,
         djost_all=djost_all,
         hs_est_mean=hs_est_mean,
         ht_est_mean=ht_est_mean,
         gst_est_all=gst_est_all,
         gst_est_all_max=gst_est_all_max,
         pop_sizes=ps,
         gst_est_all_hedrick=gst_est_all_hedrick,
         djost_est_all=djost_est_all,
         locus_names=ln,
         locus_harmonic_N=lnharm,
         npops=np,
         nloci=nl,
         pop_list=pl,
         pop_names=pn)
  } else {
    list(ht_mean=ht_mean,
         hs_mean=hs_mean,
         gst_all=gst_all,
         gst_all_max=gst_all_max,
         gst_all_hedrick=gst_all_hedrick,
         djost_all=djost_all,
         hs_est_mean=hs_est_mean,
         ht_est_mean=ht_est_mean,
         gst_est_all=gst_est_all,
         gst_est_all_max=gst_est_all_max,
         pop_sizes=ps,
         gst_est_all_hedrick=gst_est_all_hedrick,
         djost_est_all=djost_est_all,
         locus_names=ln,
         locus_harmonic_N=lnharm,
         npops=np,
         nloci=nl,
         pop_list=pl,
         pop_names=pn)
  }
}
################################################################################
# pre.div end                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# plotter, a function to create interactive plots of results from div.part     #
################################################################################

plotter<-function(x,img="1200x600"){
  image.size=img
  x=x
  spot.radius=5
  jjj<-x
  require("sendplot")
  fl_ext<-c(".tif","Dot.png","Dot.tif")
  bs_res<-list()
  lso123<-list()
  accDat<-list()
  sp.header<-list()
  pw_res<-list()
  pwso<-list()
  pw<-list()
  plot_data321<-list()
  if(jjj$plot_loci==TRUE && jjj$plot_pw==FALSE){
    bs_res<<-jjj$bs_res
    lso123<<-jjj$lso123
    accDat<<-jjj$accDat
    #Gst_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[1]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[1]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    #clean up
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
    #G'st_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[2]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[2]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
    #Djost_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[3]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[3]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
    #Fst_WC_loci
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[4]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[4]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
    }
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("bs_res", where=".GlobalEnv")==TRUE){
      rm(bs_res, pos=".GlobalEnv")
    }
    if(exists("lso123", where=".GlobalEnv")==TRUE){
      rm(lso123, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    #rm(jjj,accDat,bs_res,lso123,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==FALSE && jjj$plot_pw==TRUE){
    accDat<<-jjj$accDat
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data321<<-jjj$plot_data321
    #Gst_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[1]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[1]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",fl_ext,sep=""))
    #G'st_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[2]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[2]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
    #Djost_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[3]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[3]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",fl_ext,sep=""))
    #Fst_WC_pw
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[4]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[4]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                   fl_ext,sep=""))
    }
    
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("pw_res", where=".GlobalEnv")==TRUE){
      rm(pw_res, pos=".GlobalEnv")
    }
    if(exists("pwso", where=".GlobalEnv")==TRUE){
      rm(pwso, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    if(exists("pw", where=".GlobalEnv")==TRUE){
      rm(pw, pos=".GlobalEnv")
    }
    #rm(jjj,accDat,plot_data,pw,pw_res,pwso,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==TRUE && jjj$plot_pw==TRUE){
    bs_res<<-jjj$bs_res
    lso123<<-jjj$lso123
    accDat<<-jjj$accDat
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data321<<-jjj$plot_data321
    #Gst_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[1]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[1]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
    #G'st_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[2]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[2]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
    #Djost_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[3]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[3]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
    #Fst_WC_loci
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[4]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[4]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
    }
    #Gst_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[1]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[1]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",
                 fl_ext,sep=""))
    #G'st_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[2]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[2]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
    #Djost_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[3]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[3]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",
                 fl_ext,sep=""))
    #Fst_WC_pw
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[4]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[4]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                   fl_ext,sep=""))
    }
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("pw_res", where=".GlobalEnv")==TRUE){
      rm(pw_res, pos=".GlobalEnv")
    }
    if(exists("pwso", where=".GlobalEnv")==TRUE){
      rm(pwso, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    if(exists("pw", where=".GlobalEnv")==TRUE){
      rm(pw, pos=".GlobalEnv")
    }
    if(exists("bs_res", where=".GlobalEnv")==TRUE){
      rm(bs_res, pos=".GlobalEnv")
    }
    if(exists("lso123", where=".GlobalEnv")==TRUE){
      rm(lso123, pos=".GlobalEnv")
    }
  }
}
################################################################################
# plotter end                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# in.calc, a wrapper function for the calculation of locus informativeness     #
################################################################################
in.calc<-function(infile, outfile=NULL, gp=3, bs_locus=FALSE, bs_pairwise=FALSE, 
                  bootstraps=0, Plot=FALSE, parallel=FALSE){
  D=infile
  gp=gp
  pw=bs_pairwise
  BS=bs_locus
  NBS=bootstraps
  on=outfile
  plt=Plot
  para=parallel
  if(pw==T && NBS<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else if (BS==T && NBS<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else {
    write_res<-is.element("xlsx",installed.packages()[,1])
    if(write_res==TRUE){
      require("xlsx")
    } else {
      Warning1<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please install the package 'xlsx' if you would like your", 
                       "results written to an Excel workbook.",
                       "Alternatively, your result will automatically be written",
                       "to .txt files.",
                       "___________________________________________________________",
                       "To install 'xlsx' use:",
                       "> install.packages('xlsx', dependencies=TRUE)",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning1))
    }
    if(!is.null(on)){
      suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                             "-[diveRsity]","/",sep="")))
    }
    
    of<-paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
    # Parallel system opti
    if(para==TRUE){
      if(Sys.info()["sysname"][[1]]=="Linux"){
        para_pack_inst<-is.element(c("snow","doSNOW","foreach","iterators"),
                                   installed.packages()[,1])
      } else {
        para_pack_inst<-is.element(c("parallel","doParallel","foreach",
                                     "iterators"),installed.packages()[,1])
      }
      para_pack<-any(para_pack_inst==FALSE)
    }
    if (para == TRUE && para_pack == TRUE){
      if(Sys.info()["sysname"][[1]]=="Linux"){
        Warning3<-{paste(" "," ",
                         "[NOTE]",
                         "___________________________________________________________",
                         "Please make sure the packages 'doSNOW', 'snow', 'foreach'",
                         " and 'iterators' are installed. These are required to run",
                         " your analysis in parallel.",
                         "Your analysis will be run sequentially!",
                         "___________________________________________________________",
                         "To install these use:",
                         "> install.packages()",
                         "See:",
                         "> ?install.packages - for usage details.",
                         "___________________________________________________________",
                         sep="\n")
        }
      } else {
        Warning3<-{paste(" "," ",
                         "[NOTE]",
                         "___________________________________________________________",
                         "Please make sure the packages 'parallel', 'doParallel',",
                         "'foreach' and 'iterators' are installed. These are required",
                         " to run your analysis in parallel.",
                         "Your analysis will be run sequentially!",
                         "___________________________________________________________",
                         "To install these use:",
                         "> install.packages()",
                         "See:",
                         "> ?install.packages - for usage details.",
                         "___________________________________________________________",
                         sep="\n")
        }
      }
      cat(noquote(Warning3))
    }
    ##
    
    #source("in.bootstrap.R")
    inls2<-list(D,gp,"FALSE",0,"FALSE")
    res_out<-in.bs(inls2)[[1]]
    if(!is.null(on)){
      if(write_res==TRUE){
        write.xlsx(res_out,file=paste(of,"[In.calc].xlsx",sep=""),
                   sheetName="In_allele_stats",col.names=T,row.names=T,append=F)
      } else {
        all_out<-file(paste(of,"Allele-In[in.calc].txt",sep=""),"w")
        cat(paste(colnames(res_out),sep=""),"\n",sep="\t",file=all_out)
        for(i in 1:nrow(res_out)){
          cat(res_out[i,],"\n",sep="\t",file=all_out)
        }
        close(all_out)
      }
    }
    ######################################################################
    # overall In
    if(BS==T){
      inls1<-list(D,gp,BS,NBS,"TRUE")
      bs_sum1<-in.bs(inls1)
      if(write_res==T){
        write.xlsx(bs_sum1,file=paste(of, "[In.calc].xlsx",sep=""),
                   sheetName="Overall_Bootstrap",col.names=T,
                   row.names=T,append=T)
      } else {
        all_bs<-file(paste(of,"Overall-bootstrap[in.calc].txt",sep=""),"w")
        cat(paste(colnames(bs_sum1),sep=""),"\n",sep="\t",file=all_bs)
        for(i in 1:nrow(bs_sum1)){
          cat(bs_sum1[i,],"\n",sep="\t",file=all_bs)
        }
        close(all_bs)
      }
      loc_nms<-rownames(bs_sum1)
      if(plt==T){
        lso<-order(bs_sum1[,1],decreasing=F)
        png(filename=paste(of, on,"_In_plot.png",sep=""),width=800,height=600)
        par(mar=c(6,5,1,1))
        plot(bs_sum1[lso,1],ylim=c(0,(max(bs_sum1[,3])+0.1)),xaxt='n',
             ylab=expression('Locus '*I[n]),
             xlab="",cex.lab=1.5,cex.axis=1.3,las=1,type='n')
        points(bs_sum1[lso,1],pch=15,col='black',cex=1)
        arrows(1:nrow(bs_sum1),bs_sum1[lso,2],1:nrow(bs_sum1),bs_sum1[lso,3],
               code=3,angle=90,length=0.05,lwd=0.1)
        axis(1,at=1:nrow(bs_sum1),labels=loc_nms[lso],las=3)
        dev.off()
      }
    }
    # pairwise locus In bootstrap
    if(pw==T){
      inls<-list(D, gp, FALSE, TRUE)
      names(inls)<-c("infile","gp","bootstrap","locs")
      data<-readGenepop(inls)
      af<-data$allele_freq
      np<-data$npops
      nl<-data$nloci
      nal<-data$nalleles
      ln<-data$loci_names
      ps<-data$pop_sizes
      pl<-data$pop_list
      pwc<-combn(np,2)
      pn<-data$pop_names
      
      iv<-list()
      for(i in 1:np){
        iv[[i]]<-noquote(paste(rep(i,ps[i]),",",sep=""))
      }
      
      
      pre_data<-matrix(rep("",((nl+1)*(nl+1))),
                       ncol=(nl+1))
      pre_data[1,]<-rep("",(nl+1))
      for(i in 2:(nl+1)){
        pre_data[i,1]<-ln[(i-1)]
      }
      
      pw_data<-list()
      for (i in 1:ncol(pwc)){
        pw_data[[i]]<-data.frame(rbind(pre_data,
                                       c("POP",as.vector(rep("",nl))),
                                       cbind(iv[[pwc[1,i]]],
                                             matrix(noquote(pl[[pwc[1,i]]]),
                                                    ncol=nl)),
                                       c("POP",as.vector(rep("",nl))),
                                       cbind(iv[[pwc[2,i]]],
                                             matrix(noquote(pl[[pwc[2,i]]]),
                                                    
                                                    ncol=nl))))
      }
      pw_bs<-list()
      pw_bs_in<-list()
      pw_only<-TRUE
      pw_bs_out<-list()
      for(i in 1:ncol(pwc)){
        pw_bs_in[[i]]<-list(pw_data[[i]],gp,pw,NBS,pw_only)
      }
      if (para == TRUE && para_pack == FALSE){
        if(Sys.info()["sysname"][[1]]=="Linux"){
          library(parallel)
          core<-detectCores()
          detach(package:parallel)
          library(doSNOW)
          cl<-makeCluster(core,type="SOCK")
          registerDoSNOW(cl)
          pw_bs<-parLapply(cl,pw_bs_in,in.bs)
        } else if(Sys.info()["sysname"][[1]]=="windows"){
          library(doParallel)
          cl<-makeCluster(detectCores())
          registerDoParallel(cl)
          pw_bs<-parLapply(cl,pw_bs_in, in.bs)
        }
        stopCluster(cl)
      } else {
        pw_bs<-lapply(pw_bs_in, in.bs)
      }
      for(i in 1:ncol(pwc)){
        #  pw_bs[[i]]<-in.bs(pw_data[[i]],gp,pw,NBS)[[2]]
        pw_bs_out[[i]]<-matrix(cbind(rownames(pw_bs[[i]]),
                                     pw_bs[[i]][,1:3]),ncol=4)
      }
      pw_nms<-paste(pn[pwc[1,]],pn[pwc[2,]],sep=" vs. ")
      names(pw_bs)<-pw_nms
      hdr<-c("Loci","Actual_In","Lower_95CI","Upper_95CI")
      pw_in_bs<-matrix(rbind(hdr,c(names(pw_bs)[1],"","",""),pw_bs_out[[1]]),
                       ncol=4)
      for(j in 2:ncol(pwc)){
        pw_in_bs<-matrix(rbind(pw_in_bs,c(names(pw_bs)[j],"","",""),
                               pw_bs_out[[j]]),ncol=4)
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(pw_in_bs,file=paste(of, "[In.calc].xlsx",sep=""),
                     sheetName="Pairwise_bootstraps",col.names=F,
                     row.names=F,append=T)
        } else {
          pw_bs<-file(paste(of,"Pairwise-bootstrap[in.calc].txt",sep=""),"w")
          cat(paste(colnames(pw_in_bs),sep=""),"\n",sep="\t",file=pw_bs)
          for(i in 1:nrow(pw_in_bs)){
            cat(pw_in_bs[i,],"\n",sep="\t",file=pw_bs)
          }
          close(pw_bs)
        }
      }
    }
    
    if(BS==F && pw==F){
      list(Allele_In=res_out)
    } else if (BS==T && pw==F){
      list(Allele_In=res_out,
           l_bootstrap=bs_sum1)
    } else if (BS==F && pw==T){
      list(Allele_In=res_out,
           PW_bootstrap=pw_bs)
    } else if (BS==T && pw==T){
      list(Allele_In=res_out,
           l_bootstrap=bs_sum1,
           PW_bootstrap=pw_bs)
    }      
  }
}
################################################################################
# in.calc end                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# in.bs, a function for the bootstrap calculations of locus informativeness    #
################################################################################
in.bs<-function(x){
  D=x[[1]]
  gp=x[[2]]
  BS=x[[3]]
  NBS=x[[4]]
  pw_only=x[[5]]
  readGenepop<- function (x) {
    gp=x$gp
    infile=x$infile
    bootstrap=x$bootstrap
    locs=x$locs
    if(typeof(infile)=="list"){
      data1=infile 
    } else if (typeof(infile)=="character"){
      no_col <- max(count.fields(infile))
      suppressWarnings(flth<-length(readLines(infile)))
      data1<-matrix(NA,ncol=no_col,nrow=flth)
      datt<-list()
      for(i in 1:flth){
        datt[[i]]<-scan(infile,skip=(i-1),nlines=1,fill=T, 
                        what=list(rep("\t",no_col)),quiet = T)
        datt[[i]]<-unlist(datt[[i]])
        if(is.element(",", datt[[i]])==TRUE){
          datt[[i]]<-datt[[i]][-(which(datt[[i]] == ","))]
        }
        if(length(datt[[i]]) < no_col){
          datt[[i]]<-c(datt[[i]],rep(NA, no_col-length(datt[[i]])))
        }
        data1[i,]<-datt[[i]]
      }
      data1<-as.data.frame(data1)
    }
    data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
    raw_data<-data1
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
    pop_weights<- 1/pop_sizes
    
    n_harmonic<-npops/sum(pop_weights)
    
    N<-pop_sizes
    
    nloci<- (pop_pos[1]-2)
    if(nloci != (ncol(raw_data)-1)){
      stop("Check your input file for formatting errors!")
    }
    loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
    pop_list<-list()
    for (i in 1:npops){
      pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                     2:(nloci+1)])
    }
    
    
    
    if (gp==3) {
      plMake<-function(x){
        return(matrix(sprintf("%06g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
      }
    } else if (gp==2) {
      plMake<-function(x){
        return(matrix(sprintf("%04g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
      }
    }
    pop_list<-lapply(pop_list, plMake)
    
    
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]]=="    NA"]<-NA
    }
    
    
    if(bootstrap == T){
      bs<-function(x){
        return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
      }
      pop_list<-lapply(pop_list, bs)
    }  
    
    ###vectorize loci_pop_sizes#####################################################
    
    lps<-function(x){#
      lsp_count<-as.vector(colSums(!is.na(x)))#
      return(lsp_count)#
    }#
    pre_loci_pop_sizes<-lapply(pop_list,lps)#
    pls<-matrix(ncol=nloci,nrow=npops)#
    for(i in 1:length(pre_loci_pop_sizes)){#
      pls[i,]<-pre_loci_pop_sizes[[i]]#
    }#
    #convert pls to loci_pop_sizes format
    loci_pop_sizes<-split(pls,col(pls))
    
    
    #vectorized loci_pop_weights##################################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##############################################
    
    ###vectorize pop_alleles########################################################
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    #end vectorize pop_alleles####################################################
    
    #vectorize allele_names#######################################################
    
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    
    allele_names<-lapply(pop_alleles,alln)
    
    
    loci_combi<-allele_names[[1]]
    for(j in 1:nloci){
      for(i in 2:npops){
        loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
      }
    }
    
    #all_alleles vectorized#######################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###################################################
    
    aa<-all_alleles
    aa<-lapply(aa, FUN=`list`, npops)
    afMatrix<-function(x){
      np<-x[[2]]
      z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
      rownames(z)<-x[[1]]
      return(z)
    }
    allele_freq<-lapply(aa,afMatrix)
    
    
    #combine pop_alleles
    parbind<-function(x){
      rbind(x[[1]],x[[2]])
    }
    pa1<-lapply(pop_alleles, parbind)
    #create a function to tabulate the occurance of each allele
    afTab<-function(x){
      apply(x,2,table)
    }
    actab<-lapply(pa1, afTab)
    
    afs<-function(x){
      afsint<-function(y){
        length(na.omit(y))/2
      }
      apply(x,2,afsint)
    }
    indtyppop<-lapply(pa1,afs)
    #calculate allele frequencies
    afCalcpop<-lapply(1:length(actab), function(x){
      lapply(1:length(actab[[x]]),function(y){
        actab[[x]][[y]]/(indtyppop[[x]][y]*2)
      })
    })
    #assign allele freqs to frequency matrices
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      }
    }
    
    
    
    indtyp<-list()
    for(i in 1:nloci){
      indtyp[[i]]<-vector()
    }
    for(i in 1:npops){
      for(j in 1:nloci){
        indtyp[[j]][i]<-indtyppop[[i]][j]
      }
    }
    
    if(bootstrap==T){
      ind_vectors<-list()
      for(i in 1:npops){
        ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
      }
      
      
      pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                       ncol=(nloci+1))
      pre_data[1,]<-c("Title",rep("\t",nloci))
      for(i in 2:(nloci+1)){
        pre_data[i,1]<-loci_names[(i-1)]
      }
      pop_data<-list()
      for(i in 1:npops){
        pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                    cbind(ind_vectors[[i]],pop_list[[i]])),
                              ncol=(nloci+1))
      }
      bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
      for(i in 2:npops){
        bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
      }
      bs_data_file<-data.frame(bs_data_file)
    }
    nalleles<-vector()
    for(i in 1:nloci){
      nalleles[i]<- nrow(allele_freq[[i]])
    }
    ##############################################################################
    if(bootstrap==T){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           locs=locs,
           bs_file=bs_data_file)
    } else if(bootstrap==F){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           locs=locs)
    }
  }
  inls<-list(D, gp, FALSE, TRUE)
  names(inls)<-c("infile","gp","bootstrap","locs")
  data<-readGenepop(inls)
  af<-data$allele_freq
  np<-data$npops
  nl<-data$nloci
  nal<-data$nalleles
  ln<-data$loci_names
  ps<-data$pop_sizes
  pl<-data$pop_list
  ## Calc P[i]
  p<-list()
  for (i in 1:nl){
    p[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      p[[i]][j]<- sum(af[[i]][j,])/np
    }
  }
  exp1<-list()
  for (i in 1:nl){
    exp1[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      exp1[[i]][j]<- (-p[[i]][j]*log(p[[i]][j]))
    }
  }
  exp2_sep<-list()
  for(i in 1:nl){
    exp2_sep[[i]]<-matrix(rep(0,np*nal[i]))
    dim(exp2_sep[[i]])<-c(nal[i],np)
  }
  for(i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      for (z in 1:np){
        exp2_sep[[i]][j,z]<-(af[[i]][j,z]/np)*
          log(af[[i]][j,z])
      }
    }
  }
  ## Replace NaN's with 0.000
  for (i in 1:nl){
    exp2_sep[[i]][exp2_sep[[i]]=="NaN"]<-0
  }
  exp2<-list()
  for (i in 1:nl){
    exp2[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      exp2[[i]][j]<-sum(exp2_sep[[i]][j,])
    }
  }
  In<-list()
  for (i in 1:nl){
    In[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      In[[i]][j]<-sum(exp1[[i]][j]+exp2[[i]][j])
    }
  }
  In_sum<-vector()
  for (i in 1:nl){
    In_sum[i]<-sum(In[[i]])
  }
  results_out<-matrix(rep(NA,(nl*(max(nal)+1))),
                      nrow=nl,ncol=(max(nal)+1))
  for (i in 1:nl){
    results_out[i,1:nal[i]]<- round(In[[i]],4)
  }
  results_out[,(max(nal)+1)]<-round(In_sum,4)
  ##############################################################################
  if(BS==T){
    bs_sum<-matrix(rep(0,(nl*NBS)),ncol=nl)
    colnames(bs_sum)<-ln
    inls_bs<-list(D,gp,TRUE,TRUE)
    names(inls_bs)<-c("infile","gp","bootstrap","locs")
    for(w in 1:NBS){
      data_bs<-readGenepop(inls_bs)
      af_bs<-data_bs$allele_freq
      np_bs<-data_bs$npops
      nl_bs<-data_bs$nloci
      nal_bs<-data_bs$nalleles
      ln_bs<-data_bs$loci_names
      ## Calc P[i]
      p_bs<-list()
      for (i in 1:nl_bs){
        p_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          p_bs[[i]][j]<- sum(af_bs[[i]][j,])/np_bs
        }
      }
      exp1_bs<-list()
      for (i in 1:nl_bs){
        exp1_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          exp1_bs[[i]][j]<- (-p_bs[[i]][j]*log(p_bs[[i]][j]))
        }
      }
      exp2_sep_bs<-list()
      for(i in 1:nl_bs){
        exp2_sep_bs[[i]]<-matrix(rep(0,np_bs*nal_bs[i]))
        dim(exp2_sep_bs[[i]])<-c(nal_bs[i],np_bs)
      }
      for(i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          for (z in 1:np_bs){
            exp2_sep_bs[[i]][j,z]<-(af_bs[[i]][j,z]/np_bs)*
              log(af_bs[[i]][j,z])
          }
        }
      }
      ## Replace NaN's with 0.000
      for (i in 1:nl_bs){
        exp2_sep_bs[[i]][exp2_sep_bs[[i]]=="NaN"]<-0
      }
      exp2_bs<-list()
      for (i in 1:nl_bs){
        exp2_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          exp2_bs[[i]][j]<-sum(exp2_sep_bs[[i]][j,])
        }
      }
      In_bs<-list()
      for (i in 1:nl_bs){
        In_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          In_bs[[i]][j]<-sum(exp1_bs[[i]][j]+exp2_bs[[i]][j])
        }
      }
      In_sum_bs<-vector()
      for (i in 1:nl_bs){
        In_sum_bs[i]<-sum(In_bs[[i]])
      }
      results_out_bs<-matrix(rep(NA,(nl_bs*(max(nal_bs)+1))),
                             nrow=nl_bs,ncol=(max(nal_bs)+1))
      for (i in 1:nl_bs){
        results_out_bs[i,1:nal_bs[i]]<- In_bs[[i]]
      }
      results_out_bs[,(max(nal_bs)+1)]<-In_sum_bs
      
      bs_sum[w,]<-results_out_bs[,(max(nal_bs)+1)]
    }
    in_bs_out<-matrix(rep(0,(nl_bs*3)),ncol=3)
    colnames(in_bs_out)<-c("In","Lower_95CI","Upper_95CI")
    rownames(in_bs_out)<-ln_bs
    in_bs_out[,1]<-round(results_out[,(max(nal)+1)],4)
    for(i in 1:nl){
      in_bs_out[i,2]<-round(results_out[i,(max(nal)+1)]-(1.96*(sd(bs_sum[,i]))),
                            4)
      in_bs_out[i,3]<-round(results_out[i,(max(nal)+1)]+(1.96*(sd(bs_sum[,i]))),
                            4)
    }
  }
  colnames(results_out)<-c(paste("Allele.",1:max(nal),sep=""),"Sum")
  rownames(results_out)<-ln
  results_out[is.na(results_out)]<-""
  if(BS==T && pw_only==F){
    list(In_alleles=results_out,
         in_bs_out=in_bs_out)
  }else if(BS==F){
    list(In_alleles=results_out)
  }else if(BS==T && pw_only==T){
    return(in_bs_out)
  }
}
################################################################################
# in.bs end                                                                    #
################################################################################
################################################################################
# div.part end                                                                 #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# readGenepop.user, a usable function for basic population parameters          #
################################################################################
readGenepop.user<- function (infile=NULL, gp=3, bootstrap=FALSE) {
  gp=gp
  infile=infile
  bootstrap=bootstrap
  #locs=locs
  if(typeof(infile)=="list"){
    data1=infile 
  } else if (typeof(infile)=="character"){
    no_col <- max(count.fields(infile))
    suppressWarnings(flth<-length(readLines(infile)))
    data1<-matrix(NA,ncol=no_col,nrow=flth)
    datt<-list()
    for(i in 1:flth){
      datt[[i]]<-scan(infile,skip=(i-1),nlines=1,fill=T, 
                      what=list(rep("\t",no_col)),quiet = T)
      datt[[i]]<-unlist(datt[[i]])
      if(is.element(",", datt[[i]])==TRUE){
        datt[[i]]<-datt[[i]][-(which(datt[[i]] == ","))]
      }
      if(length(datt[[i]]) < no_col){
        datt[[i]]<-c(datt[[i]],rep(NA, no_col-length(datt[[i]])))
      }
      data1[i,]<-datt[[i]]
    }
    data1<-as.data.frame(data1)
  }
  data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
  raw_data<-data1
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
  pop_weights<- 1/pop_sizes
  
  n_harmonic<-npops/sum(pop_weights)
  
  N<-pop_sizes
  
  nloci<- (pop_pos[1]-2)
  if(nloci != (ncol(raw_data)-1)){
    stop("Check your input file for formatting errors!")
  }
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                   2:(nloci+1)])
  }
  
  
  
  if (gp==3) {
    plMake<-function(x){
      return(matrix(sprintf("%06g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
    }
  } else if (gp==2) {
    plMake<-function(x){
      return(matrix(sprintf("%04g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x)))
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  
  
  for(i in 1:npops){
    pop_list[[i]][pop_list[[i]]=="    NA"]<-NA
  }
  
  
  if(bootstrap == T){
    bs<-function(x){
      return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
    }
    pop_list<-lapply(pop_list, bs)
  }  
  
  ###vectorize loci_pop_sizes#####################################################
  
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
  pre_loci_pop_sizes<-lapply(pop_list,lps)#
  pls<-matrix(ncol=nloci,nrow=npops)#
  for(i in 1:length(pre_loci_pop_sizes)){#
    pls[i,]<-pre_loci_pop_sizes[[i]]#
  }#
  #convert pls to loci_pop_sizes format
  loci_pop_sizes<-split(pls,col(pls))
  
  
  #vectorized loci_pop_weights##################################################
  
  pre_loc_weights<- 1/pls
  loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
  loci_harm_N<-npops/colSums(pre_loc_weights)
  
  #end vectorized loci_pop_weights##############################################
  
  ###vectorize pop_alleles########################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles####################################################
  
  #vectorize allele_names#######################################################
  
  alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res<-list()
    for(i in 1:ncol(x[[1]])){
      res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    }
    return(res)
  }
  
  allele_names<-lapply(pop_alleles,alln)
  
  
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  
  #all_alleles vectorized#######################################################
  
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  
  #end all_alleles vectorized###################################################
  
  aa<-all_alleles
  aa<-lapply(aa, FUN=`list`, npops)
  afMatrix<-function(x){
    np<-x[[2]]
    z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
    rownames(z)<-x[[1]]
    return(z)
  }
  allele_freq<-lapply(aa,afMatrix)
  
  
  #combine pop_alleles
  parbind<-function(x){
    rbind(x[[1]],x[[2]])
  }
  pa1<-lapply(pop_alleles, parbind)
  #create a function to tabulate the occurance of each allele
  afTab<-function(x){
    apply(x,2,table)
  }
  actab<-lapply(pa1, afTab)
  
  afs<-function(x){
    afsint<-function(y){
      length(na.omit(y))/2
    }
    apply(x,2,afsint)
  }
  indtyppop<-lapply(pa1,afs)
  #calculate allele frequencies
  afCalcpop<-lapply(1:length(actab), function(x){
    lapply(1:length(actab[[x]]),function(y){
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
    })
  })
  #assign allele freqs to frequency matrices
  obs_count<-allele_freq
  for(i in 1:npops){
    for(j in 1:nloci){
      allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
    }
  }
  
  
  
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:npops){
    for(j in 1:nloci){
      indtyp[[j]][i]<-indtyppop[[i]][j]
    }
  }
  
  if(bootstrap==T){
    ind_vectors<-list()
    for(i in 1:npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                     ncol=(nloci+1))
    pre_data[1,]<-c("Title",rep("\t",nloci))
    for(i in 2:(nloci+1)){
      pre_data[i,1]<-loci_names[(i-1)]
    }
    pop_data<-list()
    for(i in 1:npops){
      pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                  cbind(ind_vectors[[i]],pop_list[[i]])),
                            ncol=(nloci+1))
    }
    bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
    for(i in 2:npops){
      bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
    }
    bs_data_file<-data.frame(bs_data_file)
  }
  nalleles<-vector()
  for(i in 1:nloci){
    nalleles[i]<- nrow(allele_freq[[i]])
  }
  ##############################################################################
  if(bootstrap==T){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         #locs=locs,
         bs_file=bs_data_file,
         obs_allele_num=obs_count)
  } else if(bootstrap==F){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         #locs=locs,
         obs_allele_num=obs_count)
  }
}
################################################################################
# readGenepop.user end                                                         #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# pre.divLowMemory, a low memory consumption function for locus bootstrapping  #
################################################################################
pre.divLowMemory<-function(y){
  y<-y
  ##############################################################################
  readGenepopX<- function (x) {
    gp=x$gp
    infile=x$infile
    bootstrap=x$bootstrap
    locs=x$locs
    fst=x$fst
    if(typeof(infile)=="list"){
      data1=infile 
    } else if (typeof(infile)=="character"){
      no_col <- max(count.fields(infile))
      suppressWarnings(flth<-length(readLines(infile)))
      data1<-matrix(NA,ncol=no_col,nrow=flth)
      datt<-list()
      for(i in 1:flth){
        datt[[i]]<-scan(infile,skip=(i-1),nlines=1,fill=T, 
                        what=list(rep("\t",no_col)),quiet = T)
        datt[[i]]<-unlist(datt[[i]])
        if(is.element(",", datt[[i]])==TRUE){
          datt[[i]]<-datt[[i]][-(which(datt[[i]] == ","))]
        }
        if(length(datt[[i]]) < no_col){
          datt[[i]]<-c(datt[[i]],rep(NA, no_col-length(datt[[i]])))
        }
        data1[i,]<-datt[[i]]
      }
      data1<-as.data.frame(data1)
    }
    data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
    raw_data<-data1
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,10)
    pop_weights<- 1/pop_sizes
    
    n_harmonic<-npops/sum(pop_weights)
    
    N<-pop_sizes
    
    nloci<- (pop_pos[1]-2)
    if(nloci != (ncol(raw_data)-1)){
      stop("Check your input file for formatting errors!")
    }
    loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
    pop_list<-list()
    for (i in 1:npops){
      pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                     2:(nloci+1)])
    }
    
    
    
    if (gp==3) {
      plMake<-function(x){
        suppressWarnings(return(matrix(sprintf("%06g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x))))
      }
    } else if (gp==2) {
      plMake<-function(x){
        suppressWarnings(return(matrix(sprintf("%04g",as.numeric(x)),nrow=nrow(x),ncol=ncol(x))))
      }
    }
    pop_list<-lapply(pop_list, plMake)
    
    
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]]=="    NA"]<-NA
    }
    
    
    if(bootstrap == T){
      bs<-function(x){
        return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
      }
      pop_list<-lapply(pop_list, bs)
    }  
    
    ###vectorize loci_pop_sizes#####################################################
    
    lps<-function(x){#
      lsp_count<-as.vector(colSums(!is.na(x)))#
      return(lsp_count)#
    }#
    pre_loci_pop_sizes<-lapply(pop_list,lps)#
    pls<-matrix(ncol=nloci,nrow=npops)#
    for(i in 1:length(pre_loci_pop_sizes)){#
      pls[i,]<-pre_loci_pop_sizes[[i]]#
    }#
    #convert pls to loci_pop_sizes format
    loci_pop_sizes<-split(pls,col(pls))
    
    
    #vectorized loci_pop_weights##################################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##############################################
    
    ###vectorize pop_alleles########################################################
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    #end vectorize pop_alleles####################################################
    
    #vectorize allele_names#######################################################
    
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    
    allele_names<-lapply(pop_alleles,alln)
    
    
    loci_combi<-allele_names[[1]]
    for(j in 1:nloci){
      for(i in 2:npops){
        loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
      }
    }
    
    #all_alleles vectorized#######################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###################################################
    
    aa<-all_alleles
    aa<-lapply(aa, FUN=`list`, npops)
    afMatrix<-function(x){
      np<-x[[2]]
      z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
      rownames(z)<-x[[1]]
      return(z)
    }
    allele_freq<-lapply(aa,afMatrix)
    
    
    #combine pop_alleles
    parbind<-function(x){
      rbind(x[[1]],x[[2]])
    }
    pa1<-lapply(pop_alleles, parbind)
    #create a function to tabulate the occurance of each allele
    afTab<-function(x){
      apply(x,2,table)
    }
    actab<-lapply(pa1, afTab)
    
    afs<-function(x){
      afsint<-function(y){
        length(na.omit(y))/2
      }
      apply(x,2,afsint)
    }
    indtyppop<-lapply(pa1,afs)
    #calculate allele frequencies
    afCalcpop<-lapply(1:length(actab), function(x){
      lapply(1:length(actab[[x]]),function(y){
        actab[[x]][[y]]/(indtyppop[[x]][y]*2)
      })
    })
    #assign allele freqs to frequency matrices
    obs_count<-allele_freq
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][allele_names[[i]][[j]],i]<-afCalcpop[[i]][[j]]
        obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
      }
    }
    
    
    
    indtyp<-list()
    for(i in 1:nloci){
      indtyp[[i]]<-vector()
    }
    for(i in 1:npops){
      for(j in 1:nloci){
        indtyp[[j]][i]<-indtyppop[[i]][j]
      }
    }
    
    if(bootstrap==T){
      ind_vectors<-list()
      for(i in 1:npops){
        ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
      }
      
      
      pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                       ncol=(nloci+1))
      pre_data[1,]<-c("Title",rep("\t",nloci))
      for(i in 2:(nloci+1)){
        pre_data[i,1]<-loci_names[(i-1)]
      }
      pop_data<-list()
      for(i in 1:npops){
        pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                    cbind(ind_vectors[[i]],pop_list[[i]])),
                              ncol=(nloci+1))
      }
      bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
      for(i in 2:npops){
        bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
      }
      bs_data_file<-data.frame(bs_data_file)
    }
    nalleles<-vector()
    for(i in 1:nloci){
      nalleles[i]<- nrow(allele_freq[[i]])
    }
    ##############################################################################
    if(bootstrap==T){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           locs=locs,
           #bs_file=bs_data_file,
           obs_allele_num=obs_count,
           fst=fst,
           gp=gp)
    } else if(bootstrap==F){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           locs=locs,
           obs_allele_num=obs_count,
           fst=fst,
           gp=gp)
    }
  }
  x<-readGenepopX(y)
  data1<-x #x is a readGenepop out object
  locs<-x$locs
  fst=x$fst
  ##############################################################################
  #fstWC define
  if(fst==TRUE){
    fstWC<-function(x){
      badData<-vector()
      for(i in 1:x$nloci){
        badData[i]<-is.element(0, x$indtyp[[i]])
      }
      if(length(badData)>0){
        nl<-x$nloci-length(badData)
      } else {
        nl<-x$nloci
      }
      gdData<-which(badData==FALSE)
      badData<-which(badData==TRUE)
      all_genot<-x$pop_list[[1]][,gdData]
      for(i in 2:x$npops){
        all_genot<-rbind(all_genot,x$pop_list[[i]][,gdData])
      }
      genot<-apply(all_genot,2,unique)
      genot<-lapply(genot, function(x){
        if (sum(is.na(x))>0){
          y<-which(is.na(x)==TRUE)
          x_new<-x[-y]
          return(x_new)
        } else {
          return(x)
        }
      })
      #count genotypes
      
      genoCount<-list()
      for(i in 1:ncol(all_genot)){
        genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
        for(j in 1:length(genot[[i]])){
          genoCount[[i]][,j]<-length(which(all_genot[,i]==genot[[i]][j]))
        }
        if (x$gp==3){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                          substr(genot[[i]],4,6),sep="")
        } else if (x$gp==2){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                          substr(genot[[i]],3,4),sep="")
        }
      }
      h_sum<-list()
      for(i in 1:ncol(all_genot)){
        h_sum[[i]]<-vector()
        cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
        for(j in 1:length(x$all_alleles[[gdData[i]]])){
          het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
          het_id2<-lapply(het_id1, sum)
          het_id2<-as.vector(het_id2)
          het_id3<-which(het_id2==1)
          h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
        }
      }
      indtyp_tot<-lapply(x$indtyp, sum)
      kk_hsum<-list()
      for(i in 1:ncol(all_genot)){
        kk_hsum[[i]]<-list(h_sum[[i]],indtyp_tot[[gdData[i]]])
      }
      kk_hbar<-lapply(kk_hsum, function(x){
        return(x[[1]]/x[[2]])
      })
      pdat<-list()
      for(i in 1:ncol(all_genot)){
        pdat[[i]]<-list(x$allele_freq[[gdData[i]]],x$indtyp[[gdData[i]]])
      }
      kk_p<-lapply(pdat, function(x){
        if(is.null(x[[1]])==FALSE){
          apply(x[[1]], 1, function(y){
            y*(2*x[[2]])
          })
        }
      })
      res<-matrix(0,(x$nloci+1),3)
      colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
      rownames(res)<-c(x$loci_names, "All")
      A<-vector()
      a<-vector()
      b<-vector()
      c<-vector()
      for(i in 1:ncol(all_genot)){
        kknbar<-indtyp_tot[[gdData[i]]]/x$npops
        kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
               indtyp_tot[[gdData[i]]])/(x$npops-1)
        kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
        kkptild[kkptild=="NaN"]<-NA
        kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
        kks2<-colSums(x$indtyp[[gdData[i]]]*
          (kkptild-rep(kkpbar,each = x$npops))^2)/((x$npops-1)*kknbar)
        kkA<-kkpbar*(1-kkpbar)-(x$npops-1)*kks2/x$npops
        kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
        kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
        kkc<-kk_hbar[[i]]/2
        A[i]<-sum(kkA)
        a[i]<-sum(kka)
        b[i]<-sum(kkb)
        c[i]<-sum(kkc)
        res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
        res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
        res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
      }
      res[res=="NaN"]<-NA
      res[res==0.000]<-NA
      sumA<-sum(na.omit(A))
      suma<-sum(na.omit(a))
      sumb<-sum(na.omit(b))
      sumc<-sum(na.omit(c))
      res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
      res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
      res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
      #res[is.na(res)]<-NaN
      list(Fstats=res,
           multiLoc<-res[(x$nloci+1),])
    }
    if(x$locs==TRUE){
      fstats<-fstWC(x)[[1]]
    }else if (x$locs==FALSE){
      fstats<-fstWC(x)[[2]]
    }        
  }
  ##############################################################################
  # create 'easy use' objects from data1 (readGenepop output)
  # pl = pop_list
  pl<-data1$pop_list
  # np = npops
  np<-data1$npops
  # nl = nloci
  nl<-data1$nloci
  # ps = pop sizes
  ps<-data1$pop_sizes
  # pa = pop alleles
  pa<-data1$pop_alleles
  # ant = allele names total
  ant<-data1$all_alleles
  # af = allele frequencies
  af<-data1$allele_freq
  # lnharm = locus harmonic sample size
  lnharm<-round(as.numeric(data1$loci_harm_N),4)
  # ln = locus names
  ln<-data1$loci_names
  # pn = population names
  pn<-data1$pop_names
  # ntpl = number (of individuals) typed per locus
  nt<-data1$indtyp
  ##############################################################################
  #observed heterozygosity count vectorize######################################
  
  ohcFUN<-function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data<-lapply(pa, ohcFUN)
  ohcConvert<-function(x){
    matrix(unlist(x),nrow=length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  
  
  #end observed heterozygosity count vectorize##################################
  #exhmf & exhtf vectorize######################################################
  
  square<-function(x){x^2}
  exhmf<-lapply(af, square)
  exhtf<-do.call("rbind",lapply(exhmf,function(x){
    1-colSums(x)  
  }))
  
  #end exhmf & exhtf vectorize##################################################
  #mean frequency vectorize#####################################################
  
  mf<-lapply(af,function(x){
    rowSums(x)/np  
  })
  mexhmf<-lapply(mf,square)
  ht<-sapply(mexhmf, function(x){
    1-sum(x)
  })
  ht[ht=="NaN"]<-NA
  
  #end mean frequency vectorize#################################################
  
  
  ###end locus stats legacy code
  #locus stats vectorize########################################################
  
  hs<-round(rowSums(exhtf)/np,4)
  hs_est<-round(hs*((2*lnharm)/((2*lnharm)-1)),4)
  ht_est<-round((ht + (hs_est/(2*lnharm*np))),4)
  ht_est[ht_est=="NaN"]<-NA
  hst<-(ht-hs)/(1-hs)
  dst<-ht-hs
  gst<-dst/ht
  djost<-((ht-hs)/(1-hs))*(np/(np-1))
  hst_est<-(ht_est-hs_est)/(1-hs_est)
  dst_est<-ht_est- hs_est
  gst_est<-(ht_est-hs_est)/ht_est
  gst_max<-((np-1)*(1-hs))/(np-1+hs)
  gst_est_max<-(((np-1)*(1-hs_est))/(np-1+hs_est))
  gst_hedrick<-gst/gst_max
  gst_est_hedrick<-gst_est/gst_est_max
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est<-(np/(np-1))*((ht_est-hs_est)/(1 - hs_est))
  
  #end locus stats vectorize####################################################
  # Across all loci stats #
  ht_mean<-round(mean(na.omit(ht)),4)
  hs_mean<-round(mean(hs),4)
  gst_all<-round((ht_mean-hs_mean)/ht_mean,4)
  gst_all_max<-round(((np-1)*(1-hs_mean))/(np-1+hs_mean),4)
  gst_all_hedrick<-round(gst_all/gst_all_max,4)
  djost_all<-round(((ht_mean-hs_mean)/(1-hs_mean))*(np/(np-1)),4)
  ##############################################################################
  # Across all loci estimated stats #
  hs_est_mean<-round(mean(hs_est),4)
  ht_est_mean<-round(mean(na.omit(ht_est)),4)
  gst_est_all<-round((ht_est_mean-hs_est_mean)/ht_est_mean,4)
  gst_est_all_max<-round((((np-1)*(1-hs_est_mean))/(np-1+hs_est_mean)),4)
  gst_est_all_hedrick<-round(gst_est_all/gst_est_all_max,4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  #djost_est_all<-round((np/(np-1))*((ht_est_mean-hs_est_mean)/
  #(1 - hs_est_mean)),4)
  djost_est_all<-round(1/((1/mean(na.omit(djost_est))+(var(na.omit(djost_est))*
    ((1/mean(na.omit(djost_est)))^3)))),4)
  djost_est[djost_est==0]<-NaN
  djost[djost==0]<-NaN
  ##############################################################################
  if(fst==TRUE){
    if(locs==T){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else {
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    }
  } else if(fst==FALSE){
    if(locs==T){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    } else {
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    }
  }
}
################################################################################
# end pre.divLowMemory                                                         #
################################################################################

#
#
#
#
#
#
#
#
#
################################################################################
# corPlot, plot the relationship between div.part stats and number of alleles  #
################################################################################
corPlot<-function(x,y){
  x=x
  y=y
  par(mfrow=c(2,2))
  par(mar=c(4,5,2,2))
  #Fst
  plot(y[[2]][1:(nrow(y[[2]])-1),8]~x[[16]],
       pch=16,xlab=expression(N[alleles]),ylab=expression(hat(theta)),
       ylim=c(0,1),las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),8]~x[[16]]),col="red",lwd=2)
  cor1<-cor.test(y[[2]][1:(nrow(y[[2]])-1),8],x[[16]])
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor1$estimate[[1]],3),sep=""),cex=2)
  #gst
  plot(y[[2]][1:(nrow(y[[2]])-1),4]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression(G[st]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),4]~x[[16]]),col="red",lwd=2)
  cor2<-cor.test(y[[2]][1:(nrow(y[[2]])-1),4],x[[16]])
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor2$estimate[[1]],3),sep=""),cex=2)
  #g'st
  plot(y[[2]][1:(nrow(y[[2]])-1),5]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression("G'"[st]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),5]~x[[16]]),col="red",lwd=2)
  cor3<-cor.test(y[[2]][1:(nrow(y[[2]])-1),5],x[[16]])
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor3$estimate[[1]],3),sep=""),cex=2)
  #D
  plot(y[[2]][1:(nrow(y[[2]])-1),6]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression(D[est]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),6]~x[[16]]),col="red",lwd=2)
  cor4<-cor.test(y[[2]][1:(nrow(y[[2]])-1),6],x[[16]])
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor4$estimate[[1]],3),sep=""),cex=2)
}
################################################################################
# end corPlot                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# difPlot, plot all pairwise population pairs                                  #
################################################################################
difPlot <- function (x, outfile= NULL, interactive = FALSE) {
  x=x
  on=outfile
  inta<-interactive
  #output from div.part
  require(plotrix)
  
  if(is.null(on) == TRUE && inta == TRUE){
    noter<-"To plot interactively, you need to specifiy a output folder"
    cat(noter)
  } else {
    of=paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
  }
  
  if(!exists("inta",-1)){
    inta<-FALSE
  }
  if(inta == TRUE) {
    sp.header<-list()
    colleer<-list()
    colleer<<-colorRampPalette(c("blue","white"))
    require(sendplot)
    direct<-of
    pwc<-combn(ncol(x[[3]][[1]]),2)
    pwNames<-paste(colnames(x[[3]][[1]])[pwc[1,]],
                   colnames(x[[3]][[1]])[pwc[2,]],
                   sep='vs ')
                   
    gst_lab <- as.vector(x[[3]][[4]])
    gst_lab <- na.omit(gst_lab)
    collab111<-list()
    #
    if(length(x[[3]]) > 6){
      fst_lab <- as.vector(x[[3]][[8]])
      fst_lab<-na.omit(fst_lab)
    }
    #
    gpst_lab <- as.vector(x[[3]][[5]])
    gpst_lab<-na.omit(gpst_lab)
    #
    Dest_lab <- as.vector(x[[3]][[6]])
    Dest_lab<-na.omit(Dest_lab)
    #
    
    fl_ext<-c(".tif","Dot.png","Dot.tif")
    if (length(x[[3]]) > 6){
      xy.labels <-  data.frame(pops = pwNames,
                               Nei_Gst = gst_lab,
                               Weir_Theta = fst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    } else {
      xy.labels <-  data.frame(pop = pwNames,
                               Nei_Gst = gst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    }
    #Nei Gst
    abx<-list()
    abx<<-x
    collab111 <<- c(round(min(gst_lab),3),
                    round(mean(gst_lab),3),
                    round(max(gst_lab),3))
    
    plot.call <- "image(1:nrow(abx[[3]][[4]]),1:nrow(abx[[3]][[4]]),
abx[[3]][[4]],ylab='',xlab='',main='Pairwise Gst',xaxt='n',yaxt='n',
    col = colleer(50),las=1,cex.main=3)"
  ##
    plot.extras <- "color.legend(nrow(abx[[3]][[4]])/5,
                                  nrow(abx[[3]][[4]])/3,
    nrow(abx[[3]][[4]])/4,
    nrow(abx[[3]][[4]])/1.2,
    collab111,
    rect.col=colleer(50),
    gradient='y',
    cex=3)"
  ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="Gst_matrix",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #clean up
    unlink(paste(of,"Gst_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #Fst
    if(length(x[[3]]) > 6){
      collab111 <<- c(round(min(fst_lab),3),
                      round(mean(fst_lab),3),
                      round(max(fst_lab),3))
      plot.call <- "image(1:nrow(abx[[3]][[8]]),1:nrow(abx[[3]][[8]]),
      abx[[3]][[8]],ylab = '',xlab = '',xaxt = 'n',yaxt = 'n',
      main = 'Pairwise Fst',col = colleer(50),las = 1,cex.main = 3)"
      ##
      plot.extras <- "color.legend(nrow(abx[[3]][[8]])/5,
                                    nrow(abx[[3]][[8]])/3,
      nrow(abx[[3]][[8]])/4,
      nrow(abx[[3]][[8]])/1.2,
      collab111,
      rect.col=colleer(50),
      gradient='y',
      cex=3)"
      #
      suppressWarnings(imagesend(plot.call=plot.call,
                                 x.pos=pwc[2,],
                                 y.pos=pwc[1,],
                                 xy.type="points",
                                 xy.labels = xy.labels,
                                 plot.extras=plot.extras,
                                 fname.root="Fst_matrix",
                                 dir=of,
                                 image.size="1050x800",
                                 font.size=18,
                                 spot.radius = 10,
                                 font.type ="Arial",
                                 window.size="1100x800"))
      #clean up
      unlink(paste(of,"Fst_matrix",fl_ext,sep=""))
    }
    #
    #
    #
    #
    #
    #G'st
    collab111 <<- c(round(min(gpst_lab),3),
                    round(mean(gpst_lab),3),
                    round(max(gpst_lab),3))
    plot.call <- "image(1:nrow(abx[[3]][[5]]),1:nrow(abx[[3]][[5]]),
    abx[[3]][[5]],ylab='',xlab='',xaxt='n',yaxt='n',
    main='Pairwise Gst (Hedrick)',col = colleer(50),las=1,cex.main=3)"
    
    plot.extras <- "color.legend(nrow(abx[[3]][[5]])/5,nrow(abx[[3]][[5]])/3,
    nrow(abx[[3]][[5]])/4,nrow(abx[[3]][[5]])/1.2,collab111,
    rect.col=colleer(50),gradient='y',cex=3)"
    ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="G_prime_st_matrix",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #clean up
    unlink(paste(of,"G_prime_st_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #
    #Dest
    collab111 <<- c(round(min(Dest_lab),3),
                    round(mean(Dest_lab),3),
                    round(max(Dest_lab),3))
    plot.call <- "image(1:nrow(abx[[3]][[6]]),1:nrow(abx[[3]][[6]]),
    abx[[3]][[6]],ylab='',xlab='',xaxt='n',yaxt='n',main='Pairwise D (Jost)',
    col = colleer(50),las=1,cex.main=3)"
     plot.extras <- "color.legend(nrow(abx[[3]][[6]])/5,nrow(abx[[3]][[6]])/3,
     nrow(abx[[3]][[6]])/4,nrow(abx[[3]][[6]])/1.2,collab111,
     rect.col=colleer(50),gradient='y',cex=3)"
    ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="D_matrix_",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #lean up
    
    unlink(paste(of,"D_matrix_",fl_ext,sep=""))
    
  } else {
    
    #
    if(length(x[[3]]) > 6){
      par(mfrow=c(2,2))
    } else {
      par(mfrow=c(3,1))
    }
    colleer<-colorRampPalette(c("blue","white"))
    cols<-colleer(50)
    #Gst
    image(1:nrow(x[[3]][[4]]),
          1:nrow(x[[3]][[4]]),
          x[[3]][[4]],
          ylab="population",
          xlab="population",
          main="Pairwise Gst",
          col=cols,
          las=1)
    gst<-as.vector(x[[3]][[4]])
    gst<-as.vector(na.omit(gst))
    collab111<-c(round(min(gst),3),
                 round(mean(gst),3),
                 round(max(gst),3))
    
    color.legend(nrow(x[[3]][[4]])/5,
                 nrow(x[[3]][[4]])/3,
                 nrow(x[[3]][[4]])/4,
                 nrow(x[[3]][[4]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    if(length(x[[3]]) > 6){
      #Fst
      image(1:nrow(x[[3]][[8]]),
            1:nrow(x[[3]][[8]]),
            x[[3]][[8]],
            ylab="population",
            xlab="population",
            main="Pairwise Theta",
            col = cols,
            las=1)
      fst<-as.vector(x[[3]][[8]])
      fst<-as.vector(na.omit(fst))
      collab111<-c(round(min(fst),3),round(mean(fst),3),round(max(fst),3))
      
      color.legend(nrow(x[[3]][[8]])/5,
                   nrow(x[[3]][[8]])/3,
                   nrow(x[[3]][[8]])/4,
                   nrow(x[[3]][[8]])/1.2,
                   collab111,
                   cols,
                   gradient="y")
    }
    #Hedrick's Gst
    image(1:nrow(x[[3]][[5]]),
          1:nrow(x[[3]][[5]]),
          x[[3]][[5]],
          ylab="population",
          xlab="population",
          main="Pairwise G'st",
          col = cols)
    gprimest<-as.vector(x[[3]][[5]])
    gprimest<-as.vector(na.omit(gprimest))
    collab111<-c(round(min(gprimest),3),
                 round(mean(gprimest),3),
                 round(max(gprimest),3))
    
    color.legend(nrow(x[[3]][[5]])/5,
                 nrow(x[[3]][[5]])/3,
                 nrow(x[[3]][[5]])/4,
                 nrow(x[[3]][[5]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    #Jost's D
    image(1:nrow(x[[3]][[6]]),
          1:nrow(x[[3]][[6]]),
          x[[3]][[6]],
          ylab="population",
          xlab="population",
          main="Pairwise Jost's D",
          col = cols,
          las=1)
    D<-as.vector(x[[3]][[6]])
    D<-as.vector(na.omit(D))
    collab111<-c(round(min(D),3),
                 round(mean(D),3),
                 round(max(D),3))
    
    color.legend(nrow(x[[3]][[6]])/5,
                 nrow(x[[3]][[6]])/3,
                 nrow(x[[3]][[6]])/4,
                 nrow(x[[3]][[6]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
  }
  if(exists("abx", where=".GlobalEnv")==TRUE){
    rm(abx, pos=".GlobalEnv")
  }
  if(exists("collab111", where=".GlobalEnv")==TRUE){
    rm(collab111, pos=".GlobalEnv")
  }
  if(exists("colleer", where=".GlobalEnv")==TRUE){
    rm(colleer, pos=".GlobalEnv")
  }
  if(exists("sp.header", where=".GlobalEnv")==TRUE){
    rm(sp.header, pos=".GlobalEnv")
  }
}
################################################################################
# end dif.Plot                                                                 #
################################################################################
#
#
#
################################################################################
###############################     END ALL       ##############################
################################################################################
