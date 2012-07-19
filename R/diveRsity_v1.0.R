#
# diversity.partition
# START
div.part<-function(infile,outfile=NULL, gp=3, 
                              bs_locus=FALSE, bs_pairwise=FALSE, 
                              bootstraps=0, Plot=FALSE){
  ############################ Argument definitions ############################
  D=infile
  on=outfile
  gp=gp
  bstrps=bootstraps
  bsls=bs_locus
  bspw=bs_pairwise
  plt=Plot
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
  #source("pre.div.R")
  Data<-pre.div(D,gp,F)
  # create a directory for output
  suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                         "-[diveRsity]","/",sep="")))
  of=paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
  wd<-getwd()
  write_res<-is.element("xlsx",installed.packages()[,1])
  plot_res<-is.element("sendplot",installed.packages()[,1])
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
  
  ##############################################################################
  # output file
  # multilocus stats vector
  # pre output table for global locus stats
  #standard
  pre_ot1<-cbind(Data$locus_names,round(as.numeric(Data$hst),4),
                 round(as.numeric(Data$dst),4),
                 round(as.numeric(Data$gst),4),
                 round(as.numeric(Data$gst_hedrick),4),
                 round(as.numeric(Data$djost),4))
  # Add global multi locus stats to output table
  ot1<-rbind(pre_ot1,c("Global","","",Data[[18]],Data[[20]],Data[[21]]))
  colnames(ot1)<-c("loci","H_st","D_st","G_st","G_hed_st","D_jost")
  #Estimated
  pre_ot2<-cbind(Data$locus_names,round(as.numeric(Data$locus_harmonic_N),4),
                 round(as.numeric(Data$ht_est),4),
                 round(as.numeric(Data$dst_est),4),
                 round(as.numeric(Data$gst_est),4),
                 round(as.numeric(Data$gst_est_hedrick),4),
                 round(as.numeric(Data$djost_est),4))
  ot2<-rbind(pre_ot2,c("Global","","","",Data[[24]],Data[[27]],Data[[28]]))
  colnames(ot2)<-c("loci","Harmonic_N","H_st_est","D_st_est","G_st_est",
                   "G_hed_st_est","D_Jost_est")
  
  plot_data<-c("Overall","","","",Data[[24]],Data[[27]],Data[[28]])
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
  ##############################################################################
  ############################## Bootstrapper ##################################
  ##############################################################################
  
  # Used only if bootstraps is greater than zero
  if(bsls==T){
    # Create a bootstrap output object for locus stats
    # bs_std = standard locus statistics
    bs_std<-list()
    for (i in 1:Data$nloci){
      bs_std[[i]]<-matrix(rep(0,(3*bstrps)),ncol=3,nrow=bstrps)
    }
    # bs_est = estimated locus statistics
    bs_est<-list()
    for (i in 1:Data$nloci){
      bs_est[[i]]<-matrix(rep(0,(3*bstrps)),ncol=3,nrow=bstrps)
    }
    # bs_glb = global stats 
    bs_glb<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
    ############################################################################
    # run bootstrap code
    # bstrps = 100
    # D= "Test_data.txt"
    # gp = 3
    for(w in 1:bstrps){
      data2<-pre.div(D,gp,BS=T,locus_stats=T)
      bs_glb[w,]<-c(round(data2$gst_all,4),
                    round(data2$gst_all_hedrick,4),
                    round(data2$djost_all,4),
                    round(data2$gst_est_all,4),
                    round(data2$gst_est_all_hedrick,4),
                    round(data2$djost_est_all,4))
      for(i in 1:data2$nloci){
        bs_std[[i]][w,]<-c(round(data2$gst[i],4),
                           round(data2$gst_hedrick[i],4),
                           round(data2$djost[i],4))
        bs_est[[i]][w,]<-c(round(data2$gst_est[i],4),
                           round(data2$gst_est_hedrick[i],4),
                           round(data2$djost_est[i],4))
      }
    }
    # Bootstrap output
    # bs_res = bootstrap results
    bs_res<-list()
    for(i in 1:6){
      bs_res[[i]]<-matrix(rep(0,((1+Data$nloci)*3)),ncol=3,
                          nrow=(Data$nloci+1))
    }
    sp=rep(1:3,2)
    ############################################################################
    #                 Try Manly(1997) method for CI calculation                #
    ############################################################################
    # Standard stats
    for(i in 1:3){    
      for(j in 1:Data$nloci){
        bs_res[[i]][j,1]<- as.numeric(ot1[j,(i+3)])
        bs_res[[i]][j,2]<- round(as.numeric(ot1[j,(i+3)])-
                                 (1.96*(sd(bs_std[[j]][,sp[i]]))),4)
        bs_res[[i]][j,3]<- round(as.numeric(ot1[j,(i+3)])+
                                 (1.96*(sd(bs_std[[j]][,sp[i]]))),4)
      }
      bs_res[[i]][(Data$nloci+1),1]<-as.numeric(ot1[(Data$nloci+1),(i+3)])
      bs_res[[i]][(Data$nloci+1),2]<-round(as.numeric(ot1[(Data$nloci+1),(i+3)])-
                                           (1.96*(sd(bs_glb[,i]))),4)
      bs_res[[i]][(Data$nloci+1),3]<-round(as.numeric(ot1[(Data$nloci+1),(i+3)])+
                                           (1.96*(sd(bs_glb[,i]))),4)
    }
    # estimated stats
    for(i in 4:6){
      for(j in 1:Data$nloci){
        bs_res[[i]][j,1]<- as.numeric(ot2[j,(i+1)])
        bs_res[[i]][j,2]<- round(as.numeric(ot2[j,(i+1)])-
                                 (1.96*(sd(bs_est[[j]][,sp[i]]))),4)
        bs_res[[i]][j,3]<- round(as.numeric(ot2[j,(i+1)])+
                                 (1.96*(sd(bs_est[[j]][,sp[i]]))),4)
      }
      bs_res[[i]][(Data$nloci+1),1]<-as.numeric(ot2[(Data$nloci+1),(i+1)])
      bs_res[[i]][(Data$nloci+1),2]<-round(as.numeric(ot2[(Data$nloci+1),(i+1)])-
                                           (1.96*(sd(bs_glb[,i]))),4)
      bs_res[[i]][(Data$nloci+1),3]<-round(as.numeric(ot2[(Data$nloci+1),(i+1)])+
                                           (1.96*(sd(bs_glb[,i]))),4)
    }
 ############################################################################
 #                      Alternative 95% ci calculation                      #
 ############################################################################
 #for(i in 1:3){    
 #  for(j in 1:Data$nloci){
 #    bs_res[[i]][j,1]<- round(mean(as.numeric(bs_std[[j]][,sp[i]])),4)
 #    bs_res[[i]][j,2]<- round(quantile(bs_std[[j]][,sp[i]],
 #                                      probs=c(0.025,0.975))[[1]],4)
 #    bs_res[[i]][j,3]<- round(quantile(bs_std[[j]][,sp[i]],
 #                                      probs=c(0.025,0.975))[[2]],4)
 #  }
 #  bs_res[[i]][(Data$nloci+1),1]<- round(mean(as.numeric(bs_glb[,i])),4)
 #  bs_res[[i]][(Data$nloci+1),2]<- round(quantile(bs_glb[,i],
 #                                                 probs=c(0.025,0.975))[[1]],4)
 #  bs_res[[i]][(Data$nloci+1),3]<- round(quantile(bs_glb[,i],
 #                                                 probs=c(0.025,0.975))[[2]],4)
 #}
 # estimated stats
 #for(i in 4:6){
 #  for(j in 1:Data$nloci){
 #    bs_res[[i]][j,1]<- round(mean(as.numeric(bs_est[[j]][,sp[i]])),4)
 #    bs_res[[i]][j,2]<- round(quantile(bs_est[[j]][,sp[i]],
 #                                      probs=c(0.025,0.975))[[1]],4)
 #    bs_res[[i]][j,3]<- round(quantile(bs_est[[j]][,sp[i]],
 #                                      probs=c(0.025,0.975))[[2]],4)
 #  }
 #  bs_res[[i]][(Data$nloci+1),1]<- round(mean(as.numeric(bs_glb[,i])),4)
 #  bs_res[[i]][(Data$nloci+1),2]<- round(quantile(bs_glb[,i],
 #                                                 probs=c(0.025,0.975))[[1]],4)
 #  bs_res[[i]][(Data$nloci+1),3]<- round(quantile(bs_glb[,i],
 #                                                 probs=c(0.025,0.975))[[2]],4)
 #}
   
    # bs_res list element names
    names(bs_res)<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
                     "D_Jost_est")
    bs_res1<-bs_res
    for(i in 1:6){
      dimnames(bs_res1[[i]])<-list(c(Data$locus_names,"global"),
                                   c("Actual","Lower_CI","Upper_CI"))
    }
    # bs results output object header
    hdr<-matrix(c("locus","Actual","Lower_95%CI","Upper_95%CI"),ncol=4)
    bs_out<-matrix(rbind(hdr,c(names(bs_res)[1],"","",""),
                         cbind(c(Data$locus_names,"Overall"),
                               bs_res[[1]])),ncol=4)
    for(i in 2:6){
      bs_out<-matrix(rbind(bs_out,c(names(bs_res)[i],"","",""),
                           cbind(c(Data$locus_names,"Global"),
                                 bs_res[[i]])),ncol=4)
    }
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
  # bs_res<-print(bs_res)
  ##########################################################################
  if(plot_res==TRUE){
  if(bsls==T && plt==T){    
    lso<-list()
    for(i in 1:6){
      lso[[i]]<-order(bs_res[[i]][1:Data$nloci,1],decreasing=F)
    }
    names(lso)<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
    		  "D_Jost_est")
    plot.call_loci<-list()
    plot.extras_loci<-list()
    xy.labels_loci<-list()
    y.pos_loci<-list()
    x.pos_loci=1:Data$nloci
    direct=of
    fn_pre_loci<-list()
    #Plot Gst_Nei
    plot.call_loci[[1]]=c("plot(bs_res[[4]][1:Data$nloci,1],
                          ylim=c(0,(max(bs_res[[4]][,3])+
                          min(bs_res[[4]][,3]))),xaxt='n',
                          ylab=names(bs_res)[4],type='n',
                          xlab='Loci \n (Hover over a point to see locus data)',
                          cex.lab=1.5,cex.axis=1.3,las=1)")

    plot.extras_loci[[1]]=c("points(bs_res[[4]][lso[[4]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:Data$nloci,bs_res[[4]][lso[[4]],2],
                            1:Data$nloci,bs_res[[4]][lso[[4]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=c(0,bs_res[[4]][(Data$nloci+1),2]),
                            lwd=1,lty=c(1,2),col=c('black','red'))")
  
    xy.labels_loci[[1]]=data.frame(Locus_name=Data$locus_names[lso[[4]]],
                                   Gst_Nei=round(bs_res[[4]][lso[[4]],1],4),
                                   Gst_Hedrick=round(bs_res[[5]][lso[[4]],1],4),
                                   D_jost=round(bs_res[[6]][lso[[4]],1],4))
    
    y.pos_loci[[1]]=bs_res[[4]][lso[[4]],1]
    fn_pre_loci[[1]]<-names(bs_res)[4]
    
    
    
    # Plot Gst_Hedrick
    plot.call_loci[[2]]=c("plot(bs_res[[5]][1:Data$nloci,1],
                          ylim=c(0,1),xaxt='n',ylab=names(bs_res)[5],type='n',
                          xlab='Loci \n (Hover over a point to see locus data)',
                          cex.lab=1.5,cex.axis=1.3,las=1)")

    plot.extras_loci[[2]]=c("points(bs_res[[5]][lso[[5]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:Data$nloci,bs_res[[5]][lso[[5]],2],
                            1:Data$nloci,bs_res[[5]][lso[[5]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=c(0,bs_res[[5]][(Data$nloci+1),2]),
                            lwd=1,lty=c(1,2),col=c('black','red'))")
  
    xy.labels_loci[[2]]=data.frame(Locus_name=Data$locus_names[lso[[5]]],
                                   Gst_Nei=round(bs_res[[4]][lso[[4]],1],4),
                                   Gst_Hedrick=round(bs_res[[5]][lso[[5]],1],4),
                                   D_jost=round(bs_res[[6]][lso[[5]],1],4))
    
    y.pos_loci[[2]]=bs_res[[5]][lso[[5]],1]
    fn_pre_loci[[2]]<-names(bs_res)[5]
    
    
    # Plot D_jost
    plot.call_loci[[3]]=c("plot(bs_res[[6]][1:Data$nloci,1],
                          ylim=c(0,1),xaxt='n',ylab=names(bs_res)[6],type='n',
                          xlab='Loci \n (Hover over a point to see locus data)',
                          cex.lab=1.5,cex.axis=1.3,las=1)")

    plot.extras_loci[[3]]=c("points(bs_res[[6]][lso[[6]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:Data$nloci,bs_res[[6]][lso[[6]],2],
                            1:Data$nloci,bs_res[[6]][lso[[6]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=c(0,bs_res[[6]][(Data$nloci+1),2]),
                            lwd=1,lty=c(1,2),col=c('black','red'))")
  
    xy.labels_loci[[3]]=data.frame(Locus_name=Data$locus_names[lso[[6]]],
                                   Gst_Nei=round(bs_res[[4]][lso[[4]],1],4),
                                   Gst_Hedrick=round(bs_res[[5]][lso[[5]],1],4),
                                   D_jost=round(bs_res[[6]][lso[[6]],1],4))
    
    y.pos_loci[[3]]=bs_res[[6]][lso[[6]],1]
    fn_pre_loci[[3]]<-names(bs_res)[6]
    
    }
  }
  ############################################################################
  ################################## Pairwise ################################
  ############################################################################
  if(bspw==T){
    # population pair combinations
    pw<-combn(Data$npops,2)
    # Bootstrap results data object 
    # bs_pw_glb = bootstrap pairwise global stats
    bs_pw_glb<-matrix(rep(0,(6*bstrps)),ncol=6,nrow=bstrps)
    # output results data object
    # pw_res = pairwise results
    pw_res<-list()
    for(i in 1:6){
      pw_res[[i]]<-matrix(rep(0,(ncol(pw)*3)),ncol=3)
    }
    data1<-readGenepop(D,gp,F)
    
    pairwise_combo<-combn(data1$npops,2)
    
    
    ind_vectors<-list()
    for(i in 1:data1$npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,data1$pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("",((data1$nloci+1)*(data1$nloci+1))),
                     ncol=(data1$nloci+1))
    pre_data[1,]<-rep("",(data1$nloci+1))
    for(i in 2:(data1$nloci+1)){
      pre_data[i,1]<-data1$loci_names[(i-1)]
    }
    
    pw_data<-list()
    for (i in 1:ncol(pairwise_combo)){
      pw_data[[i]]<-data.frame(rbind(pre_data,
                               c("POP",as.vector(rep("",data1$nloci))),
                               cbind(ind_vectors[[pairwise_combo[1,i]]],
                               matrix(noquote(data1$pop_list
                                              [[pairwise_combo[1,i]]]),
                                               ncol=data1$nloci)),
                               c("POP",as.vector(rep("",data1$nloci))),
                               cbind(ind_vectors[[pairwise_combo[2,i]]],
                               matrix(noquote(data1$pop_list
                                              [[pairwise_combo[2,i]]]),
                                                ncol=data1$nloci))))
    }
    pw_glb<-matrix(rep(0,(6*(ncol(pairwise_combo)))),ncol=6)
    for(i in 1:ncol(pairwise_combo)){
      dat<-pre.div(pw_data[[i]],3,F,F)
      pw_glb[i,]<-c(dat$gst_all,dat$gst_all_hedrick,
                    dat$djost_all,dat$gst_est_all,
                    dat$gst_est_all_hedrick,dat$djost_est_all)
    }
    for(i in 1:ncol(pairwise_combo)){
      for(w in 1:bstrps){
        pw_stat<-pre.div(pw_data[[i]],gp,T,T)
        bs_pw_glb[w,]<-c(pw_stat$gst_all,pw_stat$gst_all_hedrick,
                         pw_stat$djost_all,pw_stat$gst_est_all,
                         pw_stat$gst_est_all_hedrick,pw_stat$djost_est_all)
      }
      for(j in 1:6){
        pw_res[[j]][i,1]<-pw_glb[i,j]
        pw_res[[j]][i,2]<-round(pw_glb[i,j]-(1.96*(sd(bs_pw_glb[,j]))),4)
        pw_res[[j]][i,3]<-round(pw_glb[i,j]+(1.96*(sd(bs_pw_glb[,j]))),4)
      }
    }
    # pairwise comparisons
    # pw_names = pairwise population names
    pw_nms<-paste(Data$pop_names[pairwise_combo[1,]],
                  Data$pop_names[pairwise_combo[2,]],sep=" vs. ")
    pw_nms1<-paste(pairwise_combo[1,],pairwise_combo[2,],sep=" vs. ")
    names(pw_res)<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
    		     "D_Jost_est")
    pw_res1<-pw_res
    for(i in 1:6){
      dimnames(pw_res1[[i]])<-list(pw_nms,c("Actual","Lower_CI","Upper_CI"))
    }
    # bs results output object header
    hdr<-matrix(c("Pairwise","Actual","Lower_95%CI","Upper_95%CI"),ncol=4)
    pw_bs_out<-matrix(rbind(hdr,c(names(pw_res)[1],"","",""),
                            cbind(pw_nms,pw_res[[1]])),ncol=4)
    for(i in 2:6){
      pw_bs_out<-matrix(rbind(pw_bs_out,c(names(pw_res)[i],"","",""),
                              cbind(pw_nms,pw_res[[i]])),ncol=4)
    }
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
  ############################################################################
  #pw plotter
  if(plot_res==TRUE){
  if(bspw==T && plt==T){      
      pwso<-list()
      for(i in 1:6){
        pwso[[i]]<-order(pw_res[[i]][,1],decreasing=F)
      }
      names(pwso)<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
                     "D_Jost_est")
      # define plot parameters 
      plot.call_pw<-list()
      plot.extras_pw<-list()
      xy.labels_pw<-list()
      y.pos_pw<-list()
      x.pos_pw=1:ncol(pw)
      fn_pre_pw<-list()
      direct=of
      #Plot Gst_Nei
      plot.call_pw[[1]]=c("plot(pw_res[[4]][1:ncol(pw),1],
                            ylim=c(0,(max(pw_res[[4]][,3])+
                          min(pw_res[[4]][,3]))),xaxt='n',
                          ylab=names(pw_res)[4],type='n',
                          xlab='Pairwise comparisons 
                              \n (Hover over a point to see pairwise data)',
                        cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_pw[[1]]=c("points(pw_res[[4]][pwso[[4]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:ncol(pw),pw_res[[4]][pwso[[4]],2],
                            1:ncol(pw),pw_res[[4]][pwso[[4]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data[5]),
                            lwd=1,lty=2,col='red')")
  
      xy.labels_pw[[1]]=data.frame(pairwise_name=pw_nms[pwso[[4]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[4]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[4]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[4]],1],4))
      
      y.pos_pw[[1]]=pw_res[[4]][pwso[[4]],1]
      fn_pre_pw[[1]]<-names(pw_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_pw[[2]]=c("plot(pw_res[[5]][1:ncol(pw),1],
                            ylim=c(0,1),xaxt='n',ylab=names(pw_res)[5],type='n',
                          xlab='Pairwise comparisons
                              \n (Hover over a point to see locus data)',
                        cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_pw[[2]]=c("points(pw_res[[5]][pwso[[5]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:ncol(pw),pw_res[[5]][pwso[[5]],2],
                            1:ncol(pw),pw_res[[5]][pwso[[5]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data[6]),
                            lwd=1,lty=2,col='red')")
  
      xy.labels_pw[[2]]=data.frame(pairwise_name=pw_nms[pwso[[5]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[4]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[5]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[5]],1],4))
      
      y.pos_pw[[2]]=pw_res[[5]][pwso[[5]],1]
      fn_pre_pw[[2]]<-names(pw_res)[5]
      
      
      # Plot D_jost
      plot.call_pw[[3]]=c("plot(pw_res[[6]][1:ncol(pw),1],
                            ylim=c(0,1),xaxt='n',ylab=names(pw_res)[6],type='n',
                          xlab='Pairwise comparisons 
                             \n (Hover over a point to see locus data)',
                        cex.lab=1.5,cex.axis=1.3,las=1)")

      plot.extras_pw[[3]]=c("points(pw_res[[6]][pwso[[6]],1],
                              pch=15,col='black',cex=1);
                            arrows(1:ncol(pw),pw_res[[6]][pwso[[6]],2],
                            1:ncol(pw),pw_res[[6]][pwso[[6]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data[7]),
                            lwd=1,lty=2,col='red')")
    
      xy.labels_pw[[3]]=data.frame(pairwise_name=pw_nms[pwso[[6]]],
                                  Gst_Nei=round(pw_res[[4]][pwso[[4]],1],4),
                                  Gst_Hedrick=round(pw_res[[5]][pwso[[5]],1],4),
                                  D_jost=round(pw_res[[6]][pwso[[6]],1],4))
      
      y.pos_pw[[3]]=pw_res[[6]][pwso[[6]],1]
      fn_pre_pw[[3]]<-names(pw_res)[6]
    }
  }
  ############################### Bootstrap end ################################
  
  
  ################################# Plot resuts ################################
  #make necessary data available
  if(plt==T &&  plot_res==T){
  if(bsls==T && bspw==T){
    pl<-list(bs_res=bs_res,
             pw_res=pw_res,
             Data=Data,
             lso=lso,
             pwso=pwso,
             plot.call_loci=plot.call_loci,
             plot.extras_loci=plot.extras_loci,
             xy.labels_loci=xy.labels_loci,
             x.pos_loci=x.pos_loci,
             y.pos_loci=y.pos_loci,
             fn_pre_loci=fn_pre_loci,
             direct=direct,
             plot_loci=bsls,
             plot_pw=bspw,
             plot.call_pw=plot.call_pw,
             plot.extras_pw=plot.extras_pw,
             xy.labels_pw=xy.labels_pw,
             y.pos_pw=y.pos_pw,
             fn_pre_pw=fn_pre_pw,
             x.pos_pw=x.pos_pw,
             pw=pw,plot_data=plot_data)
  } else if (bsls==T & bspw==F){
    pl<-list(bs_res=bs_res,
             Data=Data,
             lso=lso,
             plot.call_loci=plot.call_loci,
             plot.extras_loci=plot.extras_loci,
             xy.labels_loci=xy.labels_loci,
             x.pos_loci=x.pos_loci,
             y.pos_loci=y.pos_loci,
             fn_pre_loci=fn_pre_loci,
             direct=direct,
             plot_loci=bsls,
             plot_pw=bspw,plot_data=plot_data)
  } else if (bsls==F && bspw==T){
    pl<-list(pw_res=pw_res,
             Data=Data,
             pwso=pwso,
             plot.call_pw=plot.call_pw,
             plot.extras_pw=plot.extras_pw,
             xy.labels_pw=xy.labels_pw,
             x.pos_pw=x.pos_pw,
             y.pos_pw=y.pos_pw,
             fn_pre_pw=fn_pre_pw,
             direct=direct,
             plot_loci=bsls,
             plot_pw=bspw,
             pw=pw,plot_data=plot_data)
  }
    #source("plotter.R")
    suppressWarnings(plotter(pl,"1000x600"))
  }
  
  
  #############################################################################
  #Data for output
  if(bspw==T && bsls==T){
    list(standard=noquote(ot1),
         estimate=noquote(ot2),
         bs_locus=bs_res1,
         bs_pairwise=pw_res1)
  } else if(bspw==T && bsls==F){
    list(standard=noquote(ot1),
         estimate=noquote(ot2),
         bs_pairwise=pw_res1)
  } else if(bspw==F && bsls==T){
    list(standard=noquote(ot1),
         estimate=noquote(ot2),
         bs_locus=bs_res1)
  }
 }
}
################################################################################
# diversity.partition
# END
#
#
#
# in.calc
# START
################################################################################
in.calc<-function(infile, outfile=NULL, gp=3, bs_locus=FALSE, bs_pairwise=FALSE, 
                  bootstraps=0, Plot=FALSE){
  D=infile
  gp=gp
  pw=bs_pairwise
  BS=bs_locus
  NBS=bootstraps
  on=outfile
  plt=Plot
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
  suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                         "-[diveRsity]","/",sep="")))
  
  of<-paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
  
  #source("in.bootstrap.R")
  res_out<-in.bs(D,gp,F,0)[[1]]
  if(write_res==TRUE){
  write.xlsx2(res_out,file=paste(of,"In_out.xlsx",sep=""),
             sheetName="In_allele_stats",col.names=T,row.names=T,append=F)
  } else {
    all_out<-file(paste(of,"Allele-In[in.calc].txt",sep=""),"w")
    cat(paste(colnames(res_out),sep=""),"\n",sep="\t",file=all_out)
    for(i in 1:nrow(res_out)){
      cat(res_out[i,],"\n",sep="\t",file=all_out)
    }
    close(all_out)
  }
  ######################################################################
  # overall In
  if(BS==T){
    bs_sum1<-in.bs(D,gp,BS,NBS)[[2]]
    if(write_res==T){
    write.xlsx(bs_sum1,file=paste(of, "In_out.xlsx",sep=""),
               sheetName="Overall_Bootstrap",col.names=T,row.names=F,append=T)
    } else {
      all_bs<-file(paste(of,"Overall-bootstrap[in.calc].txt",sep=""),"w")
      cat(paste(colnames(bs_sum1),sep=""),"\n",sep="\t",file=all_bs)
      for(i in 1:nrow(bs_sum1)){
        cat(bs_sum1[i,],"\n",sep="\t",file=all_bs)
      }
      close(all_bs)
    }
    loc_nms<-bs_sum1[,1]
    bs_sum<-matrix(as.numeric(bs_sum1[,2:4]),ncol=3)
    if(plt==T){
    lso<-order(bs_sum[,1],decreasing=F)
    png(filename=paste(of, on,"_In_plot.png",sep=""),width=800,height=600)
    par(mar=c(6,5,1,1))
    plot(bs_sum[lso,1],ylim=c(0,(max(bs_sum[,3])+0.1)),xaxt='n',
                              ylab=expression('Locus '*I[n]),
                              xlab="",cex.lab=1.5,cex.axis=1.3,las=1,type='n')
         points(bs_sum[lso,1],pch=15,col='black',cex=1)
         arrows(1:nrow(bs_sum),bs_sum[lso,2],1:nrow(bs_sum),bs_sum[lso,3],
                code=3,angle=90,length=0.05,lwd=0.1)
         axis(1,at=1:nrow(bs_sum),labels=loc_nms[lso],las=3)
         dev.off()
    }
  }
  # pairwise locus In bootstrap
  if(pw==T){
    #source("readGenepop.R")
    data<-readGenepop(D,gp,F)
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
    for(i in 1:ncol(pwc)){
      pw_bs[[i]]<-matrix(in.bs(pw_data[[i]],gp,pw,NBS)[[2]],ncol=4)
    }
    pw_nms<-paste(pn[pwc[1,]],pn[pwc[2,]],sep=" vs. ")
    names(pw_bs)<-pw_nms
    pw_bs1<-pw_bs
    hdr<-c("Loci","Actual_In","Lower_95CI","Upper_95CI")
    pw_in_bs<-matrix(rbind(hdr,c(names(pw_bs)[1],"","",""),pw_bs[[1]]),ncol=4)
    for(j in 2:ncol(pwc)){
      pw_in_bs<-matrix(rbind(pw_in_bs,c(names(pw_bs)[j],"","",""),pw_bs[[j]]),
                       ncol=4)
    }
    if(write_res==TRUE){
    write.xlsx(pw_in_bs,file=paste(of, on,"_In_out.xlsx",sep=""),
               sheetName="Pairwise_bootstraps",col.names=F,row.names=F,append=T)
    } else {
      pw_bs<-file(paste(of,"Pairwise-bootstrap[in.calc].txt",sep=""),"w")
      cat(paste(colnames(pw_in_bs),sep=""),"\n",sep="\t",file=pw_bs)
      for(i in 1:nrow(pw_in_bs)){
        cat(pw_in_bs[i,],"\n",sep="\t",file=pw_bs)
      }
      close(pw_bs)
    } 
  }
  
  if(BS==F && pw==F){
    list(Allele_In=res_out)
  } else if (BS==T && pw==F){
    list(Allele_In=res_out,
         l_bootstrap=bs_sum1)
  } else if (BS==F && pw==T){
    list(Allele_In=res_out,
         PW_bootstrap=pw_in_bs)
  } else if (BS==T && pw==T){
    list(Allele_In=res_out,
         l_bootstrap=bs_sum1,
         PW_bootstrap=pw_in_bs)
  }      
 }
}
################################################################################
# in.calc
# END
#
#
#
# in.bootsrap
# START
################################################################################
in.bs<-function(D,gp,BS,NBS){
  D=D
  gp=gp
  BS=BS
  NBS=NBS
  #source("readGenepop.R")
  data<-readGenepop(D,gp,F)
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
  for(w in 1:NBS){
    pre_data<-readGenepop(D,gp,T)
    data_bs<-readGenepop(pre_data$bs_file,gp,F)
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
  in_bs_out<-matrix(rep(0,(nl_bs*4)),ncol=4)
  colnames(in_bs_out)<-c("Locus","In","Lower_95CI","Upper_95CI")
  in_bs_out[,1]<-ln_bs
  in_bs_out[,2]<-round(results_out[,(max(nal)+1)],4)
  for(i in 1:nl){
  in_bs_out[i,3]<-round(results_out[i,(max(nal)+1)]-(1.96*(sd(bs_sum[,i]))),
                        4)
  in_bs_out[i,4]<-round(results_out[i,(max(nal)+1)]+(1.96*(sd(bs_sum[,i]))),
                        4)
    }
  }
  colnames(results_out)<-c(paste("Allele.",1:max(nal),sep=""),"Sum")
  rownames(results_out)<-ln
  results_out[is.na(results_out)]<-""
  if(BS==T){
    list(In_alleles=results_out,
         in_bs_out=in_bs_out)
  }else if(BS==F){
    list(In_alleles=results_out)
  }
}
################################################################################
# in.bootsrap
# END
#
#
#
# readGenepop
# START
################################################################################
readGenepop <- function (infile = NULL, gp = 3, bootstrap = FALSE) {
  Genepop_format=gp
  infile=infile
  bootstrap=bootstrap
  if(typeof(infile)=="list"){
    data1=infile 
  } else if (typeof(infile)=="character"){
    no_col <- max(count.fields(infile))
    data1 <- read.delim(infile,fill=T,col.names=1:no_col,header=F)
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
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,4)
  pop_weights<-vector()
  for (i in 1: npops){
    pop_weights[i]<- 1/pop_sizes[i]
  }
  n_harmonic<-npops/sum(pop_weights)
  #if (typeof(N)=="NULL"){
  N<-pop_sizes
  #}else{                      # used in sampling simulations
  #  N<-rep(N,npops)
  #}
  #pop_file_names<-paste("pop-",1:npops,sep="")
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
    2:(nloci+1)])
  }
  if (Genepop_format==3){
    for (i in 1:npops){
      pop_list[[i]]<-matrix(sprintf("%06g",as.numeric(pop_list[[i]])),
                                 nrow=pop_sizes[i],ncol=nloci)
    }
  } else if (Genepop_format==2){
    for (i in 1:npops){
      pop_list[[i]]<-matrix(sprintf("%04g",as.numeric(pop_list[[i]])),
                                           nrow=pop_sizes[i],ncol=nloci)
    }
  }
  for(i in 1:npops){
    pop_list[[i]][pop_list[[i]]=="    NA"]<-NA
  }
  if (bootstrap == T){
    for(i in 1:npops){
      pop_list[[i]]<-matrix(pop_list[[i]][sample(N[i],replace=T),],ncol=nloci)
    }
  }
  loci_pop_sizes<-list()
  for(i in 1:nloci){
    loci_pop_sizes[[i]]<-vector()
  }
  for (i in 1:nloci){
    for (j in 1:npops){
      loci_pop_sizes[[i]][j]<- length(na.omit(pop_list[[j]][,i]))
    }
  }
  loci_harm_N<-vector()
  loci_pop_weights<-list()
  for (i in 1:nloci){
    loci_pop_weights[[i]]<-vector()
  }    
  for (i in 1:npops){
    for (j in 1:nloci){
      loci_pop_weights[[j]][i]<- 1/loci_pop_sizes[[j]][i]
    }
  }
  loci_harm_N<-vector()
  for (i in 1:nloci){
    loci_harm_N[i]<- npops/sum(loci_pop_weights[[i]])
  }
  pop_alleles<-list()
  x<-seq(1,(2*npops),2)
  y<-seq(2,(2*npops),2)
  if (Genepop_format==3){
    for (i in 1:npops){
      for (j in 1:nloci){
        pop_alleles[[x[i]]]<-substr(pop_list[[i]],1,3)
        pop_alleles[[y[i]]]<-substr(pop_list[[i]],4,6)
      }
    }
  } else if (Genepop_format==2) {
    for (i in 1:npops){
      for (j in 1:nloci){
        pop_alleles[[x[i]]]<-substr(pop_list[[i]],1,2)
        pop_alleles[[y[i]]]<-substr(pop_list[[i]],3,4)
      }
    }
  }
  allele_names<-list()
  for (i in 1:npops){
    allele_names[[i]]<-list()
  }
  for (i in 1:npops){
    for (j in 1:nloci){
        allele_names[[i]][j]<-list(sort(unique(c(pop_alleles[[x[i]]][,j],
        pop_alleles[[y[i]]][,j])),decreasing=F))
    }
  }
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  all_alleles<-list()
  for (j in 1:nloci){
    all_alleles[[j]]<- sort(unique(loci_combi[[j]],decreasing=F))
  }
  allele_freq<-list()
  for (i in 1:nloci){
    allele_freq[[i]]<-matrix(rep(0,(npops*length(all_alleles[[i]]))),
                             ncol=npops)
    rownames(allele_freq[[i]])<-all_alleles[[i]]
  }
  for(i in 1:nloci){
    for(j in 1:length(all_alleles[[i]])){
      for(z in 1:npops){
        allele_freq[[i]][j,z]<-(length(which(pop_alleles[[x[z]]][,i]==
          all_alleles[[i]][j]))+length(which(pop_alleles[[y[z]]][,i]
                             ==all_alleles[[i]][j])))/
                            (2*length(na.omit(pop_alleles[[x[z]]][,i])))                                       
      }
    }
  }
  # indtyp = individuals typed (per locus)
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:nloci){
    for(j in 1:npops){
      indtyp[[i]][j]<-length(na.omit(pop_alleles[[x[j]]][,i]))
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
  bs_data_file<-data.frame(bs_data_file,ncol=(nloci+1))}
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
         nalleles=nalleles)
  }
}
################################################################################
# readGenepop
# END
#
#
#
# pre.div
# START
################################################################################
pre.div<-function(Data_file, gp, BS=F,locus_stats=T){
  D=Data_file
  gp=gp
  BS=BS
  ls=locus_stats
  # load readGenepop file reader function
  #source("readGenepop.R")
  # read genepop data file
  data1<-readGenepop(D,gp,BS)
  #data1<-readGenepop(D,3,F)
  
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
  
  # creat index vectors
  x<-seq(1,(2*np),2)
  y<-seq(2,(2*np),2)
  
  # count the number of observed hets
  # ohc = observed heterozygote count
  ohc<-pl
  for(i in 1:np){
    for(j in 1:nl){
      for(z in 1:ps[i]){
        if (is.na(pa[[x[i]]][z,j]) == TRUE) {
          ohc[[i]][z,j]<-NA
        } else if (pa[[x[i]]][z,j] != pa[[y[i]]][z,j]) {
          ohc[[i]][z,j]<-1
        } else if (pa[[x[i]]][z,j] == pa[[y[i]]][z,j]) {
          ohc[[i]][z,j]<-0
        }
      }
    }
  }
  
  # Calculate expected homozygosity
  # function for calculating expected homozygote frequency (i.e. P^2)
  square<-function(x){x^2}
  # exhf = expected homo frequency
  # create empty data matrix
  exhmf<-list()
  for (i in 1:nl){
    exhmf[[i]]<-matrix(rep(0,np*length(ant[[i]])),ncol=np,
                       nrow=length(ant[[i]]))
  }
  # calculate expected homo freq using square function
  for(i in 1:nl){
    for(z in 1:length(ant[[i]])){
      for(j in 1:np){
        exhmf[[i]][z,j]<-square(as.numeric(af[[i]][z,j]))
      }
    }
  }
  
  # use exhf to calculate expected heterozygoyte frequency
  # exhtf = expected heterozygote frequency
  exhtf<-matrix(rep(0,(np*nl)),ncol=np,nrow=nl)
  for(i in 1:nl){
    for(j in 1:np){
      exhtf[i,j]<- 1- (sum(as.numeric(na.omit(exhmf[[i]][,j]))))
    }
  }
  
  # calculate mean frequencies
  # mf = mean frequency
  mf<-list()
  for (i in 1:nl){
    mf[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      mf[[i]][j]<-sum(af[[i]][j,])/np
    }
  }
  
  # calculate mean expected homozygote frequency using "squred" function
  # mexhmf = mean expected homozygote frequency
  mexhmf<-mf
  for(i in 1:nl){
    for(j in 1:length(mf[[i]])){
      mexhmf[[i]][j]<-square(as.numeric(mf[[i]][j]))
    }
  }
  
  # claculate total expected heterozygosity
  # ht = total exp heterozygosity 
  ht<-vector()
  for(i in 1:nl){
    ht[i]<- 1-(sum(as.numeric(mexhmf[[i]])))
  }
  
  # Locus stats #
  ##############################################################################
  # Define component stats
  hs<-vector();hst<-vector();dst<-vector();gst<-vector();djost<-vector()
  hs_est<-vector();ht_est<-vector();hst_est<-vector();dst_est<-vector()
  gst_est<-vector();djost_est<-vector();gst_max<-vector() 
  gst_est_max<-vector();gst_hedrick<-vector();gst_est_hedrick<-vector()
  ##############################################################################
  for(i in 1:nl){
    hs[i]<- round((sum(exhtf[i,])/np),4)
    hs_est[i]<-round(hs[i]*((2*lnharm[i])/((2*lnharm[i])-1)),4)
    ht_est[i]<-round((ht[i]+ (hs_est[i]/(2*lnharm[i]*np))),4)
  }
  if(ls==T){
    for(i in 1:nl){
      hst[i]<-(ht[i]-hs[i])/(1-hs[i])
      dst[i]<-(ht[i]-hs[i])
      gst[i]<-(dst[i]/ht[i])
      djost[i]<-((ht[i]-hs[i])/(1-hs[i]))*(np/(np-1))
      hst_est[i]<-(ht_est[i]-hs_est[i])/(1-hs_est[i])
      dst_est[i]<-ht_est[i]- hs_est[i]
      gst_est[i]<-(ht_est[i]-hs_est[i])/ht_est[i]
      gst_max[[i]]<-((np-1)*(1-hs[i]))/(np-1+hs[i])
      gst_est_max[i]<-(((np-1)*(1-hs_est[i]))/(np-1+hs_est[i]))
      gst_hedrick[i]<-gst[i]/gst_max[i]
      gst_est_hedrick[i]<-gst_est[i]/gst_est_max[i]
    }
  }
  ##############################################################################
  for(i in 1:nl){
    djost_est[i]<-(np/(np-1))*((ht_est[i]-hs_est[i])/(1 - hs_est[i]))
  }
  
  # Across all loci stats #
  ht_mean<-round(mean(ht),4)
  hs_mean<-round(mean(hs),4)
  gst_all<-round((ht_mean-hs_mean)/ht_mean,4)
  gst_all_max<-round(((np-1)*(1-hs_mean))/(np-1+hs_mean),4)
  gst_all_hedrick<-round(gst_all/gst_all_max,4)
  djost_all<-round(((ht_mean-hs_mean)/(1-hs_mean))*(np/(np-1)),4)
  ##############################################################################
  
  
  # Across all loci estimated stats #
  hs_est_mean<-round(mean(hs_est),4)
  ht_est_mean<-round(mean(ht_est),4)
  gst_est_all<-round((ht_est_mean-hs_est_mean)/ht_est_mean,4)
  gst_est_all_max<-round((((np-1)*(1-hs_est_mean))/(np-1+hs_est_mean)),4)
  gst_est_all_hedrick<-round(gst_est_all/gst_est_all_max,4)
  #djost_est_all<-round((np/(np-1))*((ht_est_mean-hs_est_mean)/
  #(1 - hs_est_mean)),4)
  djost_est_all<-round(1/((1/mean(djost_est)+(var(djost_est)*
                           ((1/mean(djost_est))^3)))),4)
  ##############################################################################


  ##############################################################################
  if(ls==T){
    out_put<-list(hs=hs,
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
    out_put<-list(ht_mean=ht_mean,
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
    return(out_put)
}
################################################################################
# pre.div
# END
#
#
#
# plotter
# START
################################################################################
plotter<-function(x,img="1000x600"){
  image.size=img
  spot.radius=5
  jjj<-list()
  jjj<<-x
  if(!require("sendplot")){
    install.packages("sendplot",dependencies=TRUE)
  }
  require("sendplot")
  fl_ext<-c(".tif","Dot.png","Dot.tif")
  bs_res<-list()
  lso<-list()
  Data<-list()
  sp.header<-list()
  pw_res<-list()
  pwso<-list()
  pw<-list()
  plot_data<-list()
  if(jjj$plot_loci==T && jjj$plot_pw==F){
    bs_res<<-jjj$bs_res
    lso<<-jjj$lso
    Data<<-jjj$Data
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
    rm(jjj,Data,bs_res,lso,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==F && jjj$plot_pw==T){
    Data<<-jjj$Data
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data<<-jjj$plot_data
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
    
    rm(jjj,Data,plot_data,pw,pw_res,pwso,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==T && jjj$plot_pw==T){
    bs_res<<-jjj$bs_res
    lso<<-jjj$lso
    Data<<-jjj$Data
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data<<-jjj$plot_data
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
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",
    	         fl_ext,sep=""))
    rm(jjj,Data,plot_data,pw,pw_res,bs_res,lso,pwso,sp.header,pos=".GlobalEnv")
  }
}
################################################################################
# plotter
# END
################################################################################
# diveRsity
# END ALL
