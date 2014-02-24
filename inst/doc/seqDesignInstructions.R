### R code from vignette source 'seqDesignInstructions.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seqDesignInstructions.Rnw:31-34 (eval = FALSE)
###################################################
## N.pla <- 1900
## N.vax <- 1100
## N.vax.arms <- 4


###################################################
### code chunk number 2: seqDesignInstructions.Rnw:38-106 (eval = FALSE)
###################################################
## aveVElist <- list(-2, -1.5, -1, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
## aveVElist[12:19] <- lapply(aveVElist[-(1:3)], function(aveVE){ rep(aveVE, 4) })
## aveVElist[[20]] <- rep(0.5, 2)
## aveVElist[[21]] <- rep(0.5, 3)
## aveVElist[[22]] <- c(0,   0,    0,    0.4)
## aveVElist[[23]] <- c(0,   0,    0.3,  0.4)
## aveVElist[[23]] <- c(0.2, 0.2,  0.3,  0.4)
## aveVElist[[24]] <- c(0,   0,    0,    0.6)
## aveVElist[[25]] <- c(0,   0,    0.3,  0.6)
## aveVElist[[26]] <- c(0,   0,    0.45, 0.6)
## aveVElist[[27]] <- c(0.3, 0.3,  0.45, 0.6)
## aveVElist[[28]] <- c(0.3, 0.45, 0.45, 0.6)
## aveVElist[[29]] <- rep(0, 2)
## aveVElist[[30]] <- c(0.4, 0)
## aveVElist[[31]] <- c(0.4, 0.4)
## aveVElist[[32]] <- rep(0, 3)
## aveVElist[[33]] <- c(0.4, 0,    0)
## aveVElist[[34]] <- c(0.4, 0.4,  0)
## aveVElist[[35]] <- c(0.4, 0.4,  0.4)
## aveVElist[[36]] <- c(0.4, 0,    0,    0)
## aveVElist[[37]] <- c(0.4, 0.4,  0,    0)
## aveVElist[[38]] <- c(0.4, 0.4,  0.4,  0)
## aveVElist[[39]] <- rep(0.4, 4)
## 
## for (i in 1:length(aveVElist)){
##   simTrial(N=c(N.pla, rep(N.vax, length(aveVElist[[i]]))), aveVE=c(0, aveVElist[[i]]), 
##            VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
##            enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
##            visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
##            missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=30, 
##            stage1=78, saveDir="./", randomSeed=300)
##   
##   monitorTrial(dataFile=
##                paste0("simTrial_nPlac=", N.pla, "_nVacc=", 
##                       paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),
##                       "_aveVE=", paste(aveVElist[[i]], collapse="_"), "_infRate=0.04.RData"), 
##                stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=0.0106, 
##                minCnt=50, minPct=0.33, week1=26, minCnt2=2, week2=52, nonEffInterval=20, 
##                nullVE=0, altVE=0.4, highVE=0.6, alpha=0.025, estimand="combined", 
##                VEcutoffWeek=26, saveDir="./")
##   
##   censTrial(dataFile=
##             paste0("simTrial_nPlac=", N.pla, "_nVacc=", 
##                    paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),
##                    "_aveVE=", paste(aveVElist[[i]], collapse="_"), "_infRate=0.04.RData"),
##             monitorFile=
##             paste0("monitorTrial_nPlac=", N.pla, "_nVacc=", 
##                    paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),
##                    "_aveVE=", paste(aveVElist[[i]], collapse="_"), "_infRate=0.04_combined.RData"),
##             stage1=78, stage2=156, saveDir="./")
##   
##   if (i %in% 22:28){
##       rankTrial(censFile=
##                 paste0("trialDataCens_nPlac=", N.pla, "_nVacc=", 
##                        paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),
##                        "_aveVE=", paste(aveVElist[[i]], collapse="_"), "_infRate=0.04_combined.RData"),
##                 idxHighestVE=2, headHead=matrix(c(4,3), nrow=1, ncol=2), 
##                 poolHead=matrix(c(3,4,1,2), nrow=1, ncol=4), stage1=78, stage2=156, 
##                 alpha=0.025, saveDir="./")
##   }
## }
## 
## VEpowerPP(dataList=
##           as.list(paste0("simTrial_nPlac=", N.pla, "_nVacc=", N.vax, "_aveVE=", 
##                          do.call("c", aveVElist[4:11]), "_infRate=0.04.RData")),
##           VEcutoffWeek=26, stage1=78, alpha=0.025, 
##           outName=paste0("VEpwPP_nPlac=", N.pla, "_nVacc=", N.vax, "_infRate=0.04.RData"),
##           saveDir="./")  


###################################################
### code chunk number 3: seqDesignInstructions.Rnw:110-111 (eval = FALSE)
###################################################
## system.file("extdata/seqDesignReport.Rnw", package="seqDesign")


