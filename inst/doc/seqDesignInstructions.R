## ----eval=FALSE, echo=TRUE, tidy=FALSE, size='footnotesize'--------------
#  N.pla <- 1900
#  N.vax <- 1700
#  aveVElist <- list(-2, -1.5, -1, -0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, c(0,0), c(0.4,0),
#                    c(0.4,0.4), c(0.4,0.2), c(0.4,0.3), c(0.5,0.3), c(0.5,0.4), c(0.6,0.3), c(0.6,0.4))
#  infRate <- 0.04
#  estimand <- "cuminc"
#  laggedMonitoring <- TRUE
#  lagTime <- 26
#  outDir <- "./"

## ----eval=FALSE, echo=TRUE, tidy=FALSE, size='footnotesize'--------------
#  for (i in 1:length(aveVElist)){
#    simTrial(N=c(N.pla, rep(N.vax, length(aveVElist[[i]]))), aveVE=c(0, aveVElist[[i]]),
#             VEmodel="half", vePeriods=c(1,27,79), enrollPeriod=78, enrollPartial=13,
#             enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=infRate, fuTime=156,
#             visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#             missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=1000,
#             blockSize=NULL, stage1=78, saveDir=outDir, randomSeed=9)
#  
#    monitorTrial(dataFile=
#                   paste0("simTrial_nPlac=",N.pla,"_nVacc=",
#                          paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),
#                          "_aveVE=",paste(aveVElist[[i]], collapse="_"),"_infRate=",infRate,".RData"),
#                 stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL,
#                 nonEffStartMethod="FKG", nonEffInterval=20, lowerVEnoneff=0, upperVEnoneff=0.4,
#                 stage1VE=0, lowerVEuncPower=0, highVE=0.6, alphaNoneff=0.05, alphaStage1=0.05,
#                 alphaUncPower=0.05, alphaHigh=0.05, estimand=estimand,
#                 laggedMonitoring=laggedMonitoring, lagTime=lagTime, saveDir=outDir)
#  
#    censTrial(dataFile=
#                paste0("simTrial_nPlac=",N.pla, "_nVacc=",
#                       paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),"_aveVE=",
#                       paste(aveVElist[[i]], collapse="_"),"_infRate=",infRate,".RData"),
#              monitorFile=
#                paste0("monitorTrial_nPlac=", N.pla, "_nVacc=",
#                       paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),"_aveVE=",
#                       paste(aveVElist[[i]], collapse="_"),"_infRate=",infRate,"_",estimand,".RData"),
#              stage1=78, stage2=156, saveDir=outDir)
#  
#    if (i %in% 17:22){
#      rankTrial(censFile=
#                  paste0("trialDataCens_nPlac=",N.pla,"_nVacc=",
#                         paste(rep(N.vax, length(aveVElist[[i]])), collapse="_"),"_aveVE=",
#                         paste(aveVElist[[i]], collapse="_"),"_infRate=",infRate,"_",estimand,".RData"),
#                idxHighestVE=1, headHead=matrix(1:2, nrow=1, ncol=2), lowerVE=0, stage1=78, stage2=156,
#                alpha=0.05, saveDir=outDir)
#    }
#  }
#  
#  VEpowerPP(dataList=
#              as.list(paste0("trialDataCens_nPlac=", N.pla, "_nVacc=", N.vax, "_aveVE=",
#                             do.call("c", aveVElist[5:13]), "_infRate=",infRate,"_",estimand,".RData")),
#            lowerVEuncPower=0, alphaUncPower=0.05, VEcutoffWeek=26, stage1=78,
#            outName=paste0("VEpwPP_nPlac=", N.pla, "_nVacc=", N.vax, "_infRate=",infRate,".RData"),
#            saveDir=outDir)

## ----eval=FALSE, echo=TRUE, tidy=FALSE-----------------------------------
#  browseVignettes(package="seqDesign")

