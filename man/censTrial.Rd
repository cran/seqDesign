% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqDesign.R
\name{censTrial}
\alias{censTrial}
\title{Generation of Pre-Unblinded Follow-Up Data-Sets by Applying the Monitoring Outcomes}
\usage{
censTrial(dataFile, monitorFile, stage1, stage2, saveFile = NULL,
  saveDir = NULL, verbose = TRUE)
}
\arguments{
\item{dataFile}{if \code{saveDir = NULL}, a list returned by \code{simTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{simTrial}}

\item{monitorFile}{if \code{saveDir = NULL}, a list returned by \code{monitorTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{monitorTrial}}

\item{stage1}{the final week of stage 1 in a two-stage trial}

\item{stage2}{the final week of stage 2 in a two-stage trial, i.e., the maximum follow-up time}

\item{saveFile}{a character string specifying the name of the output \code{.RData} file. If \code{NULL} (default), a default file name will be used.}

\item{saveDir}{a character string specifying a path for both \code{dataFile} and \code{monitorFile}. If supplied, the output is also saved as an \code{.RData} file in this directory; otherwise the output is returned as a list.}

\item{verbose}{a logical value indicating whether information on the output directory and file name should be printed out (default is \code{TRUE})}
}
\value{
If \code{saveDir} is specified, the output list (named \code{trialListCensor}) is saved as an \code{.RData} file in \code{saveDir} (the path to \code{saveDir} is printed); otherwise it is returned. 
The output object is a list of length equal to the number of simulated trials, each of which is a \code{data.frame} with at least the variables \code{trt}, \code{entry}, \code{exit}, and \code{event} 
storing the treatment assignments, enrollment times, correctly censored study exit times, and event indicators, respectively. If available, indicators belonging to the per-protocol cohort 
(named \code{pp1}, \code{pp2}, etc.) are copied from the uncensored data-sets.
}
\description{
\code{censTrial} `correctly censors' treatment arms in data-sets generated by \code{simTrial} by including pre-unblinded follow-up data only according to the monitoring conclusions as reported by \code{monitorTrial}.
}
\details{
All time variables use week as the unit of time. Month is defined as 52/12 weeks.

The following censoring rules are applied to each data-set generated by \code{simTrial}:
\itemize{
  \item If no vaccine arm registers efficacy or high efficacy in Stage 1, the placebo arm is censored on the date when the last vaccine arm hits the harm or non-efficacy boundary.
  \item If a vaccine arm hits the harm boundary, censor the arm immediately.
  \item If a vaccine arm hits the non-efficacy boundary, censor the arm on the earliest date of the two events: (1) the last vaccine arm hits the harm or non-efficacy boundary (if applicable); and (2) all subjects in the vaccine arm have completed the final \code{stage1} visit.
}
}
\examples{
simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
                    VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
                    enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
                    infecRate=0.04, fuTime=156, 
                    visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
                    missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
                    stage1=78, randomSeed=300)

monitorData <- monitorTrial(dataFile=simData, stage1=78, stage2=156, 
                            harmMonitorRange=c(10,100), alphaPerTest=NULL, 
                            nonEffStartMethod="FKG", nonEffInterval=20, 
                            lowerVEnoneff=0, upperVEnoneff=0.4, highVE=0.7, 
                            stage1VE=0, lowerVEuncPower=0, alphaNoneff=0.05, 
                            alphaHigh=0.05, alphaStage1=0.05, 
                            alphaUncPower=0.05, estimand="cuminc", lagTime=26)

censData <- censTrial(dataFile=simData, monitorFile=monitorData, stage1=78, stage2=156)

### alternatively, to save the .RData output file (no '<-' needed):
###
### simTrial(N=c(1400, rep(1000, 2)), aveVE=seq(0, 0.4, by=0.2), VEmodel="half", 
###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=30, 
###          stage1=78, saveDir="./", randomSeed=300)
###
### monitorTrial(dataFile=
###              "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData", 
###              stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL, 
###              nonEffStartMethod="FKG", nonEffInterval=20, lowerVEnoneff=0, 
###              upperVEnoneff=0.4, highVE=0.7, stage1VE=0, lowerVEuncPower=0, 
###              alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, alphaUncPower=0.05, 
###              estimand="cuminc", lagTime=26, saveDir="./")
###
### censTrial(dataFile=
###          "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData",
###          monitorFile=
###          "monitorTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04_cuminc.RData",
###          stage1=78, stage2=156, saveDir="./")
 
}
\seealso{
\code{\link{simTrial}}, \code{\link{monitorTrial}}, and \code{\link{rankTrial}}
}
