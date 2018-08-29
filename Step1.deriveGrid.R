###############################################
## pipeline to generate grid and simulation info
## for SimClone-1000 simulations
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################


###############################################
## file paths
## PATH_PRECOMPUTED_FRACTIONABERRATED <- "allPCAWG.CNA.fAberrated.Rda"
CNAPATH <- "CNA/" ## directory containing consensus cna files (text files)
PCAWG_PURITY_PLOIDY <- "consensus.20170217.purity.ploidy.txt" ## path to summary file with purity and ploidy per sample
PCAWG_DEPTH_COVERAGE <- "all_coverage.txt" ## filename of precomputed average depths of coverage for each sample
PCAWG_CORRECTED_SIZES <- "20170406_consensus_clustering_beta_ccf_correction_batch1_2763_samples.tsv" ## file containing winner's curse-corrected cluster sizes
###############################################

###############################################
## Functions for grid generation and plotting
## give N random CCF positions for N clusters between 10% and 90% CCF from a uniform distribution
positionCCFs <- function(N)
{
    return(sort(runif(N,.1,.9),dec=T))
}

## add CCF locations to combos (the R structure containing all info for simulations)
addLocations <- function(combos)
{
    if(combos$num_subclones==0)
    {
        combos$locations <- NULL
    }
    if(combos$num_subclones>0)
    {
        combos$locationCCFs <- positionCCFs(combos$num_subclones)
    }
    combos
}

## N=number of clusters, fc=fraction clonal, nc=number of clonal mutations
## iteratively breaking the stick to assign number of mutations to clusters
stickBreaking <- function(N,fc,nc)
{
    ntot_subclone <- nc
    if(fc<1) ntot_subclone<- nc/fc*(1-fc)
    if(N==1) return(round(ntot_subclone))
    props <- runif(1,.1,.9)
    Ns <- c(ntot_subclone*props,
            if(N==2) round(ntot_subclone*(1-props)) else stickBreaking(N-1,1,round(ntot_subclone*(1-props))))
    return(round(Ns))
}

## add sizes of clusters to combos using a stick breaking approach.
addSizes <- function(combos)
{
    if(combos$num_subclones==0)
    {
        combos$sizes <- NULL
    }
    if(combos$num_subclones>0)
    {
        combos$sizes <- stickBreaking(combos$num_subclones,combos$frac_clonal, combos$num_clonal)
    }
    combos
}

## sample from the truth values for the grid combinations
sampleVec <- function(grid,index)
{
    sample(grid$values[grid$cluster==index],1)
}

## get fraction of genome aberrated for a given sample (samp=PCAWG ID)
computeFractionAberrated <- function(samp)
{
    tt <- read.table(paste0(CNAPATH,samp,
                            ".consensus.20170119.somatic.cna.annotated.txt"),## version of the cna consensus used
                     header=T,sep="\t")
    sizes <- tt$end-tt$start
    abberated <- !(tt$major_cn==1 & tt$minor_cn==1)
    fraction <- sum(as.numeric(abberated)*sizes,na.rm=T)/1000/sum(sizes/1000,na.rm=T)
    fraction
}

## get output purity for one of the grid values using empirical values and the median absolute deviation on their estimates
outputPurity <- function(grid,index)
{
    sampleIndex <- sample(which(grid$cluster==index & !fAberrated<0.1),1)
    pur <- grid$value[sampleIndex]
    mad <- purs$purity_conf_mad[sampleIndex]
    direction <- sample(c(-1,1),1)
    offset <- runif(1,-mad,mad)
    if(pur+offset*direction<0.05) direction <- 1
    if(pur+offset*direction>1) direction <- -1
    offset <- direction*runif(1,-mad,mad)
    list(sampleName=rownames(purs)[sampleIndex],
         ploidy=purs$ploidy[sampleIndex],
         wgd=purs$wgd_status[sampleIndex],
         wgdUncertain=purs$wgd_uncertain[sampleIndex],
         purityIn=pur,
         purityOut=pur+offset,
         mad=mad,
         offset=offset)
}

## add random sample IDs to combinations (allCombos)
addSampleID <- function(allCombos)
{
    for(i in 1:length(allCombos))
    {
        allCombos[[i]]$ID <- paste0("sim",paste(sample(c(letters,0:9),6,rep=T),collapse="",sep=""))
    }
    allCombos
}

## derive depth of coverage to maintain same nrpcc in diploid and aneuploid
addComputedCoverage <- function(allCombos)
{
    for(i in 1:length(allCombos))
    {
        combo <- allCombos[[i]]
        coverage <- combo$coverage/2*combo$purity$ploidy
        combo$coverageAneuploid <- coverage
        allCombos[[i]] <- combo
    }
    allCombos
}

## derives the average number of mutations that
## would be output after simulating N mutations
## and removing the mutations with less than
## alt_read_min
getNout <- function(N,purity,ccf,coverage,alt_reads_min,rep=10)
{
    N <- round(N)
    Nout <- mean(sapply(1:rep,function(x) sum(rbinom(N,coverage,purity*ccf/2)>alt_reads_min)))
}

## derives number of mutations to input to the simulator to get N called mutations
## after removing mutations with less than alt_reads_min on average
correctNdiploid <- function(N,purity,ccf,coverage,error=N/20,iterMax=20,alt_reads_min=2)
{
    Nlower <- N
    Nupper <- N*2
    coverage <- round(coverage)
    while(getNout(Nupper,purity,ccf,coverage,alt_reads_min)<N)
    {
        Nlower <- Nupper
        Nupper <- 2*Nupper
    }
    Nout <- getNout((Nupper+Nlower)/2,purity,ccf,coverage,alt_reads_min)
    count <- 0
    while(abs(N-Nout)>error & count<iterMax)
    {
        if(Nout>N)
        {
            Nupper <- Nlower+(Nupper-Nlower)/2
        }
        else
        {
            Nlower <- Nlower+(Nupper-Nlower)/2
        }
        Nout <- getNout((Nupper+Nlower)/2,purity,ccf,coverage,alt_reads_min)
        count <- count+1
    }
    return((Nupper+Nlower)/2)
}

## corrects the number of mutations to be input so that
## winner's curse effect would lead to the expected number of called mutations
addNcombo <- function(combo)
{
    combo$corrected_num_clonal <- correctNdiploid(combo$num_clonal,
                                        combo$purity$purityIn,
                                        ccf=1,
                                        combo$coverage,
                                        error=combo$num_clonal/100,
                                        iterMax=20,
                                        alt_reads_min=3)
    if(combo$num_subclones>0)
        combo$correctedSizes <- sapply(1:length(combo$sizes),function(x)
            correctNdiploid(combo$sizes[x],
                            combo$purity$purityIn,
                            ccf=combo$locationCCFs[x],
                            combo$coverage,
                            error=combo$sizes[x]/100,
                            iterMax=20,
                            alt_reads_min=3))
    combo
}

## add number of mutations to combos
addAllNCombos <- function(allCombos)
{
    for(i in 1:length(allCombos))
    {
        allCombos[[i]] <- addNcombo(allCombos[[i]])
    }
    allCombos
}

## computes the mean fraction of samples in random samplings with replacement of the total
computeBS <- function(Ntot,Nsamp,Nrep=1000)
{
    kk <- (sapply(1:Nrep,function(x) sum(1:Ntot%in%sample(1:Ntot,Nsamp,rep=T))/Ntot))
    mean(kk)
}

## returns fraction of samples in random samplings with replacement of the total
computeBS2 <- function(Ntot,Nsamp,Nrep=100)
{
    kk <- (sapply(1:Nrep,function(x) sum(1:Ntot%in%sample(1:Ntot,Nsamp,rep=T))/Ntot))
    kk
}

## computes the median, percentile 2.5% and 97.5% of the fractions
## of samples in random samplings with replacement of the total
computeBS3 <- function(Ntot,Nsamp,Nrep=100)
{
    kk <- (sapply(1:Nrep,function(x) sum(1:Ntot%in%sample(1:Ntot,Nsamp,rep=T))/Ntot))
    list(mu=median(kk),p025=quantile(kk,probs=c(0.025)),p975=quantile(kk,probs=c(.975)))
}

## plotting the PCAWG and grid values
plotHist <- function(purityGrid,
                     pursS,
                     breaks.=seq(0,1,.05),
                     xlab="PCAWG purity",
                     trace=sapply(total,function(x) allCombos[[x]]$purity),jitter=FALSE)
{
    hist(purityGrid$values,
         xlab=xlab,
         col=rgb(0,0,0,0),
         main="",
         cex.axis=1,
         cex.lab=1,
         breaks=breaks.);
    sapply(max(purityGrid$cluster):1,function(x)
        hist(purityGrid$values[purityGrid$cluster%in%1:x],breaks=breaks.,col=acols[x],add=T,xlab="",ylab=""))
    sapply(1:max(trace),function(x) rug(if(jitter) pursS[trace==x]+rnorm(sum(trace==x),sd=.1) else pursS[trace==x],col=acols[x],lwd=.2));
}

## plotting all PCAWG and grid values for the grid parameters
doThePlot <- function()
{
    pdf("simulations.grid.pdf",width=8,height=4)
    par(mfrow=c(2,2))
    par(mar=c(4,4,1,1))
    kn <- plotHist(purityGrid,pursS)
    kn <- plotHist(fcGrid,
                   fcS,
                   breaks.=seq(min(fcGrid$values),max(fcGrid$values),length.out=20),
                   xlab="PCAWG fraction clonal",
                   trace=sapply(total,function(x) allCombos[[x]]$frac_clonal))
    kn <- plotHist(ncGrid,
                   log10(Ncs),
                   breaks.=seq(min(ncGrid$values),max(ncGrid$values),length.out=20),
                   xlab="PCAWG number of clonal SNVs (log10)",
                   trace=sapply(total,function(x) allCombos[[x]]$num_clonal))
    kn <- plotHist(list(values=Ns,cluster=Ns+1),
                   nbSS,
                   breaks.=c(0,.5,1,1.5,2,2.5,3,3.5,4,5,6,7),
                   xlab="PCAWG consensus number of subclones",
                   trace=sapply(total,function(x) allCombos[[x]]$num_subclones)+1,jitter=T)
    dev.off()
}
######################################################


######################################################
## fAberrated <- sapply(rownames(purs),computeFractionAberrated)
## save(fAberrated,file=PATH_PRECOMPUTED_FRACTIONABERRATED)
load(file=PATH_PRECOMPUTED_FRACTIONABERRATED)
######################################################


######################################################
## read in PCAWG info
purs <- read.table(file=PCAWG_PURITY_PLOIDY,header=T) ## read in purity ploidy per sample
covs <- read.table(file=PCAWG_DEPTH_COVERAGE,header=T) ## read in depth of coverage
subs <- read.table(file=PCAWG_CORRECTED_SIZES,header=T) ## read in corrected cluster number of SNVs
######################################################
## process input
allSubs <- tapply(1:nrow(subs),as.factor(subs$aliquot),function(x) subs[x,]) ## combine cluster sizes per sample
Nc <- sapply(allSubs,function(x) x[order(x[,4],decreasing=T)[1],5]) ## number of clonal mutations
fc <- sapply(allSubs,function(x) ## fraction of clonal mutations
{
    if(nrow(x)==1) return(1)
    x <-  x[order(x[,4],decreasing=T),5]
    x[1]/sum(x)
})
Ns <- sapply(allSubs,nrow)-1 ##number of subclones
rownames(purs) <- purs[,1]
rownames(covs) <- covs[,1]
purs <- purs[names(allSubs),]
covs <- covs[names(allSubs),]
######################################################

######################################################
## initiate grid values with kmeans
##---------------------------------------------##
purityGrid <- kmeans(purs$purity,centers=6) ## kmeans on PCAWG purity (6 clusters)
purityGrid$values <- purs$purity
purityGrid$mad <- purs$purity_conf_mad ## annotate purities with their median absolute deviation (for later use)
##---------------------------------------------##
coverageGrid <- kmeans(covs$tum_cov,centers=1) ## one cluster for coverage (all simulations have the PCAWG average depth but not nrpcc)
##---------------------------------------------##
ncGrid <- kmeans(log10(Nc)[Ns>0],centers=4) ## cluster of logged number of mutations in samples with subclones (Ns>0)
ncCenters <- c(sort(ncGrid$centers),5) ## add a 10^5 cluster as kmeans center
for(i in 1:10) ## re-run kmeans with the 10^5 cluster forced until convergence
{
    ncGrid <- kmeans(log10(Nc)[Ns>0],centers=ncCenters,iter.max=1)
    ncCenters <- sort(ncGrid$centers)
    ncCenters[5] <- 5
    print(ncCenters)
}
ncGrid$values <- log10(Nc)[Ns>0]
##---------------------------------------------##
fcGrid <- kmeans(fc[Ns>0],centers=4) ## kmeans on fraction of clonal mutations (4 clusters)
fcCenters <- c(sort(fcGrid$centers),.995) ## add a .995 fraction clonal center to kmeans center
for(i in 1:10) ## re-run kmeans with .995 cluster forced until convergence
{
    fcGrid <- kmeans(fc[Ns>0],centers=fcCenters,iter.max=1)
    fcCenters <- sort(fcGrid$centers)
    fcCenters[5] <- .995
    print(fcCenters)
}
fcGrid$values <- fc[Ns>0]
##---------------------------------------------##
######################################################


######################################################
## inititiate list of combinations of the grid values,
## containing all the grid info and indices
allCombos <- list()
count <- 1
for(pp in 1:length(purityGrid$centers))
{
    for(cov in coverageGrid$centers)
    {
        for(nn in 0:3)
        {
            for(NN in 1:length(ncGrid$centers))
            {
                for(ff in 1:length(fcGrid$centers))
                {
                    allCombos[[count]] <- list(purity=pp,##purity
                                               num_subclones=nn,##number of subclones
                                               num_clonal=NN,##number of clonal SNVs
                                               frac_clonal=ff,##fraction clonal SNVs
                                               coverage=cov)##coverage
                    count <- count+1
                    ##increment index
                }
            }
        }
    }
}
count <- count-1
######################################################

######################################################
## add "noise" in the simulations to better hide the design
## from the methods to avoid overfitting
## 1) set1=sample 600 grid combinations with replacement
## 2) set2=sample 100 grid combinations with replacement from the remaining combinations (not in set1)
## below: double-check sampling strategy covers most of the grid
######################################################
require(parallel)
coveragePlus100. <- mclapply(0:600,function(N)
{
    bs1 <- computeBS2(600,N)
    first <- N*bs1
    remaining <- 600-first
    bs2 <- lapply(remaining,function(x) computeBS3(x,600-N+100))
    second <- remaining*sapply(bs2,function(x) x$mu)
    list(mu=mean((first+second)/600),
         p025=mean((first+remaining*sapply(bs2,function(x) x$p025))/600),
         p975=mean((first+remaining*sapply(bs2,function(x) x$p975))/600))
},mc.cores=8)
######################################################
quartz(width=4,height=4)
par(mar=c(4,4,1,1))
plot(0:600,
     coveragePlus100,
     frame=F,
     pch=19,
     col=rgb(.5,.5,.5,1),
     ylim=c(.62,.78),
     xlab="number of samples in first batch",
     ylab="total grid coverage")
points(0:600,
     coveragePlus100,type="l",lwd=2)
polygon(c(0:600,600:0),
        border=NA,
        c(sapply(coveragePlus100.,function(x) x$p025),
          (sapply(coveragePlus100.,function(x) x$p975))[601:1]),
        col=rgb(.4,.3,.8,.2))
######################################################

######################################################
## implement sampling strategy
set.seed(12345)
index <- 1:count
first <- sample(index,rep=T)
second <- sample(index[!index%in%first],100,rep=T)
total <- c(first,second) ##indices of grid combinations to keep
######################################################

######################################################
## Sampled combinations
## and their simulation values sampled from the real values
allCombos2 <- lapply(total,function(x)
{
    list(purity=outputPurity(purityGrid,allCombos[[x]]$purity),##output purity value
         num_subclones=if(allCombos[[x]]$num_subclones==3) sample(3:7,1)
                       else allCombos[[x]]$num_subclones, ## number of subclones from {1,2,[3-7]} where [3-7]~U{3,4,5,6,7}
         num_clonal=round(10^sampleVec(ncGrid,allCombos[[x]]$num_clonal)),##number of clonal mutations
         frac_clonal=sampleVec(fcGrid,allCombos[[x]]$frac_clonal),##fraction clonal
         coverage=cov)## coverage
})
######################################################

######################################################
allCombos3 <- lapply(allCombos2,function(x) addSizes(addLocations(x))) ##add cluster sizes and CCFs
######################################################
diploidIndex <- sample(1:length(allCombos3),1000-length(allCombos3),rep=F) ##sample 300 (1000-700) indices to simulate as fully diploid as well
######################################################

######################################################
## final list containing all necessary info for simulations with Simclone
toSimulate <- list(allCombos=allCombos3, ## all combination info
                   grid=allCombos,## all indices
                   purityGrid=purityGrid, ## purity grid
                   ncGrid=ncGrid, ## number of clonal SNVs grid
                   fcGrid=fcGrid,##fraction clonal grid
                   totalSample=total, ## total indices
                   diploidIndex=diploidIndex) ## diploid indices
######################################################

######################################################
## corrected number of SNVs for winner's curse effect and correct coverage of diploid for different ploidy
toSimulate$allCombos <- addComputedCoverage(addSampleID(toSimulate$allCombos)) ## add diploid coverage  and diploid sampleID
toSimulate$allCombos <- addAllNCombos(toSimulate$allCombos) ## add corrected cluster sizes
######################################################
save(toSimulate,file="toSimulate.Rda") ## save grid design
######################################################



######################################################
## load design and save objects for plotting
load(file="toSimulate.Rda")
######################################################
## extract grid parameters info
pursS <- sapply(allCombos2,function(x) x$purity$purityIn)
nbSS <- sapply(allCombos2,function(x) x$num_subclones)
Ncs <- sapply(allCombos2,function(x) x$num_clonal)
fcS <- sapply(allCombos2,function(x) x$frac_clonal)
######################################################
## list object saved for plotting
imageL <- list(doThePlot=doThePlot,
               plotHist=plotHist,
               purityGrid=purityGrid,
               pursS=pursS,
               fcGrid=fcGrid,
               fcS=fcS,
               total=total,
               ncGrid=ncGrid,
               Ncs=Ncs,
               Ns=Ns,
               acols=acols,
               nbSS=nbSS,
               allCombos=allCombos)
save(imageL,file="imageL.Rda")
######################################################
