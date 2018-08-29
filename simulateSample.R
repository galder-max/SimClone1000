###############################################
## pipeline to simulate one sample
## for SimClone-1000 simulations
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################
args <- commandArgs(TRUE)
sampleIndex <- as.numeric((args[1])) ## sample index in the grid design (toSimulate.Rda)
######################################################

######################################################
library(simclone) ## library to simulate synthetic tumour data from a tumour architecture
load(file="toSimulate.Rda") ## load grid design (run deriveGrid.R to generate)
######################################################



######################################################
## overwrite simclone functions
## historical: quick fix of an error when mutation fall between two CNA segments
## no longer required with current versions of simclone
simulate_dataset <- function (no.subclones, no.muts, no.subsamples, archive_copynumber_index = NULL,
    segmentation_type = "archive_copynumber", coverage = rep(50,
        no.subsamples), cluster_ccfs = matrix(rep(1, no.subclones *
        no.subsamples), ncol = no.subsamples), sample.cn.param = rep(-1,
        no.subsamples), sample.cn.frac = rep(0, no.subsamples),
    sample.cn.del.frac = rep(0, no.subsamples), sample.subcl.cn.frac = rep(0,
        no.subsamples), cellularity = rep(1, no.subsamples),
    minMutReads = 3, depth.per.cn = NA, copynumber = NULL)
{
    if (length(sample.cn.param) != no.subsamples) {
        stop(paste0("sample.cn.param should have the ", no.subsamples,
            " parameters (i.e. one for every sample)"))
    }
    if (length(cellularity) != no.subsamples) {
        stop(paste0("cellularity should have the ", no.subsamples,
            " parameters (i.e. one for every sample)"))
    }
    if (is.na(depth.per.cn) & length(coverage) != no.subsamples) {
        stop(paste0("coverage should have the ", no.subsamples,
            " parameters (i.e. one for every sample)"))
    }
    if (is.na(depth.per.cn) && length(depth.per.cn) != no.subsamples) {
        stop(paste0("depth.per.cn should have the ", no.subsamples,
            " parameters (i.e. one for every sample)"))
    }
    if (ncol(cluster_ccfs) != no.subsamples | nrow(cluster_ccfs) !=
        no.subclones) {
        stop(paste0("cluster_ccfs should have ", no.subsamples,
            " columns (one per sample) and ", no.subclones, " rows (one for every subclone)"))
    }
    if (!is.null(copynumber) && class(copynumber) != "list") {
        stop(paste0("Supplied copy number should be a list"))
    }
    if (is.null(copynumber)) {
        copynumber = list()
        ploidy = c()
        for (i in 1:no.subsamples) {
            index = length(copynumber) + 1
            copynumber[[index]] = simulate_copynumber(archive_sample_index = archive_copynumber_index,
                use_archive_sample_segments = segmentation_type ==
                  "archive_copynumber", use_chrom_arms_segments = segmentation_type ==
                  "chromosome_arm", use_whole_chrom_segments = segmentation_type ==
                  "whole_chromosome")
            ploidy = c(ploidy, copynumber[[index]]$ploidy)
        }
    }
    else {
        ploidy = unlist(lapply(copynumber, function(x) x$ploidy))
    }
    if (is.na(depth.per.cn)) {
        depth.per.cn = c()
        depth.per.cn.normal = c()
        for (i in 1:no.subsamples) {
            depth.per.cn = c(depth.per.cn, calc_power(cellularity[i],
                ploidy[i], coverage[i]))
            depth.per.cn.normal = c(depth.per.cn.normal, calc_power((1 -
                cellularity[i]), 2, coverage[i]))
        }
    }
    dataset = list()
    for (i in 1:no.subclones) {
        cat(paste0("Simulating subclone", i, "\n"))
        subclone = simulate.subclone(no.muts[i], no.subsamples = no.subsamples,
            depth.per.cn = depth.per.cn, depth.per.cn.normal = depth.per.cn.normal,
            mut.frac.of.cells = cluster_ccfs[i, ], copynumber = copynumber,
            mut.cn.lambda = sample.cn.param, wt.cn.lambda = sample.cn.param,
            sample.cn.frac = sample.cn.frac, sample.cn.del.frac = sample.cn.del.frac,
            sample.subcl.cn.frac = sample.subcl.cn.frac, cellularity = cellularity,
            coverage = coverage)
        dataset[[i]] = mergeColumns(subclone)
        dataset[[i]]$subcloneid = matrix(rep(i, no.muts[i]),
            ncol = 1)
    }
    dataset = appendSubcloneData(dataset)
    dataset$cellularity = cellularity
    dataset$cov = coverage
    dataset$mut.frac.of.cells = cluster_ccfs
    dataset$no.muts = no.muts
    dataset$sample.cn.param = sample.cn.param
    dataset$sample.cn.frac = sample.cn.frac
    dataset$sample.cn.del.frac = sample.cn.del.frac
    dataset$phase = "unphased"
    dataset$copynumber_profiles = copynumber
    selection = sapply(1:nrow(dataset$mutCount), function(i,
        mc, minReads) {
        any(mc[i, ] > minReads)
    }, mc = dataset$mutCount, minReads = minMutReads)
    dataset = subset.dataset(dataset, selection)
    return(dataset)
}

## find closest BB segment for a SNV position
closest_segment <- function(chrom_pos,bb)
{
    distChr <- as.numeric(bb$chr!=chrom_pos[1])*10000000
    distPosStart <- abs(bb[,2]-chrom_pos[2])/100000
    distPosEnd <- abs(bb[,3]-chrom_pos[2])/100000
    dist1 <- distChr+distPosStart
    dist2 <- distChr+distPosStart
    if(min(dist1)<min(dist2)) return(which.min(dist1))
    which.min(dist2)
}

simulate_mutation_copynumber <- function (no.muts, chrom_pos, copynumber, mcn.lambda, real_ccf)
{
    if (real_ccf < 0.95) {
        mult = rep(1, no.muts)
    }
    else {
        mult = rpois(no.muts, mcn.lambda) + 1
    }
    segment_lengths = (copynumber$endpos/1000) - (copynumber$startpos/1000)
    segment_lengths = segment_lengths/sum(segment_lengths)
    total_CN = calculate_bb_total_cn(copynumber)
    nMaj = copynumber$nMaj1_A
    nMin = copynumber$nMin1_A
    assigned_segment = rep(NA, no.muts)
    for (i in 1:no.muts) {
        err <- try(assigned_segment[i] <- which(copynumber$chr == chrom_pos$chrom[i] &
                                               copynumber$startpos < chrom_pos$pos[i]+1 & copynumber$endpos >
                                               chrom_pos$pos[i]),silent=T)
        if(inherits(err,"try-error")) assigned_segment[i] <- closest_segment(chrom_pos[i,],copynumber)## error corrected
        if (nMaj[assigned_segment[i]] < mult[i]) {
            mult[i] = nMaj[assigned_segment[i]]
        }
    }
    output = data.frame(mutation = 1:no.muts, assigned_segment = assigned_segment,
                        mut.cn = mult, wt.cn = total_CN[assigned_segment] - mult)
    return(output)
}

simulate.subclone <- function (no.muts, no.subsamples, depth.per.cn, depth.per.cn.normal,
    mut.frac.of.cells, copynumber, coverage, mut.cn.lambda = rep(-1,
        no.subsamples), wt.cn.lambda = rep(-1, no.subsamples),
    cellularity = rep(1, no.subsamples), sample.cn.frac = rep(0,
        no.subsamples), sample.cn.del.frac = rep(0, no.subsamples),
    sample.subcl.cn.frac = rep(0, no.subsamples), mut.align.bias = 0)
{
    chrom_pos = data.frame(chrom = rep(NA, no.muts), pos = rep(NA,
                                                               no.muts))
    first_sample_copynumber = copynumber[[1]]$simulated_bb
    segment_lengths = (first_sample_copynumber$endpos/1000) -
        (first_sample_copynumber$startpos/1000)
    segment_lengths = segment_lengths/sum(segment_lengths)
    for (i in 1:no.muts) {
        cn_index = sample(1:nrow(first_sample_copynumber), 1,
            prob = segment_lengths)
        chrom_pos$chrom[i] = first_sample_copynumber$chr[cn_index]
        chrom_pos$pos[i] = first_sample_copynumber$startpos[cn_index] +
            round((first_sample_copynumber$endpos[cn_index] -
                first_sample_copynumber$startpos[cn_index]) *
                runif(1))
    }
    subclone = list()
    for (i in 1:no.subsamples) {
        res = simulate_mutation_copynumber(no.muts, chrom_pos,
            copynumber[[i]]$simulated_bb, mut.cn.lambda[i], mut.frac.of.cells)
        mut.cn = res$mut.cn
        wt.cn = res$wt.cn
        mut.subcl.cn = rep(0, no.muts)
        wt.subcl.cn = rep(0, no.muts)
        muts = generate.mut(no.muts = no.muts, depth.per.cn = depth.per.cn[i],
            depth.per.cn.normal = depth.per.cn.normal[i], mut.cn = mut.cn,
            mut.subcl.cn = mut.subcl.cn, wt.cn = wt.cn, wt.subcl.cn = wt.subcl.cn,
            mut.frac.of.cells = as.numeric(mut.frac.of.cells[i]),
            cellularity = cellularity[i], mut.align.bias = mut.align.bias,
            coverage = coverage[i])
        subclone[[i]] = sim.muts2dataset(muts, cellularity[i],
            chrom_pos)
    }
    return(subclone)
}
######################################################



######################################################
## FUNCTIONS TO GENERATE INPUT TO SIMULATION FROM GRID COMBO VALUES
## create diploid BB file from template
diploidBB <- function(sex=c("female","male"),outdir,sampname)
{
    tt <- read.table(paste0("CNA/",
                            "0a6be23a-d5a0-4e95-ada2-a61b2b5d9485",## taken as template
                            ".consensus.20170119.somatic.cna.annotated.txt"),
                     header=T,sep="\t")
    tt <- tt[,c(1:6,17:26)]##extracts battenberg columns from the consensus (required by simclone)
    colnames(tt) <- gsub("battenberg_","",colnames(tt))
    ntt <- tt
    for(chr in c(1:22,"X","Y"))
    {
        wc <- which(ntt$chromosome==chr)
        ntt[wc[1],"start"] <- min(ntt[wc,"start"],na.rm=T)
        ntt[wc[1],"end"] <- max(ntt[wc,"end"],na.rm=T)
        ntt[wc[1],"total_cn"] <- 2
        ntt[wc[1],"major_cn"] <- 1
        ntt[wc[1],"minor_cn"] <- 1
        ntt[wc[1],7:16] <- rep(NA,length(7:16))
        if(length(wc)>1) ntt <- ntt[-c(wc[2:length(wc)]),]
    }
    cln <- colnames(ntt)
    levels(ntt[,"chromosome"])[levels(ntt[,"chromosome"])=="X"] <- "X"
    levels(ntt[,"chromosome"])[levels(ntt[,"chromosome"])=="Y"] <- "Y"
    ntt$chromosome=as.character(ntt$chromosome)
    ntt[,"frac1_A"] <- 1
    colnames(ntt)[1:6] <- c("chr","startpos","endpos","total_cn","nMaj1_A","nMin1_A")
    ntt <- ntt[,!duplicated(colnames(ntt))]
    if(sex[1]=="female")
    {
        ntt<- ntt[ntt[,"chr"]!="Y",]
    }
    if(sex[1]=="male")
    {
        ntt[ntt[,"chr"]=="Y","total_cn"] <- 1
        ntt[ntt[,"chr"]=="X","total_cn"] <- 1
        ntt[ntt[,"chr"]=="Y","nMaj1_A"] <- 1
        ntt[ntt[,"chr"]=="X","nMaj1_A"] <- 1
        ntt[ntt[,"chr"]=="Y","nMin1_A"] <- 0
        ntt[ntt[,"chr"]=="X","nMin1_A"] <- 0
    }
    rownames(ntt) <- NULL
    nttConsensus <- ntt
    colnames(nttConsensus) <- c("chromosome","start","end","total_cn","major_cn","minor_cn")
    write.table(nttConsensus,sep="\t",
                file=paste0(outdir,sampname,"/",sampname,"_consensus_format_subclones.txt"),
                col.names=T,row.names=F,quote=F)
    ntt
}

## get CNA file
getBB <- function(sampname,bbname,outdir)
{
    t <- read.table(paste0("CNA/",
                           bbname,
                           ".consensus.20170119.somatic.cna.annotated.txt"),
                    header=T,sep="\t")
    write.table(t,sep="\t",
                file=paste0(outdir,sampname,"/",sampname,"_consensus_format_subclones.txt"),
                col.names=T,row.names=F,quote=F)
    t <- t[,c(1:6,grep("battenberg",colnames(t)))]
    t$chromosome <- as.character(t$chromosome)
    colnames(t)[1:6] <- c("chr","startpos","endpos","total_cn","nMaj1_A","nMin1_A")
    colnames(t) <- gsub("battenberg_","",colnames(t))
    t <- t[,!duplicated(colnames(t))]
    t[,7:14] <- NA
    t[,7] <- 1
    t <- t[!is.na(t$total_cn),]
    t
}

## write BB file
writeDS_BB<- function(ds,
                      outdir,
                      sampname="test1")
{
    write.table(ds$copynumber_profiles[[1]]$simulated_bb,sep="\t",
                file=paste0(outdir,sampname,"/",sampname,"_subclones.txt"),
                col.names=T,row.names=F,quote=F)
}

## write rho psi file
write_rho_psi <- function(combo,outdir,sampname)
{
    write.table(file=paste0(outdir,sampname,"/",sampname,"_purity_ploidy.txt"),
                data.frame(rho=combo$purity,ploidy=2.0),
                sep="\t",col.names=T,row.names=F,quote=F)
}

## write VCF file
writeDS_VCF <- function(ds,
                      outdir,
                      sampname="test1")
{
    randomise <- sample(1:length(ds$chromosome))
    nds <- data.frame(chromosome=ds$chromosome[randomise],
                      position=ds$position[randomise],
                      refCount=ds$WTCount[randomise],
                      altCount=ds$mutCount[randomise])
    write.table(nds,
                sep="\t",
                file=paste0(outdir,sampname,"/",sampname,".snv.txt"),
                col.names=T,row.names=F,quote=F)
}

## save dataset
saveDS <- function(ds,
                   outdir="/srv/shared/vanloo/home/mtarabichi/simulations/",
                   sampname="test1")
{
    save(ds,file=paste0(outdir,sampname,"/",sampname,".simulated.Rda"))
}

## get sim.tree from combo
sim.tree <- function(combo)
{
    label <- sapply(0:(combo$num_subclones),function(x) paste("M:",x,":",sep="",collapse=":"))
    label <- gsub(":0","",label)
    ancestor <- c("root",rep("M:",combo$num_subclones))
    node <- 1:(combo$num_subclones+1)
    theta.S1 <- c(1,combo$locationCCFs)
    df <- data.frame(label=label,ancestor=ancestor,node=node,theta.S1=theta.S1)
    rownames(df) <- df[,1]
    df
}

## numerise chromosome names
numerise <- function(chr)
{
    chr[chr=="X"] <- 23
    chr[chr=="Y"] <- 24
    as.numeric(chr)
}

## shuffle dataset (keep relationships intact)
shuffleDS <- function(ds)
{
    shuffle <- order(numerise(ds$chromosome),ds$position,decreasing=F)
    for(elem in names(ds)[1:10])
    {
        ds[[elem]] <- matrix(data=ds[[elem]][shuffle],length(shuffle),1,dimnames=list(NULL,elem))
    }
    ds$copynumber_profiles$random_multiplicity <- ds$copynumber_profiles$random_multiplicity[shuffle]
    ds
}

## write Rho file
writeRho <- function(purityOut,outdir,sampname)
{
    fileName <- paste0(outdir,
                       sampname,
                       "/simulated_0001/simulated_0001_0001/",
                       "copynumber/battenberg/simulated_0001_0001_rho_and_psi.txt")
    t <- read.table(fileName,
                    header=T,sep="\t")
    t$rho <- c(NA,purityOut,NA)
    write.table(t,file=fileName,quote=F,col.names=T,row.names=T,sep="\t")
    NULL
}

## infer sex in CNA file
inferSex <- function(bb)
{
    return(if("Y"%in%bb[,1]) "male" else "female")
}
######################################################



######################################################
## output directory
outdir <- "/srv/shared/vanloo/home/mtarabichi/simulations1000/"
######################################################


######################################################
## samplename
sampname <- toSimulate$allCombos[[sampleIndex]]$ID
######################################################
set.seed(12345)
## create directory for sample
system(paste0("mkdir ",outdir,sampname,"/"))
if(!file.exists(paste0(outdir,sampname,"/simulated_0001/")))
{
    ## generate simulated data
    combo <- toSimulate$allCombos[[sampleIndex]]
    tree <- sim.tree(combo)
    num_muts_cluster <- round(c(combo$corrected_num_clonal,combo$correctedSizes))
    num_muts_cluster[num_muts_cluster==0] <- 1
    bbprofile <- getBB(sampname,combo$purity$sampleName,outdir)
    sex <- inferSex(bbprofile)
    bb <- list(simulated_bb = bbprofile,
               random_multiplicity = NA, #rep(1,sum(num_muts_cluster)),
               figures = NA,
               segments = NA,
               raw_data = NA,
               sex = sex,
               ploidy = combo$purity$ploidy)
    ds = simulate_dataset(no.subclones=nrow(tree),
                          no.muts=num_muts_cluster,
                          no.subsamples=1,
                          sample.cn.param=c(1),
                          cluster_ccfs=as.matrix(tree[,grepl("theta", colnames(tree))]),
                          cellularity=c(combo$purity$purityIn),
                          coverage=c(combo$coverageAneuploid),
                          copynumber=list(bb),
                          minMutReads = 2)
    ds <- shuffleDS(ds)
    writeDS_BB(ds,outdir,sampname)
    writeDS_VCF(ds,outdir,sampname)
    saveDS(ds,outdir,sampname)
    write_rho_psi(combo,outdir,sampname)
    writeDataset(outpath=paste0(outdir,sampname,"/"), # path where to write the files
                 datasetname="dataset_g", # general name of the dataset, used for naming
                 datasets=list(ds), # list of datasets
                 directory_archive=list(init_donor(paste0(outdir,sampname,"/"),1,1)),  # list of donors (init_donor output)
                 project_dir=init_donor(paste0(outdir,sampname,"/"),1,1)$donor_dir,  # return value from init_donor
                 seed=12345, # seed used
                 archive_copynumber_index=NA,  # not relevant here
                 segmentation_type=NA,  # not relevant here
                 trees=list(tree),  # list of all trees, must be in the correct format, see the build_tree function that may help
                 num_clusters=NULL) # not required when supplying trees
    writeRho(combo$purity$purityOut,outdir,sampname)
    ## ####################################################
    gc()
    ## ####################################################
    addDiploid <- sampleIndex%in%toSimulate$diploidIndex
    ## ####################################################
    rm(ds)
    gc()
    ## ####################################################
    if(addDiploid)
    {
        set.seed(12345)
        sampname <- toSimulate$allCombos[[sampleIndex]]$IDdiploid
        system(paste0("mkdir ",outdir,sampname,"/"))
        combo <- toSimulate$allCombos[[sampleIndex]]
        tree <- sim.tree(combo)
        num_muts_cluster <- round(c(combo$corrected_num_clonal,combo$correctedSizes))
        bb <- list(simulated_bb = diploidBB(sex,outdir,sampname),
                   random_multiplicity = rep(1,sum(num_muts_cluster)),
                   figures = NA,
                   segments = NA,
                   raw_data = NA,
                   sex = sex,
                   ploidy = 2)
        ds = simulate_dataset(no.subclones=nrow(tree),
                              no.muts=num_muts_cluster,
                              no.subsamples=1,
                              sample.cn.param=c(1),
                              cluster_ccfs=as.matrix(tree[,grepl("theta", colnames(tree))]),
                              cellularity=c(combo$purity$purityIn),
                              coverage=c(combo$coverage),
                              copynumber=list(bb),
                              minMutReads = 2)
        ds <- shuffleDS(ds)
        writeDS_BB(ds,outdir,sampname)
        writeDS_VCF(ds,outdir,sampname)
        saveDS(ds,outdir,sampname)
        write_rho_psi(combo,outdir,sampname)
        writeDataset(outpath=paste0(outdir,sampname,"/"), ## path where to write the files
                     datasetname="dataset_g", ## general name of the dataset, used for naming
                     datasets=list(ds), # list of datasets
                     directory_archive=list(init_donor(paste0(outdir,sampname,"/"),1,1)),  # list of donors (init_donor output)
                     project_dir=init_donor(paste0(outdir,sampname,"/"),1,1)$donor_dir,  # return value from init_donor
                     seed=12345, # seed used
                     archive_copynumber_index=NA,  # not relevant here
                     segmentation_type=NA,  # not relevant here
                     trees=list(tree),  # list of all trees, must be in the correct format, see the build_tree function that may help
                     num_clusters=NULL) # not required when supplying trees
        writeRho(combo$purity$purityOut,outdir,sampname)
        ## ####################################################
    }
}
